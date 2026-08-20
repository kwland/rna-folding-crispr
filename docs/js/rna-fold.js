/*
 * RNA secondary-structure prediction: Nussinov and Zuker.
 *
 * Two models are implemented side by side because they optimize different
 * objectives:
 *
 *   Nussinov  maximizes the COUNT of non-crossing base pairs. Every allowed
 *             pair is worth +1. Transparent, but not physical.
 *
 *   Zuker     minimizes FREE ENERGY (kcal/mol) using nearest-neighbour
 *             thermodynamics: stacking stabilizes, loops cost.
 *
 * PARAMETER PROVENANCE. Read this before quoting numbers in a report:
 *
 *   - Watson-Crick stacking energies are the published Turner/Xia nearest-
 *     neighbour values (Xia et al. 1998, Turner & Mathews 2010), in kcal/mol
 *     at 37 C. These 10 values are exact.
 *
 *   - G-U wobble stacks, loop initiation tables, and the multiloop model are
 *     SIMPLIFIED approximations of the Turner rules. They reproduce the right
 *     qualitative ordering but are not the exact published tables. Each is
 *     flagged with an APPROXIMATION comment below.
 *
 *   Consequence: predicted structures are usually reasonable, but the absolute
 *   kcal/mol will not match ViennaRNA exactly. ViennaRNA/RNAfold remains the
 *   reference implementation for any quantitative claim.
 */

window.RNA = (function () {
  "use strict";

  var CANONICAL = { AU: 1, UA: 1, GC: 1, CG: 1 };
  var WOBBLE = { GU: 1, UG: 1 };

  var MIN_LOOP = 3; // minimum unpaired bases enclosed by a hairpin
  var MAX_LEN = 300; // guard: both algorithms are polynomial, keep the UI responsive
  var MAX_INTERNAL = 30; // Turner convention: internal loops larger than this are ignored
  var EPS = 1e-7; // float tolerance for traceback equality tests

  // ---------------------------------------------------------------- sequence

  function normalize(sequence) {
    var seq = String(sequence || "")
      .replace(/\s+/g, "")
      .toUpperCase()
      .replace(/T/g, "U");
    var invalid = [];
    for (var i = 0; i < seq.length; i++) {
      var ch = seq[i];
      if (ch !== "A" && ch !== "C" && ch !== "G" && ch !== "U" && invalid.indexOf(ch) === -1) {
        invalid.push(ch);
      }
    }
    if (invalid.length) {
      throw new Error("Only A, C, G, U (and T, read as U) are supported. Found: " + invalid.join(", "));
    }
    if (seq.length > MAX_LEN) {
      throw new Error("Sequence is " + seq.length + " nt. This browser demo folds up to " + MAX_LEN + " nt.");
    }
    return seq;
  }

  function canPair(a, b, allowWobble) {
    var key = a + b;
    return CANONICAL[key] === 1 || (allowWobble !== false && WOBBLE[key] === 1);
  }

  function isWobble(a, b) {
    return WOBBLE[a + b] === 1;
  }

  function dotBracket(n, pairs) {
    var chars = new Array(n);
    for (var i = 0; i < n; i++) chars[i] = ".";
    for (var p = 0; p < pairs.length; p++) {
      chars[pairs[p][0]] = "(";
      chars[pairs[p][1]] = ")";
    }
    return chars.join("");
  }

  function matrix(n, fill) {
    var m = new Array(n);
    for (var i = 0; i < n; i++) {
      m[i] = new Array(n);
      for (var j = 0; j < n; j++) m[i][j] = fill;
    }
    return m;
  }

  // ---------------------------------------------------------------- Nussinov

  /*
   * N[i,j] = max( N[i+1,j],
   *               N[i,j-1],
   *               N[i+1,j-1] + pair(i,j),
   *               max_k N[i,k] + N[k+1,j] )
   */
  function nussinov(sequence, minLoop, allowWobble) {
    var seq = normalize(sequence);
    var n = seq.length;
    if (minLoop == null) minLoop = MIN_LOOP;
    if (allowWobble == null) allowWobble = true;
    if (n === 0) return { sequence: "", score: 0, structure: "", pairs: [], dp: [] };

    var dp = matrix(n, 0);

    for (var span = 1; span < n; span++) {
      for (var i = 0; i + span < n; i++) {
        var j = i + span;
        var best = Math.max(dp[i + 1][j], dp[i][j - 1]);

        if (j - i > minLoop && canPair(seq[i], seq[j], allowWobble)) {
          best = Math.max(best, dp[i + 1][j - 1] + 1);
        }
        for (var k = i; k < j; k++) {
          var split = dp[i][k] + dp[k + 1][j];
          if (split > best) best = split;
        }
        dp[i][j] = best;
      }
    }

    var pairs = [];
    var stack = [[0, n - 1]];
    while (stack.length) {
      var frame = stack.pop();
      var a = frame[0];
      var b = frame[1];
      if (a >= b) continue;

      if (b - a > minLoop && canPair(seq[a], seq[b], allowWobble) && dp[a][b] === dp[a + 1][b - 1] + 1) {
        pairs.push([a, b]);
        stack.push([a + 1, b - 1]);
        continue;
      }
      if (dp[a][b] === dp[a + 1][b]) {
        stack.push([a + 1, b]);
        continue;
      }
      if (dp[a][b] === dp[a][b - 1]) {
        stack.push([a, b - 1]);
        continue;
      }
      for (var s = a; s < b; s++) {
        if (dp[a][b] === dp[a][s] + dp[s + 1][b]) {
          stack.push([a, s]);
          stack.push([s + 1, b]);
          break;
        }
      }
    }

    pairs.sort(function (x, y) {
      return x[0] - y[0];
    });
    return { sequence: seq, score: dp[0][n - 1], structure: dotBracket(n, pairs), pairs: pairs, dp: dp };
  }

  // ------------------------------------------------------- Zuker energy model

  var RT37 = 0.6163; // kcal/mol at 37 C
  var LOOP_SCALE = 1.75 * RT37; // Jacobson-Stockmayer extrapolation coefficient

  /*
   * Watson-Crick nearest-neighbour stacking, kcal/mol at 37 C (Xia et al. 1998).
   * Key "WX/ZY" means the duplex 5'-W X-3' paired with 3'-Z Y-5', i.e. the pair
   * (W,Z) stacked directly on the pair (X,Y). These 10 values are exact.
   */
  var WC_STACKS = [
    ["AA/UU", -0.93],
    ["AU/UA", -1.1],
    ["UA/AU", -1.33],
    ["CU/GA", -2.08],
    ["CA/GU", -2.11],
    ["GU/CA", -2.24],
    ["GA/CU", -2.35],
    ["CG/GC", -2.36],
    ["GG/CC", -3.26],
    ["GC/CG", -3.42]
  ];

  /*
   * APPROXIMATION. Real Turner wobble parameters are context-dependent (tandem
   * G-U in particular has special-case tables). They are collapsed here to two flat
   * values that preserve the ordering: wobble stacks are weaker than any
   * Watson-Crick stack, and tandem wobbles are weaker still.
   */
  var WOBBLE_SINGLE_STACK = -1.3;
  var WOBBLE_TANDEM_STACK = -0.5;

  // Terminal A-U / G-U penalty applied at the end of a helix (Turner: +0.45).
  var TERMINAL_AU = 0.45;
  // APPROXIMATION: Turner applies a distinct A-U/G-U closure cost inside internal loops.
  var INTERNAL_TERMINAL_AU = 0.7;

  /*
   * Linear multiloop model: a + b*branches + c*unpaired.
   *
   * b is NEGATIVE, so each extra branch is rewarded, not punished. That looks
   * wrong but is the Turner 2004 convention, and it is load-bearing: with a
   * positive b the model nests tRNA's arms into a chain instead of opening the
   * four-way junction. At b = -0.9 the tRNA-Phe fold recovers the accepted
   * cloverleaf (stems of 7/4/5/5, ~86% positional agreement). See the tRNA-Phe
   * example on the site.
   *
   * APPROXIMATION: the real Turner multiloop term is not strictly linear.
   */
  var ML_CLOSE = 3.4;
  var ML_BRANCH = -0.9;
  var ML_UNPAIRED = 0.0;

  // Loop initiation tables (Turner-style, kcal/mol). Sizes beyond the table are
  // extrapolated logarithmically.
  var HAIRPIN_INIT = { 3: 5.4, 4: 5.6, 5: 5.7, 6: 5.4, 7: 6.0, 8: 5.5, 9: 6.4 };
  var BULGE_INIT = { 1: 3.8, 2: 2.8, 3: 3.2, 4: 3.6, 5: 4.0, 6: 4.4 };
  // APPROXIMATION: sizes 2-3 (1x1, 1x2 loops) use special tables in Turner; flattened here.
  var INTERNAL_INIT = { 2: 1.5, 3: 1.6, 4: 1.7, 5: 1.8, 6: 2.0 };

  var STACK = {};
  (function buildStackTable() {
    for (var s = 0; s < WC_STACKS.length; s++) {
      var parts = WC_STACKS[s][0].split("/");
      var value = WC_STACKS[s][1];
      var W = parts[0][0];
      var X = parts[0][1];
      var Z = parts[1][0];
      var Y = parts[1][1];
      // Key layout: outer pair (W,Z) then inner pair (X,Y) -> "W X Z Y".
      STACK[W + X + Z + Y] = value;
      // A duplex read from the opposite strand is the same duplex: WX/ZY == YZ/XW.
      STACK[Y + Z + X + W] = value;
    }
  })();

  function terminalPenalty(a, b) {
    return CANONICAL[a + b] === 1 && (a === "G" || a === "C") ? 0 : TERMINAL_AU;
  }

  function internalTerminalPenalty(a, b) {
    return CANONICAL[a + b] === 1 && (a === "G" || a === "C") ? 0 : INTERNAL_TERMINAL_AU;
  }

  function stackEnergy(a1, b1, a2, b2) {
    var wc = STACK[a1 + a2 + b1 + b2];
    if (wc !== undefined) return wc;
    var wobbles = (isWobble(a1, b1) ? 1 : 0) + (isWobble(a2, b2) ? 1 : 0);
    if (wobbles >= 2) return WOBBLE_TANDEM_STACK;
    if (wobbles === 1) return WOBBLE_SINGLE_STACK;
    return 0;
  }

  function extrapolate(table, size, anchor) {
    if (table[size] !== undefined) return table[size];
    return table[anchor] + LOOP_SCALE * Math.log(size / anchor);
  }

  function hairpinInit(size) {
    if (size < 3) return Infinity;
    return extrapolate(HAIRPIN_INIT, size, 9);
  }

  function bulgeInit(size) {
    return extrapolate(BULGE_INIT, size, 6);
  }

  function internalInit(size) {
    if (size < 2) return Infinity;
    return extrapolate(INTERNAL_INIT, size, 6);
  }

  /*
   * Zuker MFE folding.
   *
   *   V[i,j]  = MFE of i..j GIVEN that i pairs with j
   *   W[i,j]  = MFE of i..j with no constraint (exterior context)
   *   WM[i,j] = MFE of i..j inside a multiloop, at least one branch
   */
  function zuker(sequence, allowWobble) {
    var seq = normalize(sequence);
    var n = seq.length;
    if (allowWobble == null) allowWobble = true;
    if (n === 0) {
      return { sequence: "", energy: 0, structure: "", pairs: [], V: [], W: [] };
    }

    var V = matrix(n, Infinity);
    var W = matrix(n, 0);
    var WM = matrix(n, Infinity);

    function hairpin(i, j) {
      var size = j - i - 1;
      if (size < MIN_LOOP) return Infinity;
      return hairpinInit(size) + terminalPenalty(seq[i], seq[j]);
    }

    // Energy of the loop closed by (i,j) with inner pair (k,l): stack, bulge, or internal.
    function loopEnergy(i, j, k, l) {
      var l1 = k - i - 1;
      var l2 = j - l - 1;
      var total = l1 + l2;

      if (total === 0) {
        return stackEnergy(seq[i], seq[j], seq[k], seq[l]);
      }
      if (l1 === 0 || l2 === 0) {
        var e = bulgeInit(total);
        if (total === 1) {
          // A single bulged base does not break the helix: the stack survives.
          e += stackEnergy(seq[i], seq[j], seq[k], seq[l]);
        } else {
          e += terminalPenalty(seq[i], seq[j]) + terminalPenalty(seq[k], seq[l]);
        }
        return e;
      }
      var ie = internalInit(total);
      // APPROXIMATION: Turner's asymmetry term uses a fitted coefficient and cap.
      ie += Math.min(3.0, 0.6 * Math.abs(l1 - l2));
      ie += internalTerminalPenalty(seq[i], seq[j]) + internalTerminalPenalty(seq[k], seq[l]);
      return ie;
    }

    for (var span = 1; span < n; span++) {
      for (var i = 0; i + span < n; i++) {
        var j = i + span;

        // ---- V
        var v = Infinity;
        if (j - i - 1 >= MIN_LOOP && canPair(seq[i], seq[j], allowWobble)) {
          v = hairpin(i, j);

          for (var k = i + 1; k <= j - MIN_LOOP - 2; k++) {
            var l1 = k - i - 1;
            if (l1 > MAX_INTERNAL) break;
            for (var l = j - 1; l > k + MIN_LOOP; l--) {
              var l2 = j - l - 1;
              if (l1 + l2 > MAX_INTERNAL) break;
              if (V[k][l] === Infinity) continue;
              var cand = loopEnergy(i, j, k, l) + V[k][l];
              if (cand < v) v = cand;
            }
          }

          for (var m = i + 2; m < j - 1; m++) {
            if (WM[i + 1][m] === Infinity || WM[m + 1][j - 1] === Infinity) continue;
            var ml = WM[i + 1][m] + WM[m + 1][j - 1] + ML_CLOSE + ML_BRANCH + terminalPenalty(seq[i], seq[j]);
            if (ml < v) v = ml;
          }
        }
        V[i][j] = v;

        // ---- WM
        var wm = Infinity;
        if (V[i][j] !== Infinity) {
          wm = V[i][j] + ML_BRANCH + terminalPenalty(seq[i], seq[j]);
        }
        if (WM[i + 1][j] !== Infinity && WM[i + 1][j] + ML_UNPAIRED < wm) wm = WM[i + 1][j] + ML_UNPAIRED;
        if (WM[i][j - 1] !== Infinity && WM[i][j - 1] + ML_UNPAIRED < wm) wm = WM[i][j - 1] + ML_UNPAIRED;
        for (var q = i; q < j; q++) {
          if (WM[i][q] === Infinity || WM[q + 1][j] === Infinity) continue;
          var mm = WM[i][q] + WM[q + 1][j];
          if (mm < wm) wm = mm;
        }
        WM[i][j] = wm;

        // ---- W (exterior loop: unpaired bases are free)
        var w = 0;
        if (W[i + 1][j] < w) w = W[i + 1][j];
        if (W[i][j - 1] < w) w = W[i][j - 1];
        if (V[i][j] !== Infinity) {
          var closed = V[i][j] + terminalPenalty(seq[i], seq[j]);
          if (closed < w) w = closed;
        }
        for (var r = i; r < j; r++) {
          var sp = W[i][r] + W[r + 1][j];
          if (sp < w) w = sp;
        }
        W[i][j] = w;
      }
    }

    // ---- traceback
    var pairs = [];

    function eq(a, b) {
      if (a === Infinity || b === Infinity) return false;
      return Math.abs(a - b) < EPS;
    }

    function traceW(i, j) {
      if (j - i < MIN_LOOP + 1) return;
      if (eq(W[i][j], W[i + 1][j])) return traceW(i + 1, j);
      if (eq(W[i][j], W[i][j - 1])) return traceW(i, j - 1);
      if (V[i][j] !== Infinity && eq(W[i][j], V[i][j] + terminalPenalty(seq[i], seq[j]))) {
        pairs.push([i, j]);
        return traceV(i, j);
      }
      for (var k = i; k < j; k++) {
        if (eq(W[i][j], W[i][k] + W[k + 1][j])) {
          traceW(i, k);
          traceW(k + 1, j);
          return;
        }
      }
    }

    function traceV(i, j) {
      var v = V[i][j];
      if (v === Infinity) return;
      if (eq(v, hairpin(i, j))) return;

      for (var k = i + 1; k <= j - MIN_LOOP - 2; k++) {
        var l1 = k - i - 1;
        if (l1 > MAX_INTERNAL) break;
        for (var l = j - 1; l > k + MIN_LOOP; l--) {
          var l2 = j - l - 1;
          if (l1 + l2 > MAX_INTERNAL) break;
          if (V[k][l] === Infinity) continue;
          if (eq(v, loopEnergy(i, j, k, l) + V[k][l])) {
            pairs.push([k, l]);
            return traceV(k, l);
          }
        }
      }

      for (var m = i + 2; m < j - 1; m++) {
        if (WM[i + 1][m] === Infinity || WM[m + 1][j - 1] === Infinity) continue;
        if (eq(v, WM[i + 1][m] + WM[m + 1][j - 1] + ML_CLOSE + ML_BRANCH + terminalPenalty(seq[i], seq[j]))) {
          traceWM(i + 1, m);
          traceWM(m + 1, j - 1);
          return;
        }
      }
    }

    function traceWM(i, j) {
      if (i >= j) return;
      var wm = WM[i][j];
      if (wm === Infinity) return;
      if (V[i][j] !== Infinity && eq(wm, V[i][j] + ML_BRANCH + terminalPenalty(seq[i], seq[j]))) {
        pairs.push([i, j]);
        return traceV(i, j);
      }
      if (WM[i + 1][j] !== Infinity && eq(wm, WM[i + 1][j] + ML_UNPAIRED)) return traceWM(i + 1, j);
      if (WM[i][j - 1] !== Infinity && eq(wm, WM[i][j - 1] + ML_UNPAIRED)) return traceWM(i, j - 1);
      for (var k = i; k < j; k++) {
        if (WM[i][k] === Infinity || WM[k + 1][j] === Infinity) continue;
        if (eq(wm, WM[i][k] + WM[k + 1][j])) {
          traceWM(i, k);
          traceWM(k + 1, j);
          return;
        }
      }
    }

    traceW(0, n - 1);
    pairs.sort(function (a, b) {
      return a[0] - b[0];
    });

    return {
      sequence: seq,
      energy: W[0][n - 1],
      structure: dotBracket(n, pairs),
      pairs: pairs,
      V: V,
      W: W
    };
  }

  // ------------------------------------------- McCaskill partition function

  /*
   * Nussinov and Zuker each return ONE structure, so every base is scored a
   * hard paired or unpaired. Real RNA moves through a Boltzmann-weighted
   * ensemble of structures. McCaskill (1990) computes, in the same O(n^3)
   * dynamic-programming family as Zuker, the probability that each pair (i, j)
   * forms across that whole ensemble. The probability that base i is unpaired
   * is then 1 minus the sum of its pairing probabilities.
   *
   * Inside pass, over increasing span:
   *   QB[i][j]  partition function of i..j given that i pairs with j
   *   QM1[i][j] one multiloop branch starting exactly at i, then unpaired to j
   *   QM[i][j]  at least one multiloop branch somewhere in i..j
   *   Qf[t]     exterior partition function of the prefix 0..t-1
   *   Qr[t]     exterior partition function of the suffix t..n-1
   *
   * Outside pass, over decreasing span. A pair forms in exactly one of three
   * contexts, so its probability is the sum of three terms: nothing encloses
   * it, one closer pair encloses it alone, or it is one branch of a multiloop.
   * That third term is O(n^4) written directly; the two accumulator tables
   * below carry the enclosing sum forward instead, keeping the pass at O(n^3).
   *
   * This is a port of mccaskill.py. The Python version is checked against
   * exhaustive enumeration of every structure for short sequences, and the two
   * implementations are checked against each other.
   */

  var RT37_PF = RT37; // Boltzmann temperature factor, kcal/mol at 37 C

  function boltz(energy) {
    return energy === Infinity ? 0 : Math.exp(-energy / RT37_PF);
  }

  function mccaskill(sequence, allowWobble) {
    var seq = normalize(sequence);
    var n = seq.length;
    if (allowWobble == null) allowWobble = true;
    if (n === 0) {
      return { sequence: "", partitionFunction: 1, ensembleEnergy: 0, pairProb: [], unpaired: [] };
    }

    var i, j, k, l, h, span;

    function pairable(a, b) {
      return b - a > MIN_LOOP && canPair(seq[a], seq[b], allowWobble);
    }

    // Multiloop weights. ML_UNPAIRED is 0 in this parameter set, so zPow is all
    // ones, but the factors are kept so a changed parameter stays correct.
    var zUnpaired = Math.exp(-ML_UNPAIRED / RT37_PF);
    var zPow = new Array(n + 2);
    for (i = 0; i <= n + 1; i++) zPow[i] = Math.pow(zUnpaired, i);
    var mlOpen = Math.exp(-(ML_CLOSE + ML_BRANCH) / RT37_PF);
    var mlBranch = Math.exp(-ML_BRANCH / RT37_PF);

    // Loop weights tabulated once. Exponentials in the innermost loop dominated
    // the runtime before this.
    var term = matrix(n, 1);
    var intTerm = matrix(n, 1);
    var wHairpin = matrix(n, 0);
    var wStack = matrix(n, 0);
    var wBulge1L = matrix(n, 0);
    var wBulge1R = matrix(n, 0);
    for (i = 0; i < n; i++) {
      for (j = i + MIN_LOOP + 1; j < n; j++) {
        if (!pairable(i, j)) continue;
        term[i][j] = boltz(terminalPenalty(seq[i], seq[j]));
        intTerm[i][j] = boltz(internalTerminalPenalty(seq[i], seq[j]));
        wHairpin[i][j] = boltz(hairpinInit(j - i - 1) + terminalPenalty(seq[i], seq[j]));
        if (pairable(i + 1, j - 1)) {
          wStack[i][j] = boltz(stackEnergy(seq[i], seq[j], seq[i + 1], seq[j - 1]));
        }
        if (i + 2 < j - 1 && pairable(i + 2, j - 1)) {
          wBulge1L[i][j] = boltz(bulgeInit(1) + stackEnergy(seq[i], seq[j], seq[i + 2], seq[j - 1]));
        }
        if (i + 1 < j - 2 && pairable(i + 1, j - 2)) {
          wBulge1R[i][j] = boltz(bulgeInit(1) + stackEnergy(seq[i], seq[j], seq[i + 1], seq[j - 2]));
        }
      }
    }
    var wBulge = new Array(MAX_INTERNAL + 2);
    for (i = 2; i < MAX_INTERNAL + 2; i++) wBulge[i] = boltz(bulgeInit(i));
    var wInternal = [];
    for (i = 0; i <= MAX_INTERNAL; i++) wInternal.push(new Array(MAX_INTERNAL + 1).fill(0));
    for (var left = 1; left <= MAX_INTERNAL; left++) {
      for (var right = 1; right <= MAX_INTERNAL; right++) {
        var total = left + right;
        if (total > MAX_INTERNAL || total < 2) continue;
        wInternal[left][right] = boltz(internalInit(total) + Math.min(3.0, 0.6 * Math.abs(left - right)));
      }
    }

    // ---------------------------------------------------------- inside pass
    var QB = matrix(n, 0);
    var QM = matrix(n, 0);
    var QM1 = matrix(n, 0);

    for (span = MIN_LOOP + 1; span < n; span++) {
      for (i = 0; i + span < n; i++) {
        j = i + span;

        if (pairable(i, j)) {
          var sum = wHairpin[i][j];
          var kMax = Math.min(i + MAX_INTERNAL + 1, j - MIN_LOOP - 2);
          for (k = i + 1; k <= kMax; k++) {
            var l1 = k - i - 1;
            var lMin = Math.max(k + MIN_LOOP + 1, j - 1 - (MAX_INTERNAL - l1));
            for (l = j - 1; l >= lMin; l--) {
              var inner = QB[k][l];
              if (inner === 0) continue;
              var l2 = j - l - 1;
              if (l1 === 0) {
                if (l2 === 0) sum += inner * wStack[i][j];
                else if (l2 === 1) sum += inner * wBulge1R[i][j];
                else sum += inner * wBulge[l2] * term[i][j] * term[k][l];
              } else if (l2 === 0) {
                if (l1 === 1) sum += inner * wBulge1L[i][j];
                else sum += inner * wBulge[l1] * term[i][j] * term[k][l];
              } else {
                sum += inner * wInternal[l1][l2] * intTerm[i][j] * intTerm[k][l];
              }
            }
          }
          // Multiloop: at least two branches inside, the last one opening at h.
          if (span >= 2 * (MIN_LOOP + 2) + 1) {
            var ml = 0;
            for (h = i + 2; h < j - MIN_LOOP - 1; h++) {
              var leftBranches = QM[i + 1][h - 1];
              if (leftBranches === 0) continue;
              var rightBranch = QM1[h][j - 1];
              if (rightBranch === 0) continue;
              ml += leftBranches * rightBranch;
            }
            if (ml) sum += ml * mlOpen * term[i][j];
          }
          QB[i][j] = sum;
        }

        var q1 = j - 1 >= i ? QM1[i][j - 1] * zUnpaired : 0;
        if (QB[i][j]) q1 += QB[i][j] * mlBranch * term[i][j];
        QM1[i][j] = q1;

        var acc = 0;
        for (h = i; h < j - MIN_LOOP; h++) {
          var branch = QM1[h][j];
          if (branch === 0) continue;
          acc += (zPow[h - i] + (h - 1 >= i ? QM[i][h - 1] : 0)) * branch;
        }
        QM[i][j] = acc;
      }
    }

    // Exterior loop, as a prefix and a suffix partition function.
    var Qf = new Array(n + 1).fill(1);
    for (var t = 1; t <= n; t++) {
      var jj = t - 1;
      var vf = Qf[t - 1];
      for (h = 0; h < jj - MIN_LOOP; h++) {
        if (QB[h][jj]) vf += Qf[h] * QB[h][jj] * term[h][jj];
      }
      Qf[t] = vf;
    }
    var Qr = new Array(n + 2).fill(1);
    for (t = n - 1; t >= 0; t--) {
      var vr = Qr[t + 1];
      for (l = t + MIN_LOOP + 1; l < n; l++) {
        if (QB[t][l]) vr += QB[t][l] * term[t][l] * Qr[l + 1];
      }
      Qr[t] = vr;
    }

    var Q = Qf[n];
    if (!isFinite(Q) || Q <= 0) {
      throw new Error("The partition function overflowed. Fold a shorter sequence.");
    }

    // --------------------------------------------------------- outside pass
    var P = matrix(n, 0);
    var accUnpaired = matrix(n, 0); // enclosing pairs with everything right of q unpaired
    var accBranch = matrix(n, 0); // enclosing pairs with a further branch right of q
    var pending = [];

    for (span = n - 1; span > MIN_LOOP; span--) {
      for (var p = 0; p < pending.length; p++) {
        var ek = pending[p][0];
        var el = pending[p][1];
        var weight = pending[p][2];
        for (var q = ek + 1; q < el; q++) {
          accUnpaired[ek][q] += weight * zPow[el - 1 - q];
          if (q + 1 <= el - 1 && QM[q + 1][el - 1]) {
            accBranch[ek][q] += weight * QM[q + 1][el - 1];
          }
        }
      }
      pending = [];

      for (i = 0; i + span < n; i++) {
        j = i + span;
        var qbij = QB[i][j];
        if (qbij === 0) continue;

        var prob = (Qf[i] * qbij * term[i][j] * Qr[j + 1]) / Q;

        var kLo = Math.max(0, i - MAX_INTERNAL - 1);
        for (k = kLo; k < i; k++) {
          var li = i - k - 1;
          var lHi = Math.min(n - 1, j + 1 + (MAX_INTERNAL - li));
          for (l = j + 1; l <= lHi; l++) {
            var outer = P[k][l];
            if (outer === 0) continue;
            var ri = l - j - 1;
            var w;
            if (li === 0) {
              if (ri === 0) w = wStack[k][l];
              else if (ri === 1) w = wBulge1R[k][l];
              else w = wBulge[ri] * term[k][l] * term[i][j];
            } else if (ri === 0) {
              if (li === 1) w = wBulge1L[k][l];
              else w = wBulge[li] * term[k][l] * term[i][j];
            } else {
              w = wInternal[li][ri] * intTerm[k][l] * intTerm[i][j];
            }
            prob += (outer / QB[k][l]) * w * qbij;
          }
        }

        var multi = 0;
        for (k = 0; k < i; k++) {
          var branchAcc = accBranch[k][j];
          var openAcc = accUnpaired[k][j];
          if (branchAcc === 0 && openAcc === 0) continue;
          var leftM = k + 1 <= i - 1 ? QM[k + 1][i - 1] : 0;
          if (leftM) multi += leftM * (openAcc + branchAcc);
          if (branchAcc) multi += zPow[i - 1 - k] * branchAcc;
        }
        if (multi) prob += qbij * mlBranch * term[i][j] * multi;

        if (prob <= 0) continue;
        P[i][j] = Math.min(prob, 1);
        pending.push([i, j, (P[i][j] / qbij) * mlOpen * term[i][j]]);
      }
    }

    var unpaired = new Array(n).fill(1);
    for (i = 0; i < n; i++) {
      var used = 0;
      for (j = i + MIN_LOOP + 1; j < n; j++) used += P[i][j];
      for (k = 0; k < i - MIN_LOOP; k++) used += P[k][i];
      unpaired[i] = Math.max(0, Math.min(1, 1 - used));
    }

    return {
      sequence: seq,
      partitionFunction: Q,
      ensembleEnergy: -RT37_PF * Math.log(Q),
      pairProb: P,
      unpaired: unpaired
    };
  }

  /* Mean probability that a base is unpaired across the window [start, end). */
  function meanUnpaired(unpaired, start, end) {
    var total = 0;
    var count = 0;
    for (var i = start; i < end && i < unpaired.length; i++) {
      if (i < 0) continue;
      total += unpaired[i];
      count++;
    }
    return count ? total / count : 0;
  }

  // ------------------------------------------------------- CRISPR guide layer

  var SEED_LENGTH = 8; // PAM-proximal seed for SpCas9

  function gcPercent(seq) {
    if (!seq.length) return 0;
    var gc = 0;
    for (var i = 0; i < seq.length; i++) {
      if (seq[i] === "G" || seq[i] === "C") gc++;
    }
    return (gc / seq.length) * 100;
  }

  function pairedPositions(pairs) {
    var set = {};
    for (var i = 0; i < pairs.length; i++) {
      set[pairs[i][0]] = true;
      set[pairs[i][1]] = true;
    }
    return set;
  }

  /*
   * Turn a fold into guide-design features. The biological premise: a spacer
   * that pairs with itself, especially across the PAM-proximal seed, has less
   * of that seed available to interrogate the DNA target.
   */
  function scoreGuide(sequence, model) {
    var fold = model === "zuker" ? zuker(sequence) : nussinov(sequence);
    var seq = fold.sequence;
    var n = seq.length;
    var paired = pairedPositions(fold.pairs);

    var seedStart = Math.max(0, n - SEED_LENGTH);
    var seedBases = n - seedStart;
    var seedPaired = 0;
    for (var i = seedStart; i < n; i++) {
      if (paired[i]) seedPaired++;
    }
    var seedAccessibility = seedBases ? (seedBases - seedPaired) / seedBases : 0;

    var gc = gcPercent(seq);
    var pairCount = fold.pairs.length;
    var selfPairFraction = n ? (2 * pairCount) / n : 0;

    var warnings = [];
    var dna = seq.replace(/U/g, "T");
    if (n !== 20) warnings.push("not 20 nt");
    if (gc < 30) warnings.push("low GC");
    if (gc > 75) warnings.push("high GC");
    if (dna.indexOf("TTTT") !== -1) warnings.push("poly-T (U6 termination risk)");
    if (n >= SEED_LENGTH && seedAccessibility < 0.5) warnings.push("seed mostly paired");

    return {
      sequence: seq,
      structure: fold.structure,
      pairs: fold.pairs,
      dp: fold.dp || null,
      energy: fold.energy != null ? fold.energy : null,
      score: fold.score != null ? fold.score : pairCount,
      pairCount: pairCount,
      gc: gc,
      seedStart: seedStart,
      seedAccessibility: seedAccessibility,
      selfPairFraction: selfPairFraction,
      warnings: warnings
    };
  }

  /*
   * Per-position agreement between two structures. Returns one char per base:
   *   "=" both models agree, "1" only model A pairs it, "2" only model B pairs it.
   */
  function comparePairs(n, pairsA, pairsB) {
    function keys(pairs) {
      var m = {};
      for (var i = 0; i < pairs.length; i++) m[pairs[i][0] + ":" + pairs[i][1]] = true;
      return m;
    }
    var a = keys(pairsA);
    var b = keys(pairsB);
    var shared = 0;
    for (var key in a) {
      if (b[key]) shared++;
    }
    var union = Object.keys(a).length + Object.keys(b).length - shared;
    var onlyA = pairedPositions(pairsA);
    var onlyB = pairedPositions(pairsB);
    var track = [];
    for (var i = 0; i < n; i++) {
      var inA = !!onlyA[i];
      var inB = !!onlyB[i];
      if (inA && inB) track.push("=");
      else if (inA) track.push("1");
      else if (inB) track.push("2");
      else track.push(".");
    }
    return {
      shared: shared,
      onlyA: Object.keys(a).length - shared,
      onlyB: Object.keys(b).length - shared,
      agreement: union ? shared / union : 1,
      track: track.join("")
    };
  }

  return {
    normalize: normalize,
    canPair: canPair,
    dotBracket: dotBracket,
    nussinov: nussinov,
    zuker: zuker,
    mccaskill: mccaskill,
    meanUnpaired: meanUnpaired,
    scoreGuide: scoreGuide,
    comparePairs: comparePairs,
    gcPercent: gcPercent,
    pairedPositions: pairedPositions,
    stackEnergy: stackEnergy,
    SEED_LENGTH: SEED_LENGTH,
    MAX_LEN: MAX_LEN
  };
})();
