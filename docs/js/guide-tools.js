/*
 * CRISPR guide tooling: seed analysis, the efficiency model, PAM scanning,
 * and accuracy scoring against known structures.
 *
 * Data provenance:
 *   data/model.json       ridge weights fitted by this repository's
 *                         export_model.py on the 4,685 pooled guides.
 *   data/guides.json      Doench 2014/2016 activity via CRISPOR; accessibility
 *                         columns computed by zuker.py and mccaskill.py.
 *   data/vienna.json      ViennaRNA package output, precomputed because
 *                         ViennaRNA does not run in a browser.
 *   data/references.json  published reference structures.
 *
 * The folding here is this project's own code (js/rna-fold.js), and the model
 * below uses only features that code can reproduce in the browser.
 *
 * Seed openness is measured on the fold of the whole sgRNA, not the bare
 * spacer, because the scaffold can pair with the spacer and close the seed off.
 * A single optimal structure is a coarse thing to read openness from: several
 * structures often tie, and a different traceback lands on a different answer.
 * analyzeGuide therefore also runs the partition function (RNA.mccaskill) and
 * returns the probability each base is unpaired across the whole ensemble. The
 * efficiency model below still takes the single-structure value, because that
 * is the feature it was fitted on.
 *
 * Accuracy figures for the model live in data/model.json under "meta" and are
 * measured by cross-validation that holds out whole genes.
 */

window.GuideTools = (function () {
  "use strict";

  // Standard SpCas9 sgRNA scaffold that follows the 20 nt spacer.
  var SCAFFOLD = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC";
  var SEED_LENGTH = 8; // PAM-proximal end of the spacer
  var NFEAT = 386;
  var ORDER = "ACGT";

  var model = null;

  function toRna(seq) {
    return String(seq || "").replace(/\s+/g, "").toUpperCase().replace(/T/g, "U");
  }
  function toDna(seq) {
    return String(seq || "").replace(/\s+/g, "").toUpperCase().replace(/U/g, "T");
  }

  function gcPercent(seq) {
    var s = toDna(seq);
    if (!s.length) return 0;
    var gc = 0;
    for (var i = 0; i < s.length; i++) if (s[i] === "G" || s[i] === "C") gc++;
    return (gc / s.length) * 100;
  }

  /*
   * Openness of the seed: the share of the last 8 spacer bases left unpaired.
   * Measured on the fold of the WHOLE molecule (spacer plus scaffold), because
   * the scaffold can pair with the spacer and close the seed off.
   */
  function seedOpenness(structure, spacerLength) {
    var start = Math.max(0, spacerLength - SEED_LENGTH);
    var segment = structure.slice(start, spacerLength);
    if (!segment.length) return 0;
    var open = 0;
    for (var i = 0; i < segment.length; i++) if (segment[i] === ".") open++;
    return open / segment.length;
  }

  function designFlags(spacerRna, gc, openness) {
    var flags = [];
    var dna = toDna(spacerRna);
    if (dna.length !== 20) flags.push("not 20 letters long");
    if (gc < 30) flags.push("low G/C (under 30%)");
    if (gc > 75) flags.push("high G/C (over 75%)");
    if (dna.indexOf("TTTT") !== -1) flags.push("has a TTTT run");
    if (openness < 0.5) flags.push("seed is mostly folded up");
    return flags;
  }

  /*
   * Fold a guide and report what matters for design.
   * withScaffold defaults to true: that is the molecule that actually exists.
   */
  function analyzeGuide(spacer, options) {
    options = options || {};
    var algorithm = options.algorithm || "zuker";
    var withScaffold = options.withScaffold !== false;
    var wantEnsemble = options.ensemble !== false;

    var spacerRna = toRna(spacer);
    var full = withScaffold ? spacerRna + SCAFFOLD : spacerRna;
    var fold = algorithm === "nussinov" ? RNA.nussinov(full) : RNA.zuker(full);

    var spacerStructure = fold.structure.slice(0, spacerRna.length);
    var openness = seedOpenness(fold.structure, spacerRna.length);
    var gc = gcPercent(spacerRna);

    var spacerPairs = 0;
    for (var i = 0; i < fold.pairs.length; i++) {
      if (fold.pairs[i][0] < spacerRna.length && fold.pairs[i][1] < spacerRna.length) spacerPairs++;
    }

    /*
     * The ensemble view of the same molecule. The fold above is one structure,
     * so every base is scored a hard paired or unpaired; the partition function
     * reports the probability each base is left unpaired across every structure
     * the sequence can adopt, weighted by energy. The two often disagree.
     */
    var ensemble = null;
    if (wantEnsemble) {
      try {
        var pf = RNA.mccaskill(full);
        var seedStart = Math.max(0, spacerRna.length - SEED_LENGTH);
        ensemble = {
          unpaired: pf.unpaired.slice(0, spacerRna.length),
          seedOpenness: RNA.meanUnpaired(pf.unpaired, seedStart, spacerRna.length),
          meanUnpaired: RNA.meanUnpaired(pf.unpaired, 0, spacerRna.length),
          energy: pf.ensembleEnergy
        };
      } catch (err) {
        ensemble = null; // too long to fold without scaling; the rest still works
      }
    }

    return {
      spacer: spacerRna,
      full: full,
      algorithm: algorithm,
      structure: fold.structure,
      spacerStructure: spacerStructure,
      pairs: fold.pairs,
      energy: fold.energy != null ? fold.energy : null,
      seedOpenness: openness,
      ensemble: ensemble,
      gcPercent: gc,
      selfPairFraction: spacerRna.length ? (2 * spacerPairs) / spacerRna.length : 0,
      flags: designFlags(spacerRna, gc, openness)
    };
  }

  // ------------------------------------------------------- efficiency model

  function loadModel(url) {
    return fetch(url || "data/model.json")
      .then(function (response) {
        if (!response.ok) throw new Error("model unavailable");
        return response.json();
      })
      .then(function (json) {
        model = json;
        return json;
      })
      .catch(function () {
        model = null;
        return null;
      });
  }

  function modelMeta() {
    return model ? model.meta : null;
  }

  /*
   * Ridge regression, trained by export_model.py on 4,685 measured guides.
   *
   * Feature layout, matching model_features.py:
   *   position[p * 4 + b]        base b (ACGT) at spacer position p, 20 x 4
   *   dinucleotide[d * 16 + ...] neighbouring pair d, 19 x 16
   *   gc                         G/C fraction of the spacer
   *   structure.seedOpenness     share of the seed left unpaired in the MFE
   *                              structure of the WHOLE sgRNA
   *   structure.mfeEnergyScaled  that structure's free energy, divided by 10
   *
   * Only the two folding features the browser can compute itself are used, so a
   * prediction here is reproducible from what is on screen. The model predicts
   * an activity percentile directly, so the linear score is clamped to [0, 1]
   * rather than squashed through a sigmoid.
   */
  function predictEfficiency(spacer) {
    if (!model) return null;
    var dna = toDna(spacer);
    if (!/^[ACGT]+$/.test(dna) || dna.length !== 20) return null;

    var fold = analyzeGuide(dna, { algorithm: "zuker", withScaffold: true, ensemble: false });
    var openness = fold.seedOpenness;
    var energy = fold.energy != null ? fold.energy : 0;
    var gc = gcPercent(dna) / 100;

    var score = model.intercept;
    for (var p = 0; p < 20; p++) {
      var b = ORDER.indexOf(dna[p]);
      if (b >= 0) score += model.position[p * 4 + b];
    }
    for (var d = 0; d < 19; d++) {
      var b1 = ORDER.indexOf(dna[d]);
      var b2 = ORDER.indexOf(dna[d + 1]);
      if (b1 >= 0 && b2 >= 0) score += model.dinucleotide[d * 16 + b1 * 4 + b2];
    }
    score += model.gc * gc;
    score += model.structure.seedOpenness * openness;
    score += model.structure.mfeEnergyScaled * (energy / model.meta.energyScale);

    return {
      value: Math.max(0, Math.min(1, score)),
      openness: openness,
      energy: energy
    };
  }

  // ------------------------------------------------------------ PAM scanning

  function reverseComplement(dna) {
    var map = { A: "T", C: "G", G: "C", T: "A" };
    var out = "";
    for (var i = dna.length - 1; i >= 0; i--) out += map[dna[i]] || "N";
    return out;
  }

  /*
   * Every SpCas9 site in a target: a 20 nt spacer immediately followed by an
   * NGG PAM, on both strands. Positions are reported on the forward strand.
   */
  function findGuideSites(targetDna) {
    var dna = toDna(targetDna).replace(/[^ACGT]/g, "");
    var sites = [];
    var seen = {};

    function scan(seq, strand) {
      for (var i = 0; i + 23 <= seq.length; i++) {
        var pam = seq.slice(i + 20, i + 23);
        if (pam[1] !== "G" || pam[2] !== "G") continue;
        var spacer = seq.slice(i, i + 20);
        if (spacer.indexOf("N") !== -1) continue;
        // Report 1-based coordinates on the forward strand, the usual convention.
        var position = (strand === "+" ? i : seq.length - (i + 23)) + 1;
        var key = spacer + strand + position;
        if (seen[key]) continue;
        seen[key] = true;
        sites.push({ spacer: spacer, pam: pam, strand: strand, position: position });
      }
    }

    scan(dna, "+");
    scan(reverseComplement(dna), "-");
    return sites;
  }

  /*
   * Score every site and rank by predicted efficiency. Folding each guide with
   * its scaffold is the slow part, so callers should show a busy state.
   */
  function rankGuides(targetDna) {
    var sites = findGuideSites(targetDna);
    var rows = sites.map(function (site) {
      var info = analyzeGuide(site.spacer, { algorithm: "nussinov", withScaffold: true, ensemble: false });
      var prediction = predictEfficiency(site.spacer);
      return {
        spacer: site.spacer,
        strand: site.strand,
        position: site.position,
        pam: site.pam,
        gcPercent: info.gcPercent,
        seedOpenness: info.seedOpenness,
        efficiency: prediction ? prediction.value : null,
        // Seed openness gets its own column here, so repeating it as a flag
        // just marks almost every guide and tells the reader nothing.
        flags: info.flags.filter(function (f) {
          return f.indexOf("seed") === -1;
        })
      };
    });
    rows.sort(function (a, b) {
      return (b.efficiency || 0) - (a.efficiency || 0);
    });
    return rows;
  }

  // --------------------------------------------------------- accuracy scoring

  function pairSet(dotBracket) {
    var stack = [];
    var set = {};
    for (var i = 0; i < dotBracket.length; i++) {
      if (dotBracket[i] === "(") stack.push(i);
      else if (dotBracket[i] === ")") {
        var j = stack.pop();
        if (j !== undefined) set[j + "-" + i] = true;
      }
    }
    return set;
  }

  /*
   * F1 against a known structure: 1.0 means every real pair was found and no
   * extra ones invented. Precision and recall are over base pairs.
   */
  function f1Score(predicted, truth) {
    var a = pairSet(predicted);
    var b = pairSet(truth);
    var aKeys = Object.keys(a);
    var bKeys = Object.keys(b);
    var hit = 0;
    for (var i = 0; i < aKeys.length; i++) if (b[aKeys[i]]) hit++;
    var precision = aKeys.length ? hit / aKeys.length : bKeys.length ? 0 : 1;
    var recall = bKeys.length ? hit / bKeys.length : aKeys.length ? 0 : 1;
    if (!precision && !recall) return 0;
    return (2 * precision * recall) / (precision + recall);
  }

  function pearson(xs, ys) {
    var n = xs.length;
    if (!n) return 0;
    var mx = 0;
    var my = 0;
    var i;
    for (i = 0; i < n; i++) {
      mx += xs[i];
      my += ys[i];
    }
    mx /= n;
    my /= n;
    var num = 0;
    var dx = 0;
    var dy = 0;
    for (i = 0; i < n; i++) {
      num += (xs[i] - mx) * (ys[i] - my);
      dx += (xs[i] - mx) * (xs[i] - mx);
      dy += (ys[i] - my) * (ys[i] - my);
    }
    return dx && dy ? num / Math.sqrt(dx * dy) : 0;
  }

  function rank(values) {
    var order = values
      .map(function (v, i) {
        return { v: v, i: i };
      })
      .sort(function (a, b) {
        return a.v - b.v;
      });
    var ranks = new Array(values.length);
    var i = 0;
    while (i < order.length) {
      var j = i;
      while (j + 1 < order.length && order[j + 1].v === order[i].v) j++;
      var avg = (i + j) / 2 + 1;
      for (var k = i; k <= j; k++) ranks[order[k].i] = avg;
      i = j + 1;
    }
    return ranks;
  }

  function spearman(xs, ys) {
    return pearson(rank(xs), rank(ys));
  }

  return {
    SCAFFOLD: SCAFFOLD,
    SEED_LENGTH: SEED_LENGTH,
    toRna: toRna,
    toDna: toDna,
    gcPercent: gcPercent,
    analyzeGuide: analyzeGuide,
    loadModel: loadModel,
    modelMeta: modelMeta,
    predictEfficiency: predictEfficiency,
    reverseComplement: reverseComplement,
    findGuideSites: findGuideSites,
    rankGuides: rankGuides,
    f1Score: f1Score,
    pearson: pearson,
    spearman: spearman
  };
})();
