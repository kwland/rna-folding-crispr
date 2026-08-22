/*
 * Core UI: navigation, the Fold Lab, and the shared helpers that the other
 * chapters (js/chapters.js) build on.
 *
 * The page is one scrolling document, so every canvas is always laid out and
 * measurable. Navigation is anchor links plus a scroll spy that lights the
 * active chapter in the nav and the side rail.
 */

(function () {
  "use strict";

  var $ = function (sel) {
    return document.querySelector(sel);
  };
  var $$ = function (sel) {
    return Array.prototype.slice.call(document.querySelectorAll(sel));
  };

  var MODEL_A = "#58c4d4"; // Nussinov
  var ARC_A = "#58c4d4";
  var ARC_B = "#e3a44f";
  var SEED = "#e3a44f";

  var BASE_COLORS = { A: "#67c39a", U: "#e0a458", G: "#7fa9dc", C: "#d98ba0" };

  var CHAPTERS = ["home", "lab", "analyzer", "designer", "dataset", "findings", "check", "learn", "method", "about"];

  var state = { model: "both", nussinov: null, zuker: null, matrixCells: null };
  var viewer = null;
  var hero = null;
  var foldTimer = null;

  // ------------------------------------------------------- shared helpers

  function escapeHtml(text) {
    var div = document.createElement("div");
    div.textContent = String(text);
    return div.innerHTML;
  }

  function csvCell(value) {
    var text = String(value);
    return /[",\n]/.test(text) ? '"' + text.replace(/"/g, '""') + '"' : text;
  }

  function downloadCsv(filename, lines) {
    var blob = new Blob([lines.join("\n")], { type: "text/csv;charset=utf-8" });
    var url = URL.createObjectURL(blob);
    var link = document.createElement("a");
    link.href = url;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  }

  function setupCanvas(canvas) {
    var ratio = Math.min(window.devicePixelRatio || 1, 2);
    var rect = canvas.getBoundingClientRect();
    var w = Math.max(1, rect.width);
    var h = Math.max(1, rect.height);
    var wantW = Math.floor(w * ratio);
    var wantH = Math.floor(h * ratio);
    if (canvas.width !== wantW || canvas.height !== wantH) {
      canvas.width = wantW;
      canvas.height = wantH;
    }
    var ctx = canvas.getContext("2d");
    ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    ctx.clearRect(0, 0, w, h);
    return { ctx: ctx, w: w, h: h };
  }

  /*
   * Arc diagram. Pairs in setA arc above the backbone, setB below, so two
   * folds of the same strand can be compared at a glance.
   */
  function drawArcDiagram(canvas, sequence, setA, setB, options) {
    options = options || {};
    var env = setupCanvas(canvas);
    var ctx = env.ctx;
    var n = sequence.length;
    if (!n) return env;

    var both = !!(setA && setB);
    var pad = 24;
    var baseline = both ? env.h * 0.5 : env.h * 0.76;
    var step = n > 1 ? (env.w - pad * 2) / (n - 1) : 0;

    function xOf(i) {
      return pad + step * i;
    }
    function arc(pairs, color, up, limit) {
      if (!pairs) return;
      ctx.strokeStyle = color;
      ctx.lineWidth = n > 90 ? 1 : 1.5;
      for (var p = 0; p < pairs.length; p++) {
        var x1 = xOf(pairs[p][0]);
        var x2 = xOf(pairs[p][1]);
        var height = Math.min(limit, 12 + (x2 - x1) * 0.42);
        var dir = up ? -1 : 1;
        ctx.beginPath();
        ctx.moveTo(x1, baseline + dir * 5);
        ctx.bezierCurveTo(x1, baseline + dir * height, x2, baseline + dir * height, x2, baseline + dir * 5);
        ctx.stroke();
      }
    }

    arc(setA, options.colorA || ARC_A, true, baseline - 12);
    if (both) arc(setB, options.colorB || ARC_B, false, env.h - baseline - 12);

    ctx.strokeStyle = "rgba(132, 150, 163, 0.3)";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(xOf(0), baseline);
    ctx.lineTo(xOf(n - 1), baseline);
    ctx.stroke();

    var seedFrom = options.seedFrom != null ? options.seedFrom : -1;
    var seedTo = options.seedTo != null ? options.seedTo : -1;
    var r = n > 120 ? 1.8 : n > 90 ? 2.2 : n > 40 ? 3.8 : 6.2;

    for (var i = 0; i < n; i++) {
      var x = xOf(i);
      ctx.beginPath();
      ctx.fillStyle = BASE_COLORS[sequence[i]] || "#8496a3";
      ctx.arc(x, baseline, r, 0, Math.PI * 2);
      ctx.fill();

      if (i >= seedFrom && i < seedTo) {
        ctx.strokeStyle = SEED;
        ctx.lineWidth = 1.3;
        ctx.beginPath();
        ctx.arc(x, baseline, r + 2, 0, Math.PI * 2);
        ctx.stroke();
      }
      if (r >= 6) {
        ctx.fillStyle = "#090d11";
        ctx.font = "600 10px Inter, system-ui, sans-serif";
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(sequence[i], x, baseline + 0.5);
      }
    }

    // Marks where the spacer ends and the scaffold begins.
    if (options.divider != null && options.divider > 0 && options.divider < n) {
      var dx = xOf(options.divider) - step / 2;
      ctx.strokeStyle = "rgba(233, 238, 242, 0.24)";
      ctx.setLineDash([3, 3]);
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(dx, 8);
      ctx.lineTo(dx, env.h - 8);
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = "rgba(175, 188, 198, 0.72)";
      ctx.font = "11px Inter, system-ui, sans-serif";
      ctx.textAlign = "left";
      ctx.fillText("scaffold", dx + 5, 14);
    }
    return env;
  }

  function swatch(color, label) {
    return '<span><i class="swatch" style="background:' + color + '"></i>' + label + "</span>";
  }

  /*
   * Column chart of a per-position value in [0, 1]. Used for the probability
   * that each spacer base is left unpaired across the Boltzmann ensemble. A
   * single predicted structure can only answer that with 0 or 1; the ensemble
   * gives a real number, so the shape of this chart is the thing a single
   * structure cannot show.
   */
  function drawProfile(canvas, values, options) {
    options = options || {};
    var env = setupCanvas(canvas);
    var ctx = env.ctx;
    var n = values.length;
    if (!n) return env;

    var padLeft = 30;
    var padRight = 10;
    var padTop = 12;
    var padBottom = 20;
    var plotW = env.w - padLeft - padRight;
    var plotH = env.h - padTop - padBottom;
    var slot = plotW / n;
    var barW = Math.max(2, slot * 0.62);

    function yOf(v) {
      return padTop + (1 - v) * plotH;
    }

    // Gridlines at 0, 0.5 and 1, with the axis labelled on the left.
    ctx.font = "11px Inter, system-ui, sans-serif";
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    [0, 0.5, 1].forEach(function (v) {
      var y = yOf(v);
      ctx.strokeStyle = "rgba(233,238,242,0.09)";
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(padLeft, y);
      ctx.lineTo(env.w - padRight, y);
      ctx.stroke();
      ctx.fillStyle = "rgba(175,188,198,0.85)";
      ctx.fillText(v.toFixed(1), padLeft - 6, y);
    });

    // Shade the seed window so the assumed region is visible against the rest.
    if (options.seedFrom != null && options.seedTo != null && options.seedTo > options.seedFrom) {
      ctx.fillStyle = "rgba(227, 164, 79, 0.14)";
      ctx.fillRect(padLeft + options.seedFrom * slot, padTop, (options.seedTo - options.seedFrom) * slot, plotH);
    }

    for (var i = 0; i < n; i++) {
      var v = Math.max(0, Math.min(1, values[i]));
      var x = padLeft + i * slot + (slot - barW) / 2;
      var y = yOf(v);
      var inSeed = options.seedFrom != null && i >= options.seedFrom && i < options.seedTo;
      ctx.fillStyle = inSeed ? SEED : options.color || ARC_A;
      ctx.globalAlpha = 0.85;
      ctx.fillRect(x, y, barW, padTop + plotH - y);
      ctx.globalAlpha = 1;
    }

    // Position ticks: first, last, and every fifth in between.
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.fillStyle = "rgba(175,188,198,0.85)";
    for (var t = 0; t < n; t++) {
      var label = t + 1;
      if (label !== 1 && label !== n && label % 5 !== 0) continue;
      ctx.fillText(String(label), padLeft + t * slot + slot / 2, padTop + plotH + 5);
    }
    return env;
  }

  /*
   * Line chart for a signed per-position statistic with confidence intervals,
   * used for the correlation between unpaired probability and editing
   * efficiency at each spacer position. Unlike drawProfile the values can be
   * negative, so the axis is centred on zero and the zero line is drawn solid:
   * whether an interval crosses it is the whole question.
   */
  function drawSignedSeries(canvas, series, options) {
    options = options || {};
    var env = setupCanvas(canvas);
    var ctx = env.ctx;
    if (!series.length || !series[0].points.length) return env;

    var n = series[0].points.length;
    var padLeft = 42;
    var padRight = 14;
    var padTop = 16;
    var padBottom = 30;
    var plotW = env.w - padLeft - padRight;
    var plotH = env.h - padTop - padBottom;

    var limit = 0.05;
    series.forEach(function (s) {
      s.points.forEach(function (p) {
        [p.spearman, p.ciLow, p.ciHigh].forEach(function (v) {
          if (v != null && isFinite(v)) limit = Math.max(limit, Math.abs(v));
        });
      });
    });
    var step = 0.05;
    while (limit / step > 5) step *= 2;
    limit = Math.ceil(limit / step) * step;

    function yOf(v) {
      return padTop + (1 - (v + limit) / (2 * limit)) * plotH;
    }
    function xOf(i) {
      return padLeft + ((i + 0.5) / n) * plotW;
    }

    ctx.font = "11px Inter, system-ui, sans-serif";
    for (var g = -Math.round(limit / step); g <= Math.round(limit / step); g++) {
      var value = g * step;
      var y = yOf(value);
      var zero = g === 0;
      ctx.strokeStyle = zero ? "rgba(233,238,242,0.34)" : "rgba(233,238,242,0.08)";
      ctx.lineWidth = zero ? 1.4 : 1;
      ctx.beginPath();
      ctx.moveTo(padLeft, y);
      ctx.lineTo(env.w - padRight, y);
      ctx.stroke();
      ctx.fillStyle = "rgba(175,188,198,0.85)";
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      ctx.fillText((value >= 0 ? "+" : "") + value.toFixed(2), padLeft - 6, y);
    }

    if (options.seedFrom != null && options.seedTo != null) {
      var x0 = padLeft + (options.seedFrom / n) * plotW;
      var x1 = padLeft + (options.seedTo / n) * plotW;
      ctx.fillStyle = "rgba(227, 164, 79, 0.12)";
      ctx.fillRect(x0, padTop, x1 - x0, plotH);
      ctx.fillStyle = "rgba(227, 164, 79, 0.8)";
      ctx.textAlign = "center";
      ctx.textBaseline = "top";
      ctx.fillText("conventional seed", (x0 + x1) / 2, padTop + 2);
    }

    series.forEach(function (s, index) {
      var shift = (index - (series.length - 1) / 2) * 5;
      ctx.strokeStyle = s.color;
      ctx.lineWidth = 2;
      ctx.globalAlpha = 0.4;
      s.points.forEach(function (p, i) {
        if (p.ciLow == null || p.ciHigh == null) return;
        var x = xOf(i) + shift;
        ctx.beginPath();
        ctx.moveTo(x, yOf(p.ciLow));
        ctx.lineTo(x, yOf(p.ciHigh));
        ctx.stroke();
      });
      ctx.globalAlpha = 0.9;
      ctx.beginPath();
      s.points.forEach(function (p, i) {
        var x = xOf(i) + shift;
        var y = yOf(p.spearman);
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      });
      ctx.stroke();
      s.points.forEach(function (p, i) {
        var x = xOf(i) + shift;
        var y = yOf(p.spearman);
        var solid = p.pAdjusted != null && p.pAdjusted < 0.05;
        ctx.beginPath();
        ctx.arc(x, y, solid ? 4 : 3, 0, Math.PI * 2);
        ctx.fillStyle = solid ? s.color : "rgba(9, 13, 17, 0.92)";
        ctx.fill();
        ctx.strokeStyle = s.color;
        ctx.lineWidth = 1.6;
        ctx.stroke();
      });
      ctx.globalAlpha = 1;
    });

    ctx.fillStyle = "rgba(175,188,198,0.85)";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    for (var t = 0; t < n; t++) {
      var label = t + 1;
      if (label !== 1 && label !== n && label % 5 !== 0) continue;
      ctx.fillText(String(label), xOf(t), padTop + plotH + 6);
    }
    ctx.fillText("spacer position, 1 is the far end from the PAM", padLeft + plotW / 2, padTop + plotH + 18);
    return env;
  }

  window.FoldUI = {
    $: $,
    $$: $$,
    escapeHtml: escapeHtml,
    csvCell: csvCell,
    downloadCsv: downloadCsv,
    setupCanvas: setupCanvas,
    drawArcDiagram: drawArcDiagram,
    drawProfile: drawProfile,
    drawSignedSeries: drawSignedSeries,
    swatch: swatch,
    BASE_COLORS: BASE_COLORS,
    ARC_A: ARC_A,
    ARC_B: ARC_B,
    SEED: SEED
  };

  // ------------------------------------------------------------- navigation

  var activeTab = null;

  /*
   * Show one chapter and hide the rest.
   *
   * A canvas inside a display:none panel measures 0 wide, so anything drawn
   * while hidden ends up with a 1px backing store. Every chart here sizes
   * itself from a ResizeObserver, which fires when the panel becomes visible,
   * but that lands a frame later. Redrawing explicitly on activation avoids a
   * visible flash of empty canvas.
   */
  function activateTab(id, options) {
    options = options || {};
    if (CHAPTERS.indexOf(id) === -1) id = CHAPTERS[0];
    if (id === activeTab && !options.force) return;
    activeTab = id;

    CHAPTERS.forEach(function (cid) {
      var panel = document.getElementById(cid);
      if (!panel) return;
      var on = cid === id;
      panel.classList.toggle("is-active", on);
      panel.hidden = !on;
    });

    $$("[data-nav]").forEach(function (link) {
      var on = link.dataset.nav === id;
      link.classList.toggle("is-active", on);
      if (link.getAttribute("role") === "tab") {
        link.setAttribute("aria-selected", String(on));
        link.tabIndex = on ? 0 : -1;
      }
    });

    if (options.hash !== false && location.hash.slice(1) !== id) {
      history.pushState(null, "", "#" + id);
    }

    redrawPanel(id);
    if (window.Chapters && Chapters.onTabShown) Chapters.onTabShown(id);
  }

  /* Repaint whatever charts live in the panel that just appeared. */
  function redrawPanel(id) {
    requestAnimationFrame(function () {
      if (id === "home" && hero) hero.resize();
      if (id === "lab") {
        drawArcs();
        drawMatrix();
        if (viewer) viewer.draw();
      }
      if (window.Chapters && Chapters.redraw) Chapters.redraw(id);
    });
  }

  function initTabs() {
    $(".masthead").classList.add("is-stuck");

    /* Left and right arrows walk the strip, per the tablist pattern. */
    $$('[role="tab"]').forEach(function (tab) {
      tab.addEventListener("keydown", function (event) {
        var i = CHAPTERS.indexOf(tab.dataset.nav);
        var next = null;
        if (event.key === "ArrowRight") next = CHAPTERS[(i + 1) % CHAPTERS.length];
        if (event.key === "ArrowLeft") next = CHAPTERS[(i - 1 + CHAPTERS.length) % CHAPTERS.length];
        if (event.key === "Home") next = CHAPTERS[0];
        if (event.key === "End") next = CHAPTERS[CHAPTERS.length - 1];
        if (!next) return;
        event.preventDefault();
        activateTab(next);
        var el = document.getElementById("tab-" + next);
        if (el) el.focus();
      });
    });

    window.addEventListener("popstate", function () {
      activateTab(location.hash.slice(1) || CHAPTERS[0], { hash: false });
    });
  }

  function closeNav() {
    $(".mainnav").classList.remove("is-open");
    $(".nav-toggle").setAttribute("aria-expanded", "false");
  }

  /*
   * Busy state. Every interaction that blocks for longer than a frame flips
   * its button to a spinner and its status pill to a pulsing label, so the
   * page never looks frozen and results never appear out of nowhere.
   */
  function setBusy(button, pill, busy, label) {
    if (button) button.classList.toggle("is-loading", busy);
    if (pill) {
      pill.classList.toggle("is-busy", busy);
      if (busy) {
        pill.classList.remove("is-warn");
        pill.textContent = label || "working";
      }
    }
  }

  // ---------------------------------------------------------------- folding

  function foldFailed(status, button, error, message) {
    setBusy(button, status, false);
    error.textContent = message;
    status.textContent = "check input";
    status.classList.add("is-warn");
  }

  function runFold() {
    var error = $("#errorText");
    var status = $("#inputStatus");
    var button = $("#foldButton");
    var seq;

    error.textContent = "";

    try {
      seq = RNA.normalize($("#sequenceInput").value);
    } catch (err) {
      $("#charCount").textContent = String($("#sequenceInput").value.replace(/\s/g, "").length);
      foldFailed(status, button, error, err.message);
      return;
    }

    $("#charCount").textContent = seq.length;

    if (!seq.length) {
      state.nussinov = null;
      state.zuker = null;
      renderAll();
      setBusy(button, status, false);
      status.textContent = "empty";
      status.classList.remove("is-warn");
      return;
    }

    // Both recurrences are cubic, so a 300-letter strand is real work on the
    // main thread. Two frames: one to paint the busy state, one to fold.
    setBusy(button, status, true, "folding");
    requestAnimationFrame(function () {
      requestAnimationFrame(function () {
        try {
          state.nussinov = RNA.nussinov(seq, Number($("#loopLength").value), $("#allowWobble").checked);
          state.zuker = RNA.zuker(seq, $("#allowWobble").checked);
          renderAll();
          setBusy(button, status, false);
          status.classList.remove("is-warn");
          status.textContent = seq.length + " nt";
        } catch (err) {
          foldFailed(status, button, error, err.message);
        }
      });
    });
  }

  function scheduleFold() {
    clearTimeout(foldTimer);
    foldTimer = setTimeout(runFold, 180);
  }

  function activeFold() {
    return state.model === "zuker" ? state.zuker : state.nussinov;
  }

  function renderAll() {
    renderMetrics();
    renderDotRows();
    renderInterpretation();
    renderCompare();
    drawArcs();
    drawMatrix();
    updateViewer();
  }

  function renderMetrics() {
    var n = state.nussinov;
    var z = state.zuker;
    var metrics = $$("#metricGrid .metric");
    metrics[0].classList.toggle("is-dim", state.model === "zuker");
    metrics[0].classList.add("accent-a");
    metrics[1].classList.toggle("is-dim", state.model === "nussinov");
    metrics[1].classList.add("accent-b");

    $("#pairCount").textContent = n ? String(n.pairs.length) : "0";
    $("#mfeValue").textContent = z ? z.energy.toFixed(1) : "0.0";

    var fold = activeFold();
    var seq = fold ? fold.sequence : "";
    $("#gcPercent").textContent = seq.length ? RNA.gcPercent(seq).toFixed(0) + "%" : "0%";

    if (!fold || seq.length < RNA.SEED_LENGTH) {
      $("#seedOpen").textContent = "n/a";
    } else {
      var paired = RNA.pairedPositions(fold.pairs);
      var start = seq.length - RNA.SEED_LENGTH;
      var open = 0;
      for (var i = start; i < seq.length; i++) if (!paired[i]) open++;
      $("#seedOpen").textContent = Math.round((open / RNA.SEED_LENGTH) * 100) + "%";
    }

    var status = $("#inputStatus");
    var list = $("#warningList");
    list.innerHTML = "";
    if (!fold || !seq.length) {
      status.textContent = "ready";
      status.classList.remove("is-warn");
      return;
    }

    // Guide design rules only mean something for a 20 letter spacer. Suppress
    // them everywhere at once, or the pill counts warnings the panel denies.
    if (seq.length !== 20) {
      status.textContent = "ready";
      status.classList.remove("is-warn");
      list.innerHTML = '<span class="warn-tag ok">Guide checks apply to 20 letter spacers. Try the Analyzer.</span>';
      return;
    }

    var scored = GuideTools.analyzeGuide(seq, {
      algorithm: state.model === "zuker" ? "zuker" : "nussinov",
      withScaffold: false
    });
    status.textContent = scored.flags.length ? scored.flags.length + " to check" : "ready";
    status.classList.toggle("is-warn", scored.flags.length > 0);
    list.innerHTML = scored.flags.length
      ? scored.flags
          .map(function (w) {
            return '<span class="warn-tag">' + escapeHtml(w) + "</span>";
          })
          .join("")
      : '<span class="warn-tag ok">Nothing to flag</span>';
  }

  function renderDotRows() {
    var fold = activeFold();
    $("#normalizedSequence").textContent = fold ? fold.sequence : "";
    $("#dotBracket").textContent = state.nussinov ? state.nussinov.structure : "";
    $("#dotBracketZuker").textContent = state.zuker ? state.zuker.structure : "";
    $("#nussinovRow").style.display = state.model === "zuker" ? "none" : "flex";
    $("#zukerRow").style.display = state.model === "nussinov" ? "none" : "flex";
  }

  function renderInterpretation() {
    var el = $("#foldInterpretation");
    var n = state.nussinov;
    var z = state.zuker;
    if (!n || !n.sequence.length) {
      el.textContent = "No sequence entered yet.";
      return;
    }
    var parts = [];
    if (state.model !== "zuker") parts.push("Nussinov finds " + n.pairs.length + " pairs");
    if (state.model !== "nussinov") {
      parts.push("Zuker predicts " + z.energy.toFixed(1) + " kcal/mol with " + z.pairs.length + " pairs");
    }
    parts.push(RNA.gcPercent(n.sequence).toFixed(0) + "% G/C");
    var text = parts.join(" · ") + ".";
    if (state.model === "both" && z.pairs.length < n.pairs.length) {
      text +=
        " The Zuker-style model predicts fewer pairs because it includes energetic costs that Nussinov does not.";
    }
    el.textContent = text;
  }

  function renderCompare() {
    var card = $("#compareCard");
    if (state.model !== "both" || !state.nussinov || !state.nussinov.sequence.length) {
      card.style.display = "none";
      return;
    }
    card.style.display = "";

    var n = state.nussinov;
    var z = state.zuker;
    var cmp = RNA.comparePairs(n.sequence.length, n.pairs, z.pairs);

    $("#agreementPill").textContent = Math.round(cmp.agreement * 100) + "% agree";
    $("#agreementPill").classList.toggle("is-warn", cmp.agreement < 0.5);

    $("#compareStats").innerHTML = [
      stat("Agreed pairs", cmp.shared, "both models picked these"),
      stat("Nussinov only", cmp.onlyA, "pair counting alone"),
      stat("Zuker only", cmp.onlyB, "energy alone"),
      stat("Zuker energy", z.energy.toFixed(1), "kcal/mol")
    ].join("");

    var note;
    if (cmp.shared === 0 && n.pairs.length) {
      note =
        "They agree on no pairs here. The two objective functions favor different structures, so this is a useful warning that the result is model-sensitive.";
    } else if (cmp.onlyA > cmp.onlyB) {
      note =
        "Nussinov predicts " +
        cmp.onlyA +
        " additional pair(s). Pair counting rewards them, while the energy model may reject them when the surrounding structure is unfavorable.";
    } else if (cmp.onlyB > 0) {
      note =
        "The Zuker-style model predicts " +
        cmp.onlyB +
        " pair(s) that do not appear in the Nussinov traceback because the two methods optimize different scores.";
    } else {
      note =
        "Both land on the same shape. That is an encouraging consistency check for this input, although agreement between two simplified models does not establish the real structure.";
    }
    $("#compareNote").textContent = note;
  }

  function stat(label, value, hint) {
    return (
      '<div class="compare-stat"><span>' +
      label +
      "</span><strong>" +
      value +
      "</strong>" +
      (hint ? '<em class="stat-hint">' + hint + "</em>" : "") +
      "</div>"
    );
  }

  // ------------------------------------------------------------ diagrams

  function drawArcs() {
    var fold = state.nussinov;
    var legend = $("#arcLegend");
    if (!fold || !fold.sequence.length) {
      setupCanvas($("#rnaCanvas"));
      legend.innerHTML = "";
      return;
    }
    var seq = fold.sequence;
    var n = seq.length;
    var both = state.model === "both";
    var seedFrom = n >= RNA.SEED_LENGTH ? n - RNA.SEED_LENGTH : n;

    drawArcDiagram(
      $("#rnaCanvas"),
      seq,
      state.model === "zuker" ? state.zuker.pairs : state.nussinov.pairs,
      both ? state.zuker.pairs : null,
      { seedFrom: seedFrom, seedTo: n, colorA: state.model === "zuker" ? ARC_B : ARC_A }
    );

    var items = [];
    if (state.model !== "zuker") items.push(swatch(ARC_A, "Nussinov" + (both ? " (above)" : "")));
    if (state.model !== "nussinov") items.push(swatch(ARC_B, "Zuker" + (both ? " (below)" : "")));
    if (n >= RNA.SEED_LENGTH) items.push(swatch(SEED, "last 8 letters"));
    legend.innerHTML = items.join("");
  }

  function drawMatrix() {
    var env = setupCanvas($("#matrixCanvas"));
    var ctx = env.ctx;
    var useZuker = state.model === "zuker";

    $("#matrixTitle").textContent = useZuker ? "Zuker grid" : "Nussinov grid";
    $(".matrix-legend").firstElementChild.textContent = useZuker ? "less stable" : "fewer pairs";
    $(".matrix-legend").lastElementChild.textContent = useZuker ? "more stable" : "more pairs";

    var fold = useZuker ? state.zuker : state.nussinov;
    if (!fold || !fold.sequence.length) {
      state.matrixCells = null;
      return;
    }
    var n = fold.sequence.length;
    var grid = useZuker ? fold.W : fold.dp;
    if (!grid) return;

    var size = Math.min(env.w, env.h);
    var cell = size / n;
    var ox = (env.w - size) / 2;
    var oy = (env.h - size) / 2;

    // Nussinov counts up from 0; Zuker energies run down from 0.
    var best = 1;
    for (var a = 0; a < n; a++) {
      for (var b = a; b < n; b++) {
        var val = useZuker ? -grid[a][b] : grid[a][b];
        if (val > best) best = val;
      }
    }
    for (var i = 0; i < n; i++) {
      for (var j = 0; j < n; j++) {
        if (j < i) ctx.fillStyle = "#090d11";
        else {
          var v = useZuker ? -grid[i][j] : grid[i][j];
          ctx.fillStyle = ramp(Math.max(0, Math.min(1, v / best)));
        }
        ctx.fillRect(ox + j * cell, oy + i * cell, Math.ceil(cell), Math.ceil(cell));
      }
    }
    if (n <= 60) {
      ctx.strokeStyle = "rgba(233, 238, 242, 0.05)";
      ctx.lineWidth = 0.5;
      for (var k = 0; k <= n; k++) {
        var p = k * cell;
        ctx.beginPath();
        ctx.moveTo(ox + p, oy);
        ctx.lineTo(ox + p, oy + size);
        ctx.moveTo(ox, oy + p);
        ctx.lineTo(ox + size, oy + p);
        ctx.stroke();
      }
    }
    ctx.strokeStyle = SEED;
    ctx.lineWidth = 1.6;
    ctx.strokeRect(ox + (n - 1) * cell, oy, cell, cell);

    state.matrixCells = { n: n, cell: cell, ox: ox, oy: oy, grid: grid, useZuker: useZuker, seq: fold.sequence };
  }

  function ramp(t) {
    var stops = [
      [12, 20, 27],
      [29, 111, 124],
      [79, 216, 232]
    ];
    var seg = t < 0.5 ? 0 : 1;
    var local = t < 0.5 ? t / 0.5 : (t - 0.5) / 0.5;
    var from = stops[seg];
    var to = stops[seg + 1];
    return (
      "rgb(" +
      Math.round(from[0] + (to[0] - from[0]) * local) +
      "," +
      Math.round(from[1] + (to[1] - from[1]) * local) +
      "," +
      Math.round(from[2] + (to[2] - from[2]) * local) +
      ")"
    );
  }

  function updateViewer() {
    if (!viewer) return;
    var fold = activeFold();
    if (!fold || !fold.sequence.length) {
      viewer.setModel("", [], 0);
      return;
    }
    var seedStart =
      fold.sequence.length >= RNA.SEED_LENGTH ? fold.sequence.length - RNA.SEED_LENGTH : fold.sequence.length;
    viewer.setModel(fold.sequence, fold.pairs, seedStart);
  }

  // ------------------------------------------------------------ stack table

  function renderStackTable() {
    var rows = [
      ["A", "U", "A", "U"],
      ["A", "U", "U", "A"],
      ["U", "A", "A", "U"],
      ["C", "G", "U", "A"],
      ["C", "G", "A", "U"],
      ["G", "C", "U", "A"],
      ["G", "C", "A", "U"],
      ["C", "G", "G", "C"],
      ["G", "C", "G", "C"],
      ["G", "C", "C", "G"]
    ];
    var values = rows.map(function (r) {
      return RNA.stackEnergy(r[0], r[1], r[2], r[3]);
    });
    var min = Math.min.apply(null, values);

    $("#stackTable").innerHTML = rows
      .map(function (r, index) {
        var energy = values[index];
        var strength = energy / min;
        var duplex = "5'-" + r[0] + r[2] + "-3'\n3'-" + r[1] + r[3] + "-5'";
        return (
          '<div class="stack-cell" style="background: rgba(88, 196, 212, ' +
          (0.03 + strength * 0.12).toFixed(3) +
          ')"><div class="stack-duplex">' +
          duplex +
          '</div><div class="stack-energy" style="color: ' +
          (strength > 0.75 ? MODEL_A : "#e9eef2") +
          '">' +
          energy.toFixed(2) +
          "</div></div>"
        );
      })
      .join("");
  }

  // ---------------------------------------------------------------- reveal

  function observeReveal() {
    if (!("IntersectionObserver" in window)) return;
    var targets = $$(".card:not(.reveal), .summary-stat:not(.reveal)");
    var observer = new IntersectionObserver(
      function (entries) {
        entries.forEach(function (entry) {
          if (entry.isIntersecting) {
            entry.target.classList.add("is-visible");
            observer.unobserve(entry.target);
          }
        });
      },
      { threshold: 0.04 }
    );
    targets.forEach(function (element) {
      element.classList.add("reveal");
      observer.observe(element);
    });
  }
  window.FoldUI.observeReveal = observeReveal;
  window.FoldUI.setBusy = setBusy;

  // ------------------------------------------------------------------ init

  function randomSpacer() {
    var bases = "ACGT";
    var out = "";
    for (var i = 0; i < 20; i++) out += bases[Math.floor(Math.random() * 4)];
    return out;
  }

  function bindEvents() {
    $("#foldForm").addEventListener("submit", function (event) {
      event.preventDefault();
      runFold();
    });
    $("#sequenceInput").addEventListener("input", scheduleFold);
    $("#loopLength").addEventListener("input", function (event) {
      $("#loopLengthValue").textContent = event.target.value;
      scheduleFold();
    });
    $("#allowWobble").addEventListener("change", runFold);
    $("#clearSequence").addEventListener("click", function () {
      $("#sequenceInput").value = "";
      runFold();
    });
    $("#randomGuide").addEventListener("click", function () {
      $("#sequenceInput").value = randomSpacer();
      runFold();
    });
    $$("[data-example]").forEach(function (button) {
      button.addEventListener("click", function () {
        $("#sequenceInput").value = button.dataset.example;
        runFold();
      });
    });
    $$('input[name="model"]').forEach(function (radio) {
      radio.addEventListener("change", function () {
        state.model = radio.value;
        renderAll();
      });
    });

    $("#copyStructure").addEventListener("click", function () {
      var fold = state.nussinov;
      if (!fold || !fold.sequence.length) return;
      var text = ">fold_lab_result\n" + fold.sequence + "\n";
      if (state.model !== "zuker") text += state.nussinov.structure + "  (Nussinov, " + state.nussinov.pairs.length + " pairs)\n";
      if (state.model !== "nussinov") text += state.zuker.structure + "  (Zuker, " + state.zuker.energy.toFixed(1) + " kcal/mol)\n";
      var button = $("#copyStructure");
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text).then(
          function () {
            button.textContent = "Copied";
            setTimeout(function () {
              button.textContent = "Copy";
            }, 1400);
          },
          function () {
            $("#foldInterpretation").textContent = text;
          }
        );
      } else {
        $("#foldInterpretation").textContent = text;
      }
    });

    $("#spinToggle").addEventListener("click", function () {
      this.textContent = viewer.toggleSpin() ? "Pause spin" : "Resume spin";
    });
    $("#resetView").addEventListener("click", function () {
      viewer.resetView();
      $("#spinToggle").textContent = "Pause spin";
    });

    $("#matrixCanvas").addEventListener("mousemove", function (event) {
      var cells = state.matrixCells;
      var hint = $("#matrixHover");
      if (!cells) return;
      var rect = this.getBoundingClientRect();
      var j = Math.floor((event.clientX - rect.left - cells.ox) / cells.cell);
      var i = Math.floor((event.clientY - rect.top - cells.oy) / cells.cell);
      if (i < 0 || j < 0 || i >= cells.n || j >= cells.n || j < i) {
        hint.textContent = "hover a cell";
        return;
      }
      var v = cells.grid[i][j];
      var label = cells.useZuker ? v.toFixed(1) + " kcal/mol" : v + " pairs";
      hint.textContent =
        "[" + (i + 1) + "," + (j + 1) + "] " + cells.seq.slice(i, j + 1).slice(0, 12) + (j - i > 11 ? "…" : "") + " = " + label;
    });
    $("#matrixCanvas").addEventListener("mouseleave", function () {
      $("#matrixHover").textContent = "hover a cell";
    });

    $$('a[href^="#"]').forEach(function (link) {
      link.addEventListener("click", function (event) {
        var id = link.getAttribute("href").slice(1);
        if (!id || !document.getElementById(id)) return;
        event.preventDefault();
        closeNav();
        if (id === "top") {
          window.scrollTo({ top: 0, behavior: "smooth" });
          return;
        }
        activateTab(id);
        /* Put the panel under the header rather than leaving the reader
           looking at the hero after a click. */
        var main = document.querySelector("main");
        var top = main.getBoundingClientRect().top + window.pageYOffset - 72;
        if (window.pageYOffset < top) window.scrollTo({ top: top, behavior: "smooth" });
      });
    });

    $(".nav-toggle").addEventListener("click", function () {
      var open = $(".mainnav").classList.toggle("is-open");
      this.setAttribute("aria-expanded", String(open));
    });

    /*
     * Drive redraws from a ResizeObserver rather than measuring once at init.
     * At startup the canvas box can still be 0 wide (render blocking webfont
     * CSS), which pinned the backing store to 1px and left the hero blank until
     * the user happened to resize the window.
     */
    if ("ResizeObserver" in window) {
      [
        ["#rnaCanvas", drawArcs],
        ["#matrixCanvas", drawMatrix]
      ].forEach(function (entry) {
        var element = $(entry[0]);
        if (!element) return;
        new ResizeObserver(function () {
          entry[1]();
        }).observe(element);
      });
    }

    var resizeTimer = null;
    window.addEventListener("resize", function () {
      // Backstop for devicePixelRatio changes, which leave the CSS box alone.
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(function () {
        drawArcs();
        drawMatrix();
        if (viewer) viewer.draw();
      }, 140);
    });
  }

  function init() {
    hero = Hero.create($("#heroCanvas"));
    viewer = Viewer3D.create($("#viewerCanvas"));

    bindEvents();
    initTabs();
    renderStackTable();
    runFold();
    observeReveal();

    if (window.Chapters) window.Chapters.init();
    var start = location.hash.slice(1);
    if (start === "tool") start = "lab";
    var deepLinked = CHAPTERS.indexOf(start) !== -1;
    activateTab(deepLinked ? start : CHAPTERS[0], { hash: false, force: true });

    /*
     * A hash that matches a section id makes the browser jump to it before this
     * runs, which parks the page mid-hero. Land the panel under the header
     * instead. Anyone arriving without a hash keeps the full hero.
     */
    if ("scrollRestoration" in history) history.scrollRestoration = "manual";

    function settleScroll() {
      var main = document.querySelector("main");
      var header = parseInt(
        getComputedStyle(document.documentElement).getPropertyValue("--header-h"), 10) || 72;
      /* behavior:instant, because html sets scroll-behavior:smooth and an
         animated jump on page load both looks wrong and loses a race with the
         browser's own hash scroll. */
      window.scrollTo({ top: deepLinked ? Math.max(0, main.offsetTop - header) : 0,
                        behavior: "instant" });
    }
    /* Two passes: once now, once after the first layout, since the browser
       performs its own jump for a hash that matches an element id. */
    settleScroll();
    requestAnimationFrame(settleScroll);
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
