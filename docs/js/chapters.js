/*
 * The data-driven chapters: Guide Analyzer, Designer, Dataset, Findings, and
 * the reference checks.
 *
 * Data provenance:
 *   data/guides.json      guide activity from Doench 2014/2016 via CRISPOR;
 *                         the seedOpenness and seedEnsemble accessibility
 *                         columns are computed by this repository's own
 *                         zuker.py and mccaskill.py.
 *   data/references.json  published reference structures (textbook examples
 *                         and the accepted yeast tRNA-Phe cloverleaf).
 *   data/vienna.json      ViennaRNA package output, precomputed because
 *                         ViennaRNA does not run in the browser.
 *   data/model.json       ridge weights fitted by export_model.py.
 *   data/study.json       the four experiments in chapter 05, reduced from
 *                         analysis_outputs/study_results.json.
 *
 * All folding scored on this page is this project's own code (js/rna-fold.js).
 */

window.Chapters = (function () {
  "use strict";

  var UI = window.FoldUI;
  var $ = UI.$;
  var $$ = UI.$$;

  // A 262 bp stretch of human EMX1, the standard worked example for Cas9.
  var EMX1 =
    "GAGTCCGAGCAGAAGAAGAAGGGCTCCCATCACATCAACCGGTGGCGCATTGCCACGAAGCAGGCCAATGGGGAGGACATCGATGTCACCTCCAATGACT" +
    "AGGGTGGGCAACCACAAACCCACGAGGGCAGAGTGCTGCTTGCTGCTGGCCAGGCCCCTGCGTGGGCCCAAGCTGGACTCTGGCCACTCCCTGGCCAGGC" +
    "TTTGGGGAGGCCTGGAGTCATGGCCCCACAGGGCTTGAAGCCCGGGGCCGCCATTGACAGAG";

  var MAX_TARGET = 5000;
  var TABLE_LIMIT = 250;

  var state = { guides: [], designer: [], checkDone: false };

  function pct(x) {
    return Math.round(x * 100) + "%";
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

  function summaryStat(label, value, hint) {
    return (
      '<div class="summary-stat"><span>' +
      label +
      "</span><strong>" +
      value +
      "</strong>" +
      (hint ? '<em class="stat-hint">' + hint + "</em>" : "") +
      "</div>"
    );
  }

  // ============================================================== ANALYZER

  function analyzerModel() {
    var checked = $('input[name="analyzerModel"]:checked');
    return checked ? checked.value : "zuker";
  }

  /*
   * The analyzer folds the spacer and, when the scaffold is on, a 96-letter
   * sgRNA through both the MFE recurrence and the partition function. That is
   * long enough to be felt, so the work is deferred a frame behind a busy
   * state rather than blocking the keystroke that triggered it.
   */
  var analyzerTimer = null;

  function scheduleAnalyzer() {
    clearTimeout(analyzerTimer);
    UI.setBusy($("#analyzeButton"), $("#analyzerStatus"), true, "folding");
    analyzerTimer = setTimeout(runAnalyzer, 180);
  }

  function runAnalyzer() {
    clearTimeout(analyzerTimer);
    var error = $("#analyzerError");
    var status = $("#analyzerStatus");
    var button = $("#analyzeButton");
    error.textContent = "";

    var raw = $("#spacerInput").value;
    var dna = GuideTools.toDna(raw);
    $("#spacerCount").textContent = dna.length;

    if (!dna.length) {
      UI.setBusy(button, status, false);
      status.textContent = "empty";
      status.classList.remove("is-warn");
      $("#analyzerFlags").innerHTML = "";
      $("#analyzerInterpretation").textContent = "No spacer entered yet. A 20 letter spacer goes in the field above.";
      UI.setupCanvas($("#analyzerCanvas"));
      $("#analyzerLegend").innerHTML = "";
      return;
    }

    var info;
    try {
      info = GuideTools.analyzeGuide(raw, {
        algorithm: analyzerModel(),
        withScaffold: $("#withScaffold").checked
      });
    } catch (err) {
      UI.setBusy(button, status, false);
      error.textContent = err.message;
      status.textContent = "check input";
      status.classList.add("is-warn");
      return;
    }

    UI.setBusy(button, status, false);
    $("#seedOpennessValue").textContent = pct(info.seedOpenness);
    $("#analyzerGc").textContent = Math.round(info.gcPercent) + "%";
    $("#analyzerEnergy").textContent = info.energy != null ? info.energy.toFixed(2) : "n/a";

    var prediction = dna.length === 20 ? GuideTools.predictEfficiency(dna) : null;
    $("#predEfficiency").textContent = prediction ? pct(prediction.value) : "n/a";

    $("#analyzerSpacerSeq").textContent = info.spacer;
    $("#analyzerSpacerStruct").textContent = info.spacerStructure;

    status.textContent = info.flags.length ? info.flags.length + " to check" : "looks fine";
    status.classList.toggle("is-warn", info.flags.length > 0);

    $("#analyzerFlags").innerHTML = info.flags.length
      ? info.flags
          .map(function (f) {
            return '<span class="warn-tag">' + UI.escapeHtml(f) + "</span>";
          })
          .join("")
      : '<span class="warn-tag ok">Nothing to flag</span>';

    var spacerLen = info.spacer.length;
    UI.drawArcDiagram($("#analyzerCanvas"), info.full, info.pairs, null, {
      seedFrom: Math.max(0, spacerLen - GuideTools.SEED_LENGTH),
      seedTo: spacerLen,
      colorA: analyzerModel() === "zuker" ? UI.ARC_B : UI.ARC_A,
      divider: $("#withScaffold").checked ? spacerLen : null
    });
    $("#analyzerLegend").innerHTML = [
      UI.swatch(analyzerModel() === "zuker" ? UI.ARC_B : UI.ARC_A, analyzerModel() === "zuker" ? "Zuker pairs" : "Nussinov pairs"),
      UI.swatch(UI.SEED, "the seed, last 8 letters")
    ].join("");

    var open = Math.round(info.seedOpenness * 100);
    var verdict;
    if (open >= 75) verdict = "The model leaves most seed positions unpaired.";
    else if (open >= 50) verdict = "The model leaves at least half of the seed unpaired.";
    else if (open >= 25) verdict = "The model pairs most seed positions.";
    else verdict = "The model pairs nearly the entire seed.";

    var scaffoldNote = $("#withScaffold").checked
      ? " The standard scaffold was appended before folding."
      : " Only the bare spacer was folded; adding the scaffold can change the predicted structure.";

    var lengthNote =
      dna.length === 20
        ? ""
        : " This is not a standard 20 nt spacer, so the display is a fold diagnostic only and no efficiency estimate is shown.";
    $("#analyzerInterpretation").textContent =
      open + "% of the seed is free. " + verdict + scaffoldNote + lengthNote;

    renderEnsemble(info);
  }

  /*
   * The ensemble panel. Everything here comes from the partition function, so
   * the numbers are probabilities rather than counts taken off one structure.
   */
  function renderEnsemble(info) {
    var status = $("#ensembleStatus");
    var canvas = $("#ensembleCanvas");
    var ensemble = info.ensemble;

    if (!ensemble) {
      status.textContent = "unavailable";
      status.classList.add("is-warn");
      ["#ensembleSeed", "#ensembleSingle", "#ensembleMean", "#ensembleEnergy"].forEach(function (id) {
        $(id).textContent = "n/a";
      });
      UI.setupCanvas(canvas);
      $("#ensembleLegend").innerHTML = "";
      $("#ensembleInterpretation").textContent =
        "The partition function did not run for this sequence.";
      return;
    }

    status.textContent = "computed";
    status.classList.remove("is-warn");

    var spacerLen = info.spacer.length;
    var seedFrom = Math.max(0, spacerLen - GuideTools.SEED_LENGTH);

    $("#ensembleSeed").textContent = pct(ensemble.seedOpenness);
    $("#ensembleSingle").textContent = pct(info.seedOpenness);
    $("#ensembleMean").textContent = pct(ensemble.meanUnpaired);
    $("#ensembleEnergy").textContent = ensemble.energy.toFixed(2);

    UI.drawProfile(canvas, ensemble.unpaired, {
      seedFrom: seedFrom,
      seedTo: spacerLen,
      color: UI.ARC_A
    });
    $("#ensembleLegend").innerHTML = [
      UI.swatch(UI.ARC_A, "probability the base is unpaired"),
      UI.swatch(UI.SEED, "the seed, last 8 letters")
    ].join("");

    var gap = ensemble.seedOpenness - info.seedOpenness;
    var comparison;
    if (Math.abs(gap) < 0.05) {
      comparison =
        "Both views agree here, so the single structure is a fair summary of the ensemble for this guide.";
    } else if (gap > 0) {
      comparison =
        "The ensemble reports a more open seed than the single structure does, by " +
        Math.round(Math.abs(gap) * 100) +
        " points. Some of the pairs in that one structure are formed only part of the time.";
    } else {
      comparison =
        "The ensemble reports a less open seed than the single structure does, by " +
        Math.round(Math.abs(gap) * 100) +
        " points. Bases left free in that one structure spend some of their time paired in others.";
    }

    var lowest = 0;
    for (var i = 1; i < ensemble.unpaired.length; i++) {
      if (ensemble.unpaired[i] < ensemble.unpaired[lowest]) lowest = i;
    }
    var context = $("#withScaffold").checked
      ? "folded as the full sgRNA"
      : "folded as a bare spacer, which is not the molecule that exists in a cell";

    $("#ensembleInterpretation").textContent =
      "Measured across the whole ensemble " +
      context +
      ", the seed is " +
      pct(ensemble.seedOpenness) +
      " open and the least available position is " +
      (lowest + 1) +
      " at " +
      pct(ensemble.unpaired[lowest]) +
      " unpaired. " +
      comparison;
  }

  function renderModelCard() {
    var meta = GuideTools.modelMeta();
    var statusPill = $("#modelStatus");
    if (!meta) {
      statusPill.textContent = "unavailable";
      statusPill.classList.add("is-warn");
      $("#modelStats").innerHTML = "";
      $("#modelNote").textContent =
        "The efficiency model did not load, so the predicted number is hidden. Everything else on this page still works.";
      return;
    }
    statusPill.textContent = "loaded";
    statusPill.classList.remove("is-warn");

    var delta = meta.heldOutSpearman - meta.heldOutSpearmanSequenceOnly;
    $("#modelStats").innerHTML = [
      stat("Held-out correlation", "ρ " + meta.heldOutSpearman.toFixed(3), "within gene, whole genes held out"),
      stat("Sequence only", "ρ " + meta.heldOutSpearmanSequenceOnly.toFixed(3), "without folding features"),
      stat("G/C baseline", "ρ " + meta.baselineGcSpearman.toFixed(3), "within gene, G/C alone"),
      stat("Weights", meta.params.toLocaleString(), "ridge regression")
    ].join("");

    $("#modelNote").textContent =
      "This is a ridge regression fitted on " +
      meta.trainedOn +
      ". Its inputs are position-specific bases, neighboring pairs, G/C content, and two folding " +
      "features of the whole sgRNA that this page recomputes on demand: how much of the seed " +
      "the minimum-free-energy structure leaves unpaired, and that structure's free energy. " +
      "Accuracy is measured by " +
      meta.validation +
      ", across " +
      meta.nGenes +
      " genes. Guides targeting the same gene are not independent, so a random split would " +
      "flatter the model. Held out that way it reaches ρ ≈ " +
      meta.heldOutSpearman.toFixed(3) +
      " for ranking guides within a target, against ρ ≈ " +
      meta.baselineGcSpearman.toFixed(3) +
      " for G/C alone. Adding the folding features to the sequence-only model changes that by " +
      (Math.abs(delta) < 0.0005 ? "less than 0.001" : (delta >= 0 ? "+" : "") + delta.toFixed(3)) +
      ". Structure adds nothing here once sequence is accounted for. That is the study's main " +
      "finding, not a shortcoming of this model.";
  }

  function initAnalyzer() {
    $("#analyzerForm").addEventListener("submit", function (e) {
      e.preventDefault();
      runAnalyzer();
    });
    $("#spacerInput").addEventListener("input", scheduleAnalyzer);
    $("#withScaffold").addEventListener("change", runAnalyzer);
    $$('input[name="analyzerModel"]').forEach(function (r) {
      r.addEventListener("change", runAnalyzer);
    });
    $$("[data-spacer]").forEach(function (b) {
      b.addEventListener("click", function () {
        $("#spacerInput").value = b.dataset.spacer;
        runAnalyzer();
      });
    });
    $("#analyzerClear").addEventListener("click", function () {
      $("#spacerInput").value = "";
      $("#spacerInput").focus();
      runAnalyzer();
    });
    if ("ResizeObserver" in window) {
      new ResizeObserver(function () {
        runAnalyzer();
      }).observe($("#analyzerCanvas"));
    }
    runAnalyzer();
  }

  // ============================================================== DESIGNER

  function renderDesigner(rows) {
    var body = $("#designerBody");
    if (!rows.length) {
      body.innerHTML = '<tr><td colspan="9" class="table-empty">No NGG sites found in that sequence.</td></tr>';
      $("#designerStats").innerHTML = "";
      return;
    }
    var clean = rows.filter(function (r) {
      return !r.flags.length;
    }).length;

    $("#designerStats").innerHTML = [
      stat("Guide sites", rows.length, "NGG PAM, both strands"),
      stat("Best predicted", rows[0].efficiency != null ? pct(rows[0].efficiency) : "n/a", "top ranked guide"),
      stat("Clean guides", clean, "nothing flagged")
    ].join("");

    body.innerHTML = rows
      .map(function (r, i) {
        var flags = r.flags.length
          ? r.flags
              .map(function (f) {
                return '<span class="warn-tag tiny">' + UI.escapeHtml(f) + "</span>";
              })
              .join(" ")
          : '<span class="flag-clean">clean</span>';
        return (
          "<tr><td>" +
          (i + 1) +
          '</td><td class="mono">' +
          r.spacer +
          "</td><td>" +
          r.strand +
          "</td><td>" +
          r.position +
          "</td><td>" +
          Math.round(r.gcPercent) +
          "%</td><td>" +
          pct(r.seedOpenness) +
          '</td><td class="num-strong">' +
          (r.efficiency != null ? pct(r.efficiency) : "n/a") +
          "</td><td>" +
          flags +
          '</td><td><button class="chip tiny" type="button" data-open-guide="' +
          r.spacer +
          '">open</button></td></tr>'
        );
      })
      .join("");

    $$("[data-open-guide]").forEach(function (b) {
      b.addEventListener("click", function () {
        $("#spacerInput").value = b.dataset.openGuide;
        runAnalyzer();
        var target = document.getElementById("analyzer");
        window.scrollTo({ top: target.getBoundingClientRect().top + window.pageYOffset - 76, behavior: "smooth" });
      });
    });
  }

  function skeletonRows(columns, rows) {
    var cells = "";
    for (var c = 0; c < columns; c++) {
      cells += '<td><span class="skeleton skeleton-line"></span></td>';
    }
    var out = "";
    for (var r = 0; r < rows; r++) out += "<tr>" + cells + "</tr>";
    return out;
  }

  function renderDesignerSkeleton() {
    $("#designerStats").innerHTML = "";
    $("#designerBody").innerHTML = skeletonRows(9, 6);
  }

  function runDesigner() {
    var error = $("#designerError");
    var status = $("#designerStatus");
    error.textContent = "";

    var dna = GuideTools.toDna($("#targetInput").value).replace(/[^ACGT]/g, "");
    if (!dna.length) {
      error.textContent = "No target DNA entered yet.";
      return;
    }
    if (dna.length > MAX_TARGET) {
      error.textContent = "That is " + dna.length.toLocaleString() + " bp. This demo handles up to " + MAX_TARGET.toLocaleString() + ".";
      return;
    }

    var button = $("#designerButton");
    UI.setBusy(button, status, true, "folding");
    renderDesignerSkeleton();

    // Folding every candidate takes a moment, so let the busy state paint
    // before the main thread is taken.
    setTimeout(function () {
      var t0 = performance.now();
      state.designer = GuideTools.rankGuides(dna);
      var ms = Math.round(performance.now() - t0);
      renderDesigner(state.designer);
      UI.setBusy(button, status, false);
      status.textContent = state.designer.length + " found · " + ms + " ms";
    }, 30);
  }

  function initDesigner() {
    $("#designerForm").addEventListener("submit", function (e) {
      e.preventDefault();
      runDesigner();
    });
    $("#targetInput").addEventListener("input", function () {
      $("#targetCount").textContent = GuideTools.toDna($("#targetInput").value).replace(/[^ACGT]/g, "").length;
    });
    $("#designerExample").addEventListener("click", function () {
      $("#targetInput").value = EMX1;
      $("#targetCount").textContent = EMX1.length;
      runDesigner();
    });
    $("#designerClear").addEventListener("click", function () {
      $("#targetInput").value = "";
      $("#targetCount").textContent = "0";
      $("#designerBody").innerHTML = '<tr><td colspan="9" class="table-empty">Paste a target sequence and press Find guides.</td></tr>';
      $("#designerStats").innerHTML = "";
      $("#designerStatus").textContent = "ready";
    });
    $("#designerCsv").addEventListener("click", function () {
      if (!state.designer.length) return;
      var lines = ["rank,spacer_dna,strand,position,pam,gc_percent,seed_openness,predicted_efficiency,flags"];
      state.designer.forEach(function (r, i) {
        lines.push(
          [
            i + 1,
            r.spacer,
            r.strand,
            r.position,
            r.pam || "",
            r.gcPercent.toFixed(0),
            r.seedOpenness.toFixed(3),
            r.efficiency != null ? r.efficiency.toFixed(3) : "",
            r.flags.join(";") || "clean"
          ]
            .map(UI.csvCell)
            .join(",")
        );
      });
      UI.downloadCsv("designed_guides.csv", lines);
    });
  }

  // =============================================================== DATASET

  function binned(guides) {
    // Openness is a fraction of 8 letters, so it lands on 0, 1/8, 2/8 ...
    var bins = {};
    guides.forEach(function (g) {
      var key = Math.round(g.seedOpenness * 8) / 8;
      if (!bins[key]) bins[key] = { sum: 0, n: 0 };
      bins[key].sum += g.activity;
      bins[key].n++;
    });
    return Object.keys(bins)
      .map(Number)
      .sort(function (a, b) {
        return a - b;
      })
      .map(function (k) {
        return { openness: k, mean: bins[k].sum / bins[k].n, n: bins[k].n };
      });
  }

  function drawDatasetChart() {
    var canvas = $("#datasetChart");
    if (!canvas) return;
    var env = UI.setupCanvas(canvas);
    var ctx = env.ctx;
    if (!state.guides.length) return;

    var rows = binned(state.guides);
    var padL = 58;
    var padR = 22;
    var padT = 22;
    var padB = 42;
    var plotW = env.w - padL - padR;
    var plotH = env.h - padT - padB;

    ctx.font = "11px Inter, system-ui, sans-serif";
    ctx.textBaseline = "middle";

    // y axis 0..100%
    for (var t = 0; t <= 5; t++) {
      var y = padT + plotH - (plotH * t) / 5;
      ctx.strokeStyle = "rgba(233,238,242,0.06)";
      ctx.beginPath();
      ctx.moveTo(padL, y);
      ctx.lineTo(padL + plotW, y);
      ctx.stroke();
      ctx.fillStyle = "rgba(132,150,163,0.85)";
      ctx.textAlign = "right";
      ctx.fillText(t * 20 + "%", padL - 10, y);
    }

    var barW = plotW / rows.length;
    rows.forEach(function (r, i) {
      var h = plotH * Math.max(0, Math.min(1, r.mean));
      var x = padL + i * barW + barW * 0.22;
      var w = barW * 0.56;
      var y = padT + plotH - h;
      var grad = ctx.createLinearGradient(0, y, 0, padT + plotH);
      grad.addColorStop(0, "#58c4d4");
      grad.addColorStop(1, "rgba(88,196,212,0.22)");
      ctx.fillStyle = grad;
      ctx.fillRect(x, y, w, h);

      ctx.fillStyle = "rgba(233,238,242,0.95)";
      ctx.textAlign = "center";
      ctx.font = "11px Inter, system-ui, sans-serif";
      ctx.fillText(Math.round(r.mean * 100) + "%", x + w / 2, y - 9);
      ctx.fillStyle = "rgba(132,150,163,0.9)";
      ctx.fillText(Math.round(r.openness * 100) + "%", x + w / 2, padT + plotH + 14);
      ctx.fillStyle = "rgba(132,150,163,0.7)";
      ctx.font = "10px Inter, system-ui, sans-serif";
      ctx.fillText("n=" + r.n, x + w / 2, padT + plotH + 27);
    });

    ctx.fillStyle = "rgba(132,150,163,0.95)";
    ctx.font = "11px Inter, system-ui, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Seed openness", padL + plotW / 2, env.h - 6);
    ctx.save();
    ctx.translate(13, padT + plotH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText("Mean efficiency", 0, 0);
    ctx.restore();
  }

  function renderDataset() {
    var guides = state.guides;
    if (!guides.length) {
      $("#datasetCount").textContent = "unavailable";
      $("#datasetFinding").textContent =
        "The dataset did not load, so this chapter is empty. Everything else on this page still works.";
      $("#datasetBody").innerHTML = '<tr><td colspan="6" class="table-empty">Dataset unavailable.</td></tr>';
      return;
    }

    var openness = guides.map(function (g) {
      return g.seedOpenness;
    });
    var activity = guides.map(function (g) {
      return g.activity;
    });
    var gc = guides.map(function (g) {
      return g.gcPercent;
    });

    var rhoOpen = GuideTools.spearman(openness, activity);
    var rhoGc = GuideTools.spearman(gc, activity);
    var genes = {};
    guides.forEach(function (g) {
      genes[g.gene] = true;
    });

    $("#datasetCount").textContent = guides.length.toLocaleString() + " guides";
    $("#datasetStats").innerHTML = [
      summaryStat("Guides", guides.length.toLocaleString(), "with lab measured efficiency"),
      summaryStat("Genes", Object.keys(genes).length.toLocaleString(), "across two screens"),
      summaryStat("Openness vs efficiency", "ρ " + rhoOpen.toFixed(2), "little association"),
      summaryStat("G/C vs efficiency", "ρ " + rhoGc.toFixed(2), "weak but positive")
    ].join("");

    $("#datasetFinding").textContent =
      "Each bar shows the mean measured efficiency for guides at one level of seed openness. Across " +
      guides.length.toLocaleString() +
      " guides, the relationship is nearly flat (Spearman ρ ≈ " +
      rhoOpen.toFixed(2) +
      "). Seed openness therefore adds little information on its own in this pooled dataset. G/C content has a slightly stronger but still weak relationship (ρ ≈ " +
      rhoGc.toFixed(2) +
      "). One caveat applies to both numbers: they pool guides across 18 genes that differ a lot in how editable they are. Measured within each gene, which is the comparison a designer actually makes when choosing among guides for one target, G/C drops to about ρ = -0.02. Its apparent usefulness here is a difference between genes rather than between guides. Chapter 05 works through this properly.";

    $("#datasetSource").textContent =
      "This table contains " +
      guides.length.toLocaleString() +
      " guides with measured editing efficiency, pooled from Doench 2014 and 2016 through CRISPOR / Haeussler 2016. Activity is represented as a within-dataset percentile. Seed openness is the fraction of the last 8 spacer bases left unpaired in the minimum-free-energy structure of the full sgRNA, precomputed with this project's Python folder (zuker.py) because folding all " +
      guides.length.toLocaleString() +
      " guides in the browser would be slow. The table shows the first " +
      TABLE_LIMIT +
      " rows.";

    $("#datasetBody").innerHTML = guides
      .slice(0, TABLE_LIMIT)
      .map(function (g) {
        return (
          '<tr><td class="mono">' +
          UI.escapeHtml(g.spacer) +
          "</td><td>" +
          UI.escapeHtml(g.gene) +
          "</td><td>" +
          Math.round(g.gcPercent) +
          "%</td><td>" +
          pct(g.seedOpenness) +
          '</td><td class="num-strong">' +
          pct(g.activity) +
          '</td><td><button class="chip tiny" type="button" data-open-guide="' +
          UI.escapeHtml(g.spacer) +
          '">open</button></td></tr>'
        );
      })
      .join("");

    $$("#datasetBody [data-open-guide]").forEach(function (b) {
      b.addEventListener("click", function () {
        $("#spacerInput").value = b.dataset.openGuide;
        runAnalyzer();
        var target = document.getElementById("analyzer");
        window.scrollTo({ top: target.getBoundingClientRect().top + window.pageYOffset - 76, behavior: "smooth" });
      });
    });

    drawDatasetChart();
  }

  function skeletonStats(count) {
    var out = "";
    for (var i = 0; i < count; i++) {
      out +=
        '<div class="summary-stat" aria-hidden="true">' +
        '<span class="skeleton skeleton-line" style="width: 60%"></span>' +
        '<strong><span class="skeleton skeleton-line" style="width: 45%; height: 28px; margin-top: 8px"></span></strong>' +
        "</div>";
    }
    return out;
  }

  function initDataset() {
    $("#datasetCsv").addEventListener("click", function () {
      if (!state.guides.length) return;
      var lines = ["spacer,gene,gc_percent,seed_openness,activity_percentile"];
      state.guides.forEach(function (g) {
        lines.push([g.spacer, g.gene, g.gcPercent, g.seedOpenness.toFixed(3), g.activity.toFixed(3)].map(UI.csvCell).join(","));
      });
      UI.downloadCsv("guide_dataset.csv", lines);
    });

    if ("ResizeObserver" in window) {
      new ResizeObserver(function () {
        drawDatasetChart();
      }).observe($("#datasetChart"));
    }

    $("#datasetStats").innerHTML = skeletonStats(4);

    fetch("data/guides.json")
      .then(function (r) {
        if (!r.ok) throw new Error("no dataset");
        return r.json();
      })
      .then(function (data) {
        // Accept either the bare array or the {source, guides} wrapper.
        state.guides = Array.isArray(data) ? data : data.guides || [];
        renderDataset();
        UI.observeReveal();
      })
      .catch(function () {
        state.guides = [];
        renderDataset();
      });
  }

  // ========================================================= CHECK OUR WORK

  function drawViennaScatter(points, r) {
    var canvas = $("#viennaScatter");
    if (!canvas) return;
    var env = UI.setupCanvas(canvas);
    var ctx = env.ctx;
    if (!points.length) return;

    var padL = 58;
    var padR = 20;
    var padT = 20;
    var padB = 44;
    var plotW = env.w - padL - padR;
    var plotH = env.h - padT - padB;

    var xs = points.map(function (p) {
      return p.vienna;
    });
    var ys = points.map(function (p) {
      return p.mine;
    });
    var xMin = Math.min.apply(null, xs) - 1;
    var xMax = Math.max.apply(null, xs) + 1;
    var yMin = Math.min.apply(null, ys) - 1;
    var yMax = Math.max.apply(null, ys) + 1;

    function px(v) {
      return padL + ((v - xMin) / (xMax - xMin)) * plotW;
    }
    function py(v) {
      return padT + plotH - ((v - yMin) / (yMax - yMin)) * plotH;
    }

    ctx.font = "11px Inter, system-ui, sans-serif";
    ctx.textBaseline = "middle";
    for (var t = 0; t <= 4; t++) {
      var gx = padL + (plotW * t) / 4;
      var gy = padT + (plotH * t) / 4;
      ctx.strokeStyle = "rgba(233,238,242,0.05)";
      ctx.beginPath();
      ctx.moveTo(gx, padT);
      ctx.lineTo(gx, padT + plotH);
      ctx.moveTo(padL, gy);
      ctx.lineTo(padL + plotW, gy);
      ctx.stroke();

      ctx.fillStyle = "rgba(132,150,163,0.85)";
      ctx.textAlign = "center";
      ctx.fillText(Math.round(xMin + ((xMax - xMin) * t) / 4), gx, padT + plotH + 15);
      ctx.textAlign = "right";
      ctx.fillText(Math.round(yMax - ((yMax - yMin) * t) / 4), padL - 10, gy);
    }

    points.forEach(function (p) {
      ctx.fillStyle = "rgba(88,196,212,0.85)";
      ctx.beginPath();
      ctx.arc(px(p.vienna), py(p.mine), 3.4, 0, Math.PI * 2);
      ctx.fill();
    });

    ctx.fillStyle = "rgba(132,150,163,0.95)";
    ctx.textAlign = "center";
    ctx.fillText("ViennaRNA energy (kcal/mol)", padL + plotW / 2, env.h - 8);
    ctx.save();
    ctx.translate(13, padT + plotH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText("This model, kcal/mol", 0, 0);
    ctx.restore();

    ctx.fillStyle = "#58c4d4";
    ctx.textAlign = "left";
    ctx.font = "12px Inter, system-ui, sans-serif";
    ctx.fillText("r = " + r.toFixed(2), padL + 10, padT + 12);
  }

  function runCheck() {
    if (!state.checkDone) {
      $("#checkStats").innerHTML = skeletonStats(4);
    }
    if (state.checkDone) return;
    state.checkDone = true;

    Promise.all([
      fetch("data/references.json").then(function (r) {
        return r.ok ? r.json() : null;
      }),
      fetch("data/vienna.json").then(function (r) {
        return r.ok ? r.json() : null;
      })
    ])
      .then(function (results) {
        var refs = results[0];
        var vienna = results[1];
        if (!refs || !vienna) throw new Error("no reference data");

        var scored = refs.map(function (ref) {
          var mz = RNA.zuker(ref.seq);
          var mn = RNA.nussinov(ref.seq);
          var v = vienna.references[ref.id];
          return {
            ref: ref,
            zukerF1: GuideTools.f1Score(mz.structure, ref.structure),
            nussinovF1: GuideTools.f1Score(mn.structure, ref.structure),
            viennaF1: v ? GuideTools.f1Score(v.structure, ref.structure) : 0,
            zukerStruct: mz.structure,
            nussinovStruct: mn.structure,
            viennaStruct: v ? v.structure : "",
            zukerEnergy: mz.energy,
            viennaEnergy: v ? v.energy : null
          };
        });

        var avg = function (key) {
          return scored.reduce(function (s, r) {
            return s + r[key];
          }, 0) / scored.length;
        };

        $("#checkStats").innerHTML = [
          summaryStat("Answer key shapes", scored.length, "solved in a lab"),
          summaryStat("Zuker-style", avg("zukerF1").toFixed(2), "mean pair F1"),
          summaryStat("Nussinov", avg("nussinovF1").toFixed(2), "mean pair F1"),
          summaryStat("ViennaRNA", avg("viennaF1").toFixed(2), "reference software")
        ].join("");

        $("#checkCards").innerHTML = scored
          .map(function (s) {
            function row(label, struct, f1, cls) {
              return (
                '<div class="check-row"><div class="check-row-head"><span class="dot-tag ' +
                (cls || "") +
                '">' +
                label +
                '</span><span class="check-f1' +
                (f1 >= 0.99 ? " is-perfect" : f1 < 0.6 ? " is-poor" : "") +
                '">F1 ' +
                f1.toFixed(2) +
                "</span></div><code>" +
                UI.escapeHtml(struct) +
                "</code></div>"
              );
            }
            return (
              '<article class="card check-card"><div class="card-head"><div><p class="kicker">' +
              s.ref.seq.length +
              " letters</p><h3>" +
              UI.escapeHtml(s.ref.name) +
              '</h3></div></div><p class="guide-notes">' +
              UI.escapeHtml(s.ref.note) +
              '</p><div class="check-scroll"><div class="check-row"><div class="check-row-head"><span class="dot-tag">Known shape</span></div><code>' +
              UI.escapeHtml(s.ref.structure) +
              "</code></div>" +
              row("Zuker-style", s.zukerStruct, s.zukerF1, "tag-b") +
              row("Nussinov", s.nussinovStruct, s.nussinovF1, "tag-a") +
              row("ViennaRNA", s.viennaStruct, s.viennaF1, "") +
              "</div></article>"
            );
          })
          .join("");

        // 30 real sgRNAs folded by this project's Zuker model and by ViennaRNA.
        var points = vienna.guideSample.map(function (g) {
          var mz = RNA.zuker(g.full);
          return { mine: mz.energy, vienna: g.energy, f1: GuideTools.f1Score(mz.structure, g.structure) };
        });
        var r = GuideTools.pearson(
          points.map(function (p) {
            return p.vienna;
          }),
          points.map(function (p) {
            return p.mine;
          })
        );
        var meanF1 =
          points.reduce(function (s, p) {
            return s + p.f1;
          }, 0) / points.length;

        $("#viennaCorr").textContent = "r = " + r.toFixed(2) + " · " + points.length + " guides";
        drawViennaScatter(points, r);

        $("#viennaNote").textContent =
          "Each point is one sgRNA evaluated by ViennaRNA and by the simplified Zuker-style model. Their energies correlate at r = " +
          r.toFixed(2) +
          ", and the predicted structures recover " +
          Math.round(meanF1 * 100) +
          "% of ViennaRNA's pairs on average (pair F1). The browser model produces systematically less-negative energies because its loop and wobble terms are simplified. The methods are useful to compare, but they are not interchangeable. The ViennaRNA values are precomputed with the ViennaRNA package because it does not run in this browser demo.";

        if ("ResizeObserver" in window) {
          new ResizeObserver(function () {
            drawViennaScatter(points, r);
          }).observe($("#viennaScatter"));
        }
        UI.observeReveal();
      })
      .catch(function () {
        $("#checkStats").innerHTML = "";
        $("#viennaCorr").textContent = "unavailable";
        $("#viennaNote").textContent =
          "The reference data did not load, so this chapter is empty. Everything else on this page still works.";
      });
  }

  function initCheck() {
    // Folding 30 sgRNAs with the scaffold is about a second of work, so hold
    // off until the reader is actually heading for this chapter.
    var section = document.getElementById("check");
    if (!section || !("IntersectionObserver" in window)) {
      runCheck();
      return;
    }
    var observer = new IntersectionObserver(
      function (entries) {
        if (entries[0].isIntersecting) {
          runCheck();
          observer.disconnect();
        }
      },
      { rootMargin: "300px" }
    );
    observer.observe(section);

    // Safety net. IntersectionObserver only delivers on a rendering frame, so
    // in a tab that is not painting it may never fire and the chapter would sit
    // on "computing" forever. Compute anyway after a few seconds; runCheck
    // guards itself, so whichever path wins, the work happens once.
    setTimeout(function () {
      if (!state.checkDone) {
        runCheck();
        observer.disconnect();
      }
    }, 5000);
  }

  // ------------------------------------------------------------------ init

  // ============================================================== FINDINGS

  /*
   * Chapter 05 renders data/study.json, written by export_site_data.py from the
   * full analysis. Nothing is recomputed in the browser: the cross-validation
   * behind these numbers holds out whole genes and takes about twenty-five
   * minutes, so the site reads the result rather than reproducing it.
   */
  var study = null;

  function signed(value, places) {
    if (value == null) return "n/a";
    return (value >= 0 ? "+" : "") + value.toFixed(places == null ? 3 : places);
  }

  function interval(low, high) {
    if (low == null || high == null) return "n/a";
    return "[" + signed(low) + ", " + signed(high) + "]";
  }

  function crossesZero(low, high) {
    return low == null || high == null || (low <= 0 && high >= 0);
  }

  function renderMeasures() {
    var doench = study.datasets.doench;
    var body = $("#measuresBody");
    body.innerHTML = doench.correlations
      .map(function (row) {
        var solid = !crossesZero(row.withinLow, row.withinHigh);
        return (
          "<tr><td>" +
          UI.escapeHtml(row.label) +
          '</td><td class="mono">' +
          signed(row.spearman) +
          '</td><td class="mono' +
          (solid ? " num-strong" : "") +
          '">' +
          signed(row.withinGene) +
          '</td><td class="mono subtle">' +
          interval(row.withinLow, row.withinHigh) +
          "</td></tr>"
        );
      })
      .join("");

    var best = doench.correlations.reduce(function (a, b) {
      return Math.abs(b.withinGene) > Math.abs(a.withinGene) ? b : a;
    });
    var ensemble = doench.correlations.filter(function (r) {
      return r.feature === "seed_ensemble_full";
    })[0];
    var nussinov = doench.correlations.filter(function (r) {
      return r.feature === "seed_nussinov";
    })[0];

    $("#measuresFinding").textContent =
      "Measured on the full sgRNA across the whole ensemble, the seed correlates with efficiency at " +
      signed(ensemble.withinGene) +
      ", against " +
      signed(nussinov.withinGene) +
      " for the original pair-counting measure. The sharper instrument does recover more, and in the " +
      "direction the hypothesis predicts, but a correlation of " +
      ensemble.withinGene.toFixed(3) +
      " accounts for under half a percent of the variation in activity. The strongest measure in the " +
      "table is " +
      best.label.toLowerCase() +
      " at " +
      signed(best.withinGene) +
      ". Experiment 03 tests whether any of this survives once sequence is taken into account.";
  }

  function positionScreen() {
    var checked = $('input[name="positionScreen"]:checked');
    return checked ? checked.value : "doench";
  }

  function renderPositions() {
    var key = positionScreen();
    var block = study.datasets[key];
    var canvas = $("#positionCanvas");

    UI.drawSignedSeries(
      canvas,
      [
        { color: UI.ARC_A, points: block.positions.spacer },
        { color: UI.ARC_B, points: block.positions.full }
      ],
      { seedFrom: 12, seedTo: 20 }
    );
    $("#positionLegend").innerHTML = [
      UI.swatch(UI.ARC_A, "spacer folded alone"),
      UI.swatch(UI.ARC_B, "spacer folded in the full sgRNA"),
      UI.swatch(UI.SEED, "the conventional seed, positions 13 to 20")
    ].join("");

    var full = block.positions.full;
    var strongest = full.reduce(function (a, b) {
      return Math.abs(b.spearman) > Math.abs(a.spearman) ? b : a;
    });
    var seedMean = 0;
    var restMean = 0;
    full.forEach(function (p) {
      if (p.position >= 13) seedMean += Math.abs(p.spearman) / 8;
      else restMean += Math.abs(p.spearman) / 12;
    });
    var where =
      seedMean > restMean
        ? "the conventional seed carries more of the signal than the rest of the spacer"
        : "the signal sits outside the conventional seed, in the half furthest from the PAM";

    $("#positionFinding").textContent =
      "In " +
      block.label +
      ", folded as the full sgRNA, the strongest position is " +
      strongest.position +
      " at " +
      signed(strongest.spearman) +
      ". Averaged over the window, " +
      where +
      " (" +
      seedMean.toFixed(3) +
      " inside against " +
      restMean.toFixed(3) +
      " outside). Switching between the two screens shows the difficulty: they place their strongest " +
      "positions at opposite ends of the spacer, which is not what one mechanism acting on one region " +
      "would produce.";
  }

  function renderIncremental() {
    var order = ["doench", "crisprscan"];
    var body = $("#incrementalBody");
    body.innerHTML = order
      .map(function (key) {
        var block = study.datasets[key];
        var full = block.models.filter(function (m) {
          return m.name === "structure_all";
        })[0];
        var solid = !crossesZero(full.deltaLow, full.deltaHigh);
        return (
          "<tr><td>" +
          UI.escapeHtml(block.label) +
          '</td><td class="mono">' +
          block.baseline.spearman.toFixed(3) +
          '</td><td class="mono">' +
          full.spearman.toFixed(3) +
          '</td><td class="mono' +
          (solid ? " num-strong" : "") +
          '">' +
          signed(full.delta) +
          '</td><td class="mono subtle">' +
          interval(full.deltaLow, full.deltaHigh) +
          "</td></tr>"
        );
      })
      .join("");

    var human = study.datasets.doench;
    var fish = study.datasets.crisprscan;
    var humanDelta = human.models.filter(function (m) {
      return m.name === "structure_all";
    })[0];
    var fishDelta = fish.models.filter(function (m) {
      return m.name === "structure_all";
    })[0];

    function verdict(model) {
      return crossesZero(model.deltaLow, model.deltaHigh)
        ? "an interval that includes zero"
        : "an interval that excludes zero";
    }

    // The ViennaRNA comparison, when present, is what settles the zebrafish case.
    var viennaFish = null;
    if (study.validation && study.validation.datasets && study.validation.datasets.crisprscan) {
      viennaFish = study.validation.datasets.crisprscan.incremental.models.filter(function (m) {
        return m.name === "vienna_structure_all";
      })[0];
    }

    var sentence =
      "In human cells the sequence-only model reaches " +
      human.baseline.spearman.toFixed(3) +
      " on guides it never saw, and folding changes that by " +
      signed(humanDelta.delta) +
      ", " + verdict(humanDelta) +
      ". In zebrafish the baseline is higher at " +
      fish.baseline.spearman.toFixed(3) +
      " and folding changes it by " +
      signed(fishDelta.delta) +
      ", " + verdict(fishDelta) + ".";

    if (viennaFish) {
      sentence +=
        " Those folding features come from this project's simplified energy model. Recomputing " +
        "them with ViennaRNA under standard Turner parameters raises the zebrafish gain to " +
        signed(viennaFish.delta) +
        " " + interval(viennaFish.deltaLow, viennaFish.deltaHigh) +
        ", so the approximations were hiding a real effect in that screen. It is carried by " +
        "overall folding stability rather than by seed accessibility, and it runs opposite to " +
        "the original idea: more stable folding accompanies higher activity, not less.";
    }
    $("#incrementalFinding").textContent = sentence;
  }

  function renderTransfer() {
    var body = $("#transferBody");
    body.innerHTML = study.transfer
      .map(function (row) {
        var parts = row.transfer.split("->");
        var solid = !crossesZero(row.deltaLow, row.deltaHigh);
        return (
          "<tr><td>" +
          UI.escapeHtml(parts[0].trim()) +
          "</td><td>" +
          UI.escapeHtml(parts[1].trim()) +
          '</td><td class="mono">' +
          row.baseline.toFixed(3) +
          '</td><td class="mono">' +
          row.structure.toFixed(3) +
          '</td><td class="mono' +
          (solid ? " num-strong" : "") +
          '">' +
          signed(row.delta) +
          "</td></tr>"
        );
      })
      .join("");

    var toFish = study.transfer.filter(function (r) {
      return r.transfer.indexOf("-> CRISPRscan") !== -1;
    })[0];
    var toHuman = study.transfer.filter(function (r) {
      return r.transfer.indexOf("CRISPRscan ->") === 0;
    })[0];

    $("#transferFinding").textContent =
      "Accuracy falls sharply whenever a model leaves the screen it was trained on. A model trained on " +
      "human-cell data reaches " +
      toFish.baseline.toFixed(3) +
      " on zebrafish guides, and in the other direction only " +
      toHuman.baseline.toFixed(3) +
      ". The zebrafish folding gain does not survive the crossing either. Guide activity depends heavily " +
      "on the system it is measured in, which is the practical reason a rule found in one screen should " +
      "not be applied to another without checking.";
  }

  /*
   * The ViennaRNA comparison. Present only when validate_vienna.py has been run,
   * so the card stays hidden rather than showing empty tables on a checkout that
   * has not produced it.
   */
  function renderValidation() {
    var card = $("#viennaCard");
    if (!card) return;
    var validation = study.validation;
    if (!validation || !validation.datasets || !validation.datasets.doench) {
      card.hidden = true;
      return;
    }
    card.hidden = false;
    $("#viennaValidationStatus").textContent = "checked";

    var block = validation.datasets.doench;
    var byLabel = {};
    block.correlations.forEach(function (row) {
      byLabel[row.label] = row;
    });

    $("#viennaBody").innerHTML = block.agreement
      .map(function (row) {
        var match = byLabel[row.label] || {};
        return (
          "<tr><td>" + UI.escapeHtml(row.label) +
          '</td><td class="mono">' + row.pearson.toFixed(2) +
          '</td><td class="mono">' + signed(match.custom) +
          '</td><td class="mono">' + signed(match.vienna) +
          "</td></tr>"
        );
      })
      .join("");

    var worst = block.agreement.reduce(function (a, b) {
      return b.pearson < a.pearson ? b : a;
    });
    var best = block.agreement.reduce(function (a, b) {
      return b.pearson > a.pearson ? b : a;
    });
    var shifts = block.correlations.map(function (row) {
      return Math.abs(row.vienna - row.custom);
    });
    var largest = Math.max.apply(null, shifts);

    $("#viennaValidationFinding").textContent =
      "The two implementations agree most closely on " +
      best.label.toLowerCase() +
      " (r = " + best.pearson.toFixed(2) + ") and least closely on " +
      worst.label.toLowerCase() +
      " (r = " + worst.pearson.toFixed(2) +
      "). Ensemble measures agree better than single-structure ones, which is expected: " +
      "one optimal structure can flip between two near-tied alternatives when a parameter " +
      "changes, while an average over all structures moves smoothly. What matters for the " +
      "conclusions is that swapping in the reference parameters shifts no correlation with " +
      "activity by more than " + largest.toFixed(3) + ".";

    var body = $("#viennaSensitivityBody");
    if (body) {
      body.innerHTML = block.sensitivity
        .map(function (row) {
          return (
            "<tr><td>" + UI.escapeHtml(row.setting) +
            '</td><td class="mono">' + signed(row.seedEnsemble) +
            '</td><td class="mono">' + signed(row.meanUnpaired) +
            '</td><td class="mono">' + signed(row.ensembleEnergy) +
            "</td></tr>"
          );
        })
        .join("");
    }

    var global37 = block.sensitivity[0];
    var local = block.sensitivity[block.sensitivity.length - 1];
    $("#viennaSensitivityFinding").textContent =
      "Folding is temperature dependent, so the same quantities were recomputed at 25 and " +
      "42 degrees, and again with long-range pairs forbidden, which approximates a molecule " +
      "folding as it is transcribed. Seed accessibility moves from " +
      signed(global37.seedEnsemble) + " to " + signed(local.seedEnsemble) +
      " across those settings. None of them turns the result into a usable predictor.";
  }

  function renderStudy() {
    renderMeasures();
    renderPositions();
    renderIncremental();
    renderTransfer();
    renderValidation();
    $("#studySource").textContent =
      "All figures come from analysis_outputs/study_results.json, produced by run_study.py over " +
      study.datasets.doench.n.toLocaleString() +
      " guides in " +
      study.datasets.doench.nGenes +
      " genes and " +
      study.datasets.crisprscan.n.toLocaleString() +
      " guides in " +
      study.datasets.crisprscan.nGenes +
      " genes. Accuracy is measured within each gene by cross-validation that holds out whole genes, " +
      "and every interval comes from resampling genes rather than individual guides.";
  }

  function initFindings() {
    var status = $("#findingsStatus");
    if (!status) return;

    UI.$$('input[name="positionScreen"]').forEach(function (input) {
      input.addEventListener("change", function () {
        if (study) renderPositions();
      });
    });

    if ("ResizeObserver" in window) {
      var canvas = $("#positionCanvas");
      if (canvas) {
        new ResizeObserver(function () {
          if (study) renderPositions();
        }).observe(canvas);
      }
    }

    fetch("data/study.json")
      .then(function (response) {
        if (!response.ok) throw new Error("study data unavailable");
        return response.json();
      })
      .then(function (json) {
        study = json;
        status.textContent = "loaded";
        status.classList.remove("is-warn");
        renderStudy();
      })
      .catch(function () {
        status.textContent = "unavailable";
        status.classList.add("is-warn");
        $("#measuresBody").innerHTML =
          '<tr><td colspan="4" class="table-empty">Study results did not load.</td></tr>';
      });
  }

  // ------------------------------------------------------------------ init

  function init() {
    initAnalyzer();
    initDesigner();
    initDataset();
    initFindings();
    initCheck();

    GuideTools.loadModel("data/model.json").then(function () {
      renderModelCard();
      runAnalyzer();
    });
  }

  /*
   * Tab hooks. Each panel's charts already resize themselves from a
   * ResizeObserver when the panel stops being display:none, but that lands a
   * frame later; redrawing here keeps a freshly opened tab from flashing an
   * empty canvas. Every call is guarded because a reader can open any tab
   * before its data has finished loading.
   */
  function redraw(id) {
    try {
      if (id === "analyzer") runAnalyzer();
      else if (id === "dataset" && state.guides.length) drawDatasetChart();
      else if (id === "findings" && study) renderPositions();
    } catch (err) {
      /* a panel that is not ready yet simply redraws when its data arrives */
    }
  }

  /* The Check panel folds 30 sgRNAs on demand. It used to start when the
     reader scrolled near it; with tabs, opening the tab is that signal. */
  function onTabShown(id) {
    if (id === "check") runCheck();
  }

  return { init: init, runAnalyzer: runAnalyzer, redraw: redraw, onTabShown: onTabShown };
})();
