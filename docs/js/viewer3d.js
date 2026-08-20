/*
 * 3D viewer for a predicted RNA secondary structure.
 *
 * There is no external 3D library here. The pipeline is:
 *
 *   1. Relax the molecule with a force-directed layout in 3D. Backbone bonds
 *      and base pairs are springs; every base repels every other base. Stems
 *      settle into ladders and loops open into rings on their own.
 *   2. Project with a perspective divide, sort back-to-front, and draw.
 *
 * This is a SCHEMATIC, not a tertiary-structure prediction. It shows the
 * topology of the predicted secondary structure (which stems and loops exist),
 * laid out in 3D so stems read as helices. It is not an atomic model and
 * carries no claim about real 3D coordinates.
 */

window.Viewer3D = (function () {
  "use strict";

  var BASE_COLORS = {
    A: "#67c39a",
    U: "#e0a458",
    G: "#7fa9dc",
    C: "#d98ba0"
  };

  var BOND_LENGTH = 13;
  var PAIR_LENGTH = 15;

  function layout(sequence, pairs, seedStart) {
    var n = sequence.length;
    var nodes = new Array(n);

    // Seed on a loose spiral: closer to the relaxed answer than a straight
    // line, so the layout converges in fewer iterations.
    for (var i = 0; i < n; i++) {
      var t = n > 1 ? i / (n - 1) : 0;
      var angle = t * Math.PI * 4;
      nodes[i] = {
        x: Math.cos(angle) * 30 + (Math.random() - 0.5) * 4,
        y: (t - 0.5) * BOND_LENGTH * n * 0.5,
        z: Math.sin(angle) * 30 + (Math.random() - 0.5) * 4,
        base: sequence[i],
        index: i,
        isSeed: i >= seedStart,
        paired: false
      };
    }

    var pairOf = {};
    for (var p = 0; p < pairs.length; p++) {
      pairOf[pairs[p][0]] = pairs[p][1];
      pairOf[pairs[p][1]] = pairs[p][0];
      nodes[pairs[p][0]].paired = true;
      nodes[pairs[p][1]].paired = true;
    }

    // Springs that hold a stem's two strands in register, giving ladders
    // instead of twisted rope.
    var stackBonds = [];
    for (var s = 0; s < pairs.length; s++) {
      var a = pairs[s][0];
      var b = pairs[s][1];
      if (pairOf[a + 1] === b - 1 && a + 1 < b - 1) {
        stackBonds.push([a, b - 1]);
        stackBonds.push([a + 1, b]);
      }
    }

    var iterations = n > 120 ? 220 : 340;

    for (var step = 0; step < iterations; step++) {
      var cooling = 1 - step / iterations;

      // repulsion (O(n^2); n is bounded by the 300 nt fold limit)
      for (var u = 0; u < n; u++) {
        for (var v = u + 1; v < n; v++) {
          var dx = nodes[v].x - nodes[u].x;
          var dy = nodes[v].y - nodes[u].y;
          var dz = nodes[v].z - nodes[u].z;
          var d2 = dx * dx + dy * dy + dz * dz;
          if (d2 > 3600 || d2 < 1e-6) continue;
          var d = Math.sqrt(d2);
          var force = (260 / d2) * cooling;
          var fx = (dx / d) * force;
          var fy = (dy / d) * force;
          var fz = (dz / d) * force;
          nodes[u].x -= fx;
          nodes[u].y -= fy;
          nodes[u].z -= fz;
          nodes[v].x += fx;
          nodes[v].y += fy;
          nodes[v].z += fz;
        }
      }

      applySprings(nodes, backboneBonds(n), BOND_LENGTH, 0.45);
      applySprings(nodes, pairs, PAIR_LENGTH, 0.55);
      applySprings(nodes, stackBonds, Math.sqrt(BOND_LENGTH * BOND_LENGTH + PAIR_LENGTH * PAIR_LENGTH), 0.25);
    }

    center(nodes);
    return { nodes: nodes, pairs: pairs.slice() };
  }

  function backboneBonds(n) {
    var bonds = [];
    for (var i = 0; i + 1 < n; i++) bonds.push([i, i + 1]);
    return bonds;
  }

  function applySprings(nodes, bonds, rest, k) {
    for (var i = 0; i < bonds.length; i++) {
      var a = nodes[bonds[i][0]];
      var b = nodes[bonds[i][1]];
      var dx = b.x - a.x;
      var dy = b.y - a.y;
      var dz = b.z - a.z;
      var d = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-6;
      var shift = ((d - rest) / d) * k * 0.5;
      var sx = dx * shift;
      var sy = dy * shift;
      var sz = dz * shift;
      a.x += sx;
      a.y += sy;
      a.z += sz;
      b.x -= sx;
      b.y -= sy;
      b.z -= sz;
    }
  }

  function center(nodes) {
    var cx = 0;
    var cy = 0;
    var cz = 0;
    for (var i = 0; i < nodes.length; i++) {
      cx += nodes[i].x;
      cy += nodes[i].y;
      cz += nodes[i].z;
    }
    cx /= nodes.length;
    cy /= nodes.length;
    cz /= nodes.length;
    for (var j = 0; j < nodes.length; j++) {
      nodes[j].x -= cx;
      nodes[j].y -= cy;
      nodes[j].z -= cz;
    }
  }

  function radius(nodes) {
    var max = 1;
    for (var i = 0; i < nodes.length; i++) {
      var d = Math.sqrt(nodes[i].x * nodes[i].x + nodes[i].y * nodes[i].y + nodes[i].z * nodes[i].z);
      if (d > max) max = d;
    }
    return max;
  }

  function create(canvas) {
    var ctx = canvas.getContext("2d");
    var model = null;
    var rotX = -0.25;
    var rotY = 0.4;
    var autoRotate = true;
    var dragging = false;
    var lastX = 0;
    var lastY = 0;
    var zoom = 1;
    var raf = null;
    var reduced = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

    function setModel(sequence, pairs, seedStart) {
      if (!sequence || !sequence.length) {
        model = null;
        draw();
        return;
      }
      model = layout(sequence, pairs, seedStart == null ? sequence.length : seedStart);
      model.radius = radius(model.nodes);
      draw();
    }

    function project(node, cx, cy, scale, dist) {
      var cosY = Math.cos(rotY);
      var sinY = Math.sin(rotY);
      var x = node.x * cosY - node.z * sinY;
      var z = node.x * sinY + node.z * cosY;
      var cosX = Math.cos(rotX);
      var sinX = Math.sin(rotX);
      var y = node.y * cosX - z * sinX;
      var z2 = node.y * sinX + z * cosX;
      var f = dist / (dist + z2);
      return { x: cx + x * f * scale, y: cy + y * f * scale, z: z2, f: f };
    }

    function draw() {
      var ratio = Math.min(window.devicePixelRatio || 1, 2);
      var rect = canvas.getBoundingClientRect();
      var w = Math.max(1, rect.width);
      var h = Math.max(1, rect.height);

      // Assigning canvas.width reallocates the backing store and is far too
      // expensive to do on every animation frame, so only touch it on a resize.
      var wantW = Math.floor(w * ratio);
      var wantH = Math.floor(h * ratio);
      if (canvas.width !== wantW || canvas.height !== wantH) {
        canvas.width = wantW;
        canvas.height = wantH;
      }
      ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
      ctx.clearRect(0, 0, w, h);

      if (!model) {
        ctx.fillStyle = "rgba(175, 188, 198, 0.6)";
        ctx.font = "13px Inter, system-ui, sans-serif";
        ctx.textAlign = "center";
        ctx.fillText("Type a strand to build the 3D model", w / 2, h / 2);
        return;
      }

      var cx = w / 2;
      var cy = h / 2;
      var dist = model.radius * 3.2;
      var scale = (Math.min(w, h) / (model.radius * 2.6)) * zoom;

      var pts = new Array(model.nodes.length);
      for (var i = 0; i < model.nodes.length; i++) {
        pts[i] = project(model.nodes[i], cx, cy, scale, dist);
      }

      // Collect every drawable, then paint far-to-near so depth reads correctly.
      var items = [];
      for (var b = 0; b + 1 < model.nodes.length; b++) {
        items.push({ kind: "backbone", a: b, b: b + 1, z: (pts[b].z + pts[b + 1].z) / 2 });
      }
      for (var p = 0; p < model.pairs.length; p++) {
        var pa = model.pairs[p][0];
        var pb = model.pairs[p][1];
        items.push({ kind: "pair", a: pa, b: pb, z: (pts[pa].z + pts[pb].z) / 2 });
      }
      for (var nIdx = 0; nIdx < model.nodes.length; nIdx++) {
        items.push({ kind: "base", a: nIdx, z: pts[nIdx].z });
      }
      items.sort(function (m, n) {
        return n.z - m.z;
      });

      for (var it = 0; it < items.length; it++) {
        var item = items[it];
        if (item.kind === "backbone") {
          var p1 = pts[item.a];
          var p2 = pts[item.b];
          var depth = clamp((p1.f + p2.f) / 2, 0.35, 1.4);
          ctx.strokeStyle = "rgba(132, 150, 163, " + (0.22 * depth).toFixed(3) + ")";
          ctx.lineWidth = 3.4 * depth;
          ctx.lineCap = "round";
          ctx.beginPath();
          ctx.moveTo(p1.x, p1.y);
          ctx.lineTo(p2.x, p2.y);
          ctx.stroke();
        } else if (item.kind === "pair") {
          var q1 = pts[item.a];
          var q2 = pts[item.b];
          var pd = clamp((q1.f + q2.f) / 2, 0.35, 1.4);
          ctx.strokeStyle = "rgba(227, 164, 79, " + (0.5 * pd).toFixed(3) + ")";
          ctx.lineWidth = 1.8 * pd;
          ctx.setLineDash([4 * pd, 3 * pd]);
          ctx.beginPath();
          ctx.moveTo(q1.x, q1.y);
          ctx.lineTo(q2.x, q2.y);
          ctx.stroke();
          ctx.setLineDash([]);
        } else {
          var node = model.nodes[item.a];
          var pt = pts[item.a];
          var d = clamp(pt.f, 0.3, 1.5);
          var r = 6.2 * d;
          var color = BASE_COLORS[node.base] || "#8496a3";

          var grad = ctx.createRadialGradient(pt.x - r * 0.35, pt.y - r * 0.35, r * 0.15, pt.x, pt.y, r);
          grad.addColorStop(0, "#e9eef2");
          grad.addColorStop(0.35, color);
          grad.addColorStop(1, shade(color, -0.45));
          ctx.fillStyle = grad;
          ctx.globalAlpha = clamp(d, 0.35, 1);
          ctx.beginPath();
          ctx.arc(pt.x, pt.y, r, 0, Math.PI * 2);
          ctx.fill();

          if (node.isSeed) {
            ctx.strokeStyle = "rgba(227, 164, 79, " + (0.9 * d).toFixed(3) + ")";
            ctx.lineWidth = 1.8 * d;
            ctx.beginPath();
            ctx.arc(pt.x, pt.y, r + 2.2, 0, Math.PI * 2);
            ctx.stroke();
          }
          ctx.globalAlpha = 1;

          if (r > 6 && model.nodes.length <= 80) {
            ctx.fillStyle = "rgba(9, 13, 17, 0.88)";
            ctx.font = "600 " + Math.round(8.5 * d) + "px Inter, system-ui, sans-serif";
            ctx.textAlign = "center";
            ctx.textBaseline = "middle";
            ctx.fillText(node.base, pt.x, pt.y + 0.5);
          }
        }
      }
    }

    function tick() {
      if (autoRotate && !dragging && !reduced) {
        rotY += 0.0045;
        draw();
      }
      raf = requestAnimationFrame(tick);
    }

    canvas.addEventListener("pointerdown", function (event) {
      dragging = true;
      autoRotate = false;
      lastX = event.clientX;
      lastY = event.clientY;
      canvas.setPointerCapture(event.pointerId);
      canvas.style.cursor = "grabbing";
    });

    canvas.addEventListener("pointermove", function (event) {
      if (!dragging) return;
      rotY += (event.clientX - lastX) * 0.008;
      rotX += (event.clientY - lastY) * 0.008;
      rotX = clamp(rotX, -1.4, 1.4);
      lastX = event.clientX;
      lastY = event.clientY;
      draw();
    });

    function endDrag() {
      dragging = false;
      canvas.style.cursor = "grab";
    }
    canvas.addEventListener("pointerup", endDrag);
    canvas.addEventListener("pointercancel", endDrag);

    canvas.addEventListener(
      "wheel",
      function (event) {
        event.preventDefault();
        zoom = clamp(zoom * (event.deltaY > 0 ? 0.92 : 1.08), 0.4, 4);
        draw();
      },
      { passive: false }
    );

    // The spin loop re-measures every frame, but a paused viewer would keep a
    // stale 1px backing store if it was first drawn before layout settled.
    if ("ResizeObserver" in window) {
      new ResizeObserver(function () {
        draw();
      }).observe(canvas);
    }

    tick();

    return {
      setModel: setModel,
      draw: draw,
      resetView: function () {
        rotX = -0.25;
        rotY = 0.4;
        zoom = 1;
        autoRotate = true;
        draw();
      },
      toggleSpin: function () {
        autoRotate = !autoRotate;
        return autoRotate;
      },
      isSpinning: function () {
        return autoRotate;
      },
      destroy: function () {
        if (raf) cancelAnimationFrame(raf);
      }
    };
  }

  function clamp(value, min, max) {
    return value < min ? min : value > max ? max : value;
  }

  function shade(hex, amount) {
    var num = parseInt(hex.slice(1), 16);
    var r = (num >> 16) & 255;
    var g = (num >> 8) & 255;
    var b = num & 255;
    r = Math.round(clamp(r + r * amount, 0, 255));
    g = Math.round(clamp(g + g * amount, 0, 255));
    b = Math.round(clamp(b + b * amount, 0, 255));
    return "rgb(" + r + ", " + g + ", " + b + ")";
  }

  return { create: create, BASE_COLORS: BASE_COLORS };
})();
