/*
 * Hero background: a rotating DNA double helix with a CRISPR guide RNA strand
 * bound across its middle, drawn live on a canvas.
 *
 * Rendered rather than photographed so it stays sharp at any viewport size,
 * ships no image bytes, and carries no licensing question. Build points in 3D,
 * spin about the helix axis, perspective-divide, sort back-to-front, draw with
 * depth-scaled alpha.
 *
 * Reading as "DNA" depends almost entirely on two things: a chunky
 * diameter-to-length ratio, and clearly visible base-pair rungs. A long thin
 * spiral with faint rungs just reads as a squiggle. Hence the proportions
 * below. They are deliberately stockier than real B-DNA.
 *
 * Respects prefers-reduced-motion (renders one static frame) and pauses when
 * the tab is hidden or the hero scrolls out of view.
 */

window.Hero = (function () {
  "use strict";

  var BASE_COLORS = {
    A: "#67c39a",
    U: "#e0a458",
    G: "#7fa9dc",
    C: "#d98ba0"
  };

  var STRAND_A = "#58c4d4";
  var STRAND_B = "#e3a44f";
  var GUIDE = "#e3a44f";

  function clamp(v, min, max) {
    return v < min ? min : v > max ? max : v;
  }

  function create(canvas) {
    var ctx = canvas.getContext("2d");
    var reduced = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    var phase = 0;
    var raf = null;
    var visible = true;
    var width = 0;
    var height = 0;

    // Helix axis runs along X, across the hero.
    var TURNS = 3.2;
    var SAMPLES = 200;
    var RADIUS = 88; // fat enough to read as a helix rather than a thread
    var LENGTH = 560;
    var RUNG_EVERY = 4;
    // Real strands are not diametrically opposed; the offset makes the groove.
    var MINOR_GROOVE = 2.3;

    // The guide RNA is confined to the middle of the duplex, like a guide bound
    // to its target site. Letting it run the full length just made a tangle.
    var GUIDE_FROM = 0.3;
    var GUIDE_TO = 0.78;

    var particles = [];
    for (var i = 0; i < 46; i++) {
      particles.push({
        x: (Math.random() - 0.5) * LENGTH * 1.4,
        y: (Math.random() - 0.5) * 320,
        z: (Math.random() - 0.5) * 280,
        r: Math.random() * 1.5 + 0.4
      });
    }

    // Fixed sequence so the rung colours are stable between frames.
    var bases = [];
    var alphabet = ["A", "U", "G", "C"];
    for (var s = 0; s < SAMPLES; s++) {
      bases.push(alphabet[Math.floor(Math.random() * 4)]);
    }

    function resize() {
      var ratio = Math.min(window.devicePixelRatio || 1, 2);
      var rect = canvas.getBoundingClientRect();
      width = Math.max(1, rect.width);
      height = Math.max(1, rect.height);
      canvas.width = Math.floor(width * ratio);
      canvas.height = Math.floor(height * ratio);
      ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    }

    function project(p, cx, cy, scale, dist, tilt) {
      // spin about the helix axis (X), then tilt the whole assembly
      var cosP = Math.cos(phase);
      var sinP = Math.sin(phase);
      var y = p.y * cosP - p.z * sinP;
      var z = p.y * sinP + p.z * cosP;

      var cosT = Math.cos(tilt);
      var sinT = Math.sin(tilt);
      var x2 = p.x * cosT - y * sinT;
      var y2 = p.x * sinT + y * cosT;

      var f = dist / (dist + z);
      return { x: cx + x2 * f * scale, y: cy + y2 * f * scale, z: z, f: f };
    }

    function helixPoint(t, offset, radius) {
      var angle = t * TURNS * Math.PI * 2 + offset;
      return {
        x: (t - 0.5) * LENGTH,
        y: Math.cos(angle) * radius,
        z: Math.sin(angle) * radius
      };
    }

    function draw() {
      ctx.clearRect(0, 0, width, height);

      var cx = width * 0.63;
      var cy = height * 0.5;
      var scale = clamp(Math.min(width / 1000, height / 560), 0.52, 1.2);
      var dist = 520;
      var tilt = -0.3;

      var items = [];

      for (var i = 0; i < particles.length; i++) {
        var pp = project(particles[i], cx, cy, scale, dist, tilt);
        items.push({ kind: "particle", p: pp, r: particles[i].r, z: pp.z });
      }

      var prevA = null;
      var prevB = null;
      var prevG = null;

      for (var s = 0; s < SAMPLES; s++) {
        var t = s / (SAMPLES - 1);
        var a = project(helixPoint(t, 0, RADIUS), cx, cy, scale, dist, tilt);
        var b = project(helixPoint(t, Math.PI + MINOR_GROOVE, RADIUS), cx, cy, scale, dist, tilt);

        if (prevA) items.push({ kind: "strand", from: prevA, to: a, color: STRAND_A, w: 3.4, z: (prevA.z + a.z) / 2 });
        if (prevB) items.push({ kind: "strand", from: prevB, to: b, color: STRAND_B, w: 3.4, z: (prevB.z + b.z) / 2 });

        // Base-pair rungs: the ladder is what makes the shape legible as DNA.
        if (s % RUNG_EVERY === 0) {
          items.push({ kind: "rung", from: a, to: b, base: bases[s], z: (a.z + b.z) / 2 });
        }

        if (t >= GUIDE_FROM && t <= GUIDE_TO) {
          var g = project(helixPoint(t, 1.5, RADIUS * 1.3), cx, cy, scale, dist, tilt);
          if (prevG) items.push({ kind: "guide", from: prevG, to: g, color: GUIDE, w: 2.2, z: (prevG.z + g.z) / 2 });
          prevG = g;
        } else {
          prevG = null;
        }

        prevA = a;
        prevB = b;
      }

      items.sort(function (m, n) {
        return n.z - m.z;
      });

      ctx.lineCap = "round";

      for (var k = 0; k < items.length; k++) {
        var item = items[k];
        var depth = clamp(item.p ? item.p.f : (item.from.f + item.to.f) / 2, 0.4, 1.6);
        var alpha = clamp((depth - 0.4) / 1.0, 0.06, 1);

        if (item.kind === "particle") {
          ctx.fillStyle = "rgba(175, 188, 198, " + (alpha * 0.45).toFixed(3) + ")";
          ctx.beginPath();
          ctx.arc(item.p.x, item.p.y, item.r * depth, 0, Math.PI * 2);
          ctx.fill();
          continue;
        }

        if (item.kind === "rung") {
          var color = BASE_COLORS[item.base];
          ctx.strokeStyle = hexToRgba(color, alpha * 0.72);
          ctx.lineWidth = 2.2 * depth;
          ctx.beginPath();
          ctx.moveTo(item.from.x, item.from.y);
          ctx.lineTo(item.to.x, item.to.y);
          ctx.stroke();

          // Nodes where the rung meets each backbone give the strands substance.
          ctx.fillStyle = hexToRgba(color, alpha * 0.9);
          ctx.beginPath();
          ctx.arc(item.from.x, item.from.y, 1.9 * depth, 0, Math.PI * 2);
          ctx.arc(item.to.x, item.to.y, 1.9 * depth, 0, Math.PI * 2);
          ctx.fill();
          continue;
        }

        // Glow via a wide translucent underlay plus a bright core stroke.
        // ctx.shadowBlur would look marginally softer but forces a full blur
        // pass per segment, which drops the hero to single-digit FPS.
        var isGuide = item.kind === "guide";
        var core = isGuide ? alpha * 0.5 : alpha * 0.95;
        var halo = isGuide ? alpha * 0.1 : alpha * 0.18;

        ctx.strokeStyle = hexToRgba(item.color, halo);
        ctx.lineWidth = item.w * depth * 4;
        ctx.beginPath();
        ctx.moveTo(item.from.x, item.from.y);
        ctx.lineTo(item.to.x, item.to.y);
        ctx.stroke();

        ctx.strokeStyle = hexToRgba(item.color, core);
        ctx.lineWidth = item.w * depth;
        ctx.beginPath();
        ctx.moveTo(item.from.x, item.from.y);
        ctx.lineTo(item.to.x, item.to.y);
        ctx.stroke();
      }
    }

    function tick() {
      if (visible) {
        phase += 0.0032;
        draw();
      }
      raf = requestAnimationFrame(tick);
    }

    resize();
    draw();

    if (!reduced) {
      tick();
      document.addEventListener("visibilitychange", function () {
        visible = !document.hidden;
      });
      if ("IntersectionObserver" in window) {
        new IntersectionObserver(
          function (entries) {
            visible = entries[0].isIntersecting && !document.hidden;
          },
          { threshold: 0 }
        ).observe(canvas);
      }
    }

    /*
     * Measure from a ResizeObserver rather than trusting the size at init.
     * At startup the canvas box can still be 0-wide (render-blocking webfont
     * CSS, 100svh not yet resolved), which sized the backing store to 1px and
     * left the hero blank until the user happened to resize the window.
     * The observer also fires immediately on observe, so it doubles as setup.
     */
    if ("ResizeObserver" in window) {
      new ResizeObserver(function () {
        resize();
        draw();
      }).observe(canvas);
    }

    var resizeTimer = null;
    window.addEventListener("resize", function () {
      // Backstop for devicePixelRatio changes, which leave the CSS box (and so
      // the ResizeObserver) untouched.
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(function () {
        resize();
        draw();
      }, 120);
    });

    return { draw: draw, resize: resize };
  }

  function hexToRgba(hex, alpha) {
    var num = parseInt(hex.slice(1), 16);
    return "rgba(" + ((num >> 16) & 255) + ", " + ((num >> 8) & 255) + ", " + (num & 255) + ", " + alpha.toFixed(3) + ")";
  }

  return { create: create };
})();
