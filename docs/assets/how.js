/* How Recore works: seven illustrated steps for one empty cell, driven by what
   the search finds, whether the neighbours agree, and the matrix policy. The
   decisions follow evidence_state and matrix_policy_corrects in src/pipeline.rs. */
(function () {
  "use strict";

  var stage = document.querySelector("[data-stage]");
  if (!stage) return;

  var NS = "http://www.w3.org/2000/svg";
  var reduceMotion = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  // `type` is the rescue type's name; `written` is how rescue_type reads in
  // orthology_evidence.tsv (tools/check_site.py holds both to the code).
  var SCENARIOS = [
    { id: "annotated", label: "The whole gene, annotated", type: "DROPPED", written: "DROPPED", identity: 99.6, coverage: 100 },
    { id: "unannotated", label: "The whole gene, nothing annotated", type: "MISSING_CDS", written: "MISSING_CDS", identity: 99.4, coverage: 100 },
    { id: "split", label: "The same gene, in another cluster", type: "SPLIT", written: "SPLIT:cluster_Y", identity: 99.2, coverage: 100 },
    { id: "stop", label: "The gene, with a stop codon", type: "PSEUDOGENE", written: "PSEUDOGENE:1stops", identity: 99.1, coverage: 100 },
    { id: "pieces", label: "The gene in two pieces", type: "FRAGMENTED", written: "FRAGMENTED:2frags:96%cov", identity: 99.0, coverage: 96 },
    { id: "paralog", label: "A related gene", type: "PARALOG", written: "PARALOG:cluster_Z", identity: 92.4, coverage: 97 },
    { id: "none", label: "Nothing", type: "", written: "", identity: 0, coverage: 0 }
  ];

  var CORRECTS = {
    strict: ["STRICT_PRESENT", "ALIAS_CLUSTER"],
    balanced: ["STRICT_PRESENT", "ALIAS_CLUSTER", "SEQUENCE_ONLY"],
    sensitive: ["STRICT_PRESENT", "ALIAS_CLUSTER", "SEQUENCE_ONLY", "PARTIAL_HOMOLOGY", "PSEUDOGENE"]
  };

  var state = { step: 0, scenario: "annotated", synteny: "agree", policy: "balanced" };
  var timers = [];
  var playTimer = null;

  function scenario() {
    return SCENARIOS.filter(function (s) {
      return s.id === state.scenario;
    })[0];
  }

  function evidenceState() {
    var sc = scenario();
    if (sc.id === "none") return "ABSENT";
    if (state.synteny === "differ") return "REJECTED_CONTEXT";
    var agree = state.synteny === "agree";
    switch (sc.id) {
      case "annotated":
        return agree ? "STRICT_PRESENT" : "SEQUENCE_ONLY";
      case "unannotated":
        return agree ? "STRICT_PRESENT" : "MISSING_CDS_UNVALIDATED";
      case "split":
        return "ALIAS_CLUSTER";
      case "stop":
        return "PSEUDOGENE";
      case "pieces":
        return "PARTIAL_HOMOLOGY";
      default:
        return agree ? "STRICT_PRESENT" : "HOMOLOG_PARALOG";
    }
  }

  function rejected() {
    var s = evidenceState();
    return s === "REJECTED_CONTEXT" || (s === "HOMOLOG_PARALOG" && state.synteny !== "agree");
  }

  function outcome() {
    var s = evidenceState();
    if (s === "ABSENT") return "absent";
    return CORRECTS[state.policy].indexOf(s) !== -1 ? "putback" : "kept";
  }

  function matrixAction() {
    if (outcome() === "putback") return "correct_presence";
    return rejected() ? "audit_only" : "audit_only_policy";
  }

  /* ── SVG helpers ── */

  function el(tag, attrs, parent) {
    var node = document.createElementNS(NS, tag);
    Object.keys(attrs || {}).forEach(function (k) {
      node.setAttribute(k, attrs[k]);
    });
    if (parent) parent.appendChild(node);
    return node;
  }

  function text(parent, x, y, content, attrs) {
    var a = attrs || {};
    a.x = x;
    a.y = y;
    var t = el("text", a, parent);
    t.textContent = content;
    return t;
  }

  // A gene drawn as an arrow pointing right.
  function gene(parent, x, y, w, cls, label) {
    var h = 30;
    var g = el("g", { class: "gene " + cls }, parent);
    el("path", { d: "M" + x + " " + y + "h" + (w - 12) + "l12 " + h / 2 + "l-12 " + h / 2 + "h-" + (w - 12) + "z" }, g);
    if (label) text(g, x + (w - 6) / 2, y + h / 2 + 5, label, { class: "gene-label", "text-anchor": "middle" });
    return g;
  }

  // Elements enter with a CSS animation after `delay` seconds.
  function enter(node, kind, delay) {
    node.classList.add("in-" + kind);
    node.style.animationDelay = (delay || 0) + "s";
    return node;
  }

  function later(ms, fn) {
    timers.push(setTimeout(fn, reduceMotion ? 0 : ms));
  }

  function animate(node, frames, options) {
    if (reduceMotion || !node.animate) {
      var last = frames[frames.length - 1];
      Object.keys(last).forEach(function (k) {
        if (k !== "offset" && k !== "easing") node.style[k] = last[k];
      });
      return;
    }
    node.animate(frames, Object.assign({ fill: "forwards", easing: "cubic-bezier(.4,0,.2,1)" }, options));
  }

  function countUp(node, target, suffix, ms) {
    if (reduceMotion) {
      node.textContent = target.toFixed(1) + suffix;
      return;
    }
    var start = performance.now();
    function tick(now) {
      var t = Math.min(1, (now - start) / ms);
      node.textContent = (target * (1 - Math.pow(1 - t, 3))).toFixed(1) + suffix;
      if (t < 1 && node.isConnected) requestAnimationFrame(tick);
    }
    requestAnimationFrame(tick);
  }

  // Rendered width of a text node, which grows with the phone font sizes.
  function measure(node, fallback) {
    try {
      var n = node.getComputedTextLength();
      return n > 0 ? n : fallback;
    } catch (e) {
      return fallback;
    }
  }

  function badge(parent, x, y, label, cls) {
    var g = el("g", { class: "badge-svg " + (cls || "") }, parent);
    var rect = el("rect", { y: y - 19, height: 38, rx: 19 }, g);
    var t = text(g, x, y + 6, label, { "text-anchor": "middle" });
    var w = Math.max(90, measure(t, label.length * 10.5) + 28);
    rect.setAttribute("x", x - w / 2);
    rect.setAttribute("width", w);
    return g;
  }

  /* ── The seven scenes ── */

  function sceneMatrix(svg) {
    var cols = ["A", "B", "X", "C", "D", "E", "F"];
    var x0 = 186;
    var y0 = 88;
    var cw = 58;
    var ch = 42;
    text(svg, x0 + (cols.length * cw) / 2, 40, "genes", { class: "axis", "text-anchor": "middle" });
    text(svg, 70, y0 + (5 * ch) / 2, "genomes", { class: "axis", "text-anchor": "middle", transform: "rotate(-90 70 " + (y0 + (5 * ch) / 2) + ")" });
    cols.forEach(function (c, i) {
      var t = text(svg, x0 + i * cw + cw / 2 - 3, y0 - 14, c === "X" ? "gene X" : c, { class: "col-label" + (c === "X" ? " is-x" : ""), "text-anchor": "middle" });
      enter(t, "fade", 0.1);
    });
    var empty = { "2,2": "target", "4,6": "gap", "1,5": "gap", "3,0": "gap" };
    for (var r = 0; r < 5; r++) {
      var label = text(svg, x0 - 18, y0 + r * ch + ch / 2 + 2, "G" + (r + 1), { class: "row-label" + (r === 2 ? " is-x" : ""), "text-anchor": "end" });
      enter(label, "fade", 0.1 + r * 0.06);
      for (var c = 0; c < cols.length; c++) {
        var kind = empty[r + "," + c];
        var cell = el("rect", {
          x: x0 + c * cw,
          y: y0 + r * ch,
          width: cw - 8,
          height: ch - 8,
          rx: 6,
          class: "m-cell " + (kind === "target" ? "m-target" : kind ? "m-empty" : "m-present")
        }, svg);
        enter(cell, "pop", 0.15 + r * 0.07 + c * 0.02);
      }
    }
    var cx = x0 + 2 * cw + (cw - 8) / 2;
    var cy = y0 + 2 * ch + (ch - 8) / 2;
    enter(text(svg, cx, cy + 7, "?", { class: "m-question", "text-anchor": "middle" }), "pop", 1.0);
    var note = text(svg, x0 + (cols.length * cw) / 2 - 4, y0 + 5 * ch + 34, "gene X in G3: absent?", { class: "callout-text", "text-anchor": "middle" });
    enter(note, "fade", 1.3);
  }

  function sceneSequence(svg) {
    var file = el("g", { class: "file" }, svg);
    el("path", { d: "M60 92h92l26 26v150h-118z" }, file);
    el("path", { class: "file-fold", d: "M152 92v26h26" }, file);
    for (var i = 0; i < 6; i++) {
      el("rect", { class: "file-line" + (i === 2 ? " is-x" : ""), x: 76, y: 134 + i * 18, width: i % 2 ? 66 : 84, height: 6, rx: 3 }, file);
    }
    text(svg, 119, 300, "pan_genome_reference.fa", { class: "small-label", "text-anchor": "middle" });
    enter(file, "fade", 0);

    var moving = el("g", {}, svg);
    gene(moving, 300, 150, 300, "gx", "gene X");
    animate(moving, [{ transform: "translate(-200px, 40px) scale(.5)", opacity: 0 }, { transform: "translate(0, 0) scale(1)", opacity: 1 }], { duration: 1100, delay: 300 });

    var bases = ["a", "c", "g", "t"];
    var strip = el("g", { class: "bases" }, svg);
    for (var k = 0; k < 36; k++) {
      var b = el("rect", { class: "base base-" + bases[(k * 7 + 3) % 4], x: 300 + k * 8.3, y: 206, width: 6.5, height: 22, rx: 2 }, strip);
      enter(b, "pop", 1.2 + k * 0.025);
    }
    enter(text(svg, 450, 262, "gene X, 1,281 bp", { class: "small-label", "text-anchor": "middle" }), "fade", 2.0);
    var arrow = el("path", { class: "flow-arrow", d: "M190 180 C 230 180, 250 170, 290 168", pathLength: 1 }, svg);
    enter(arrow, "draw", 0.4);
  }

  function sceneSearch(svg) {
    var sc = scenario();
    text(svg, 330, 34, "genome G3's assembly", { class: "axis", "text-anchor": "middle" });
    var cell = el("rect", { class: "bact", x: 120, y: 56, width: 420, height: 270, rx: 135 }, svg);
    enter(cell, "fade", 0);
    var cx = 330;
    var cy = 191;
    var r = 92;
    el("circle", { class: "chromosome", cx: cx, cy: cy, r: r }, svg);
    // Genes along the chromosome, as coloured arcs.
    var palette = ["ga", "gb", "gc", "gd", "gp", "ga", "gc", "gb", "gd", "gp", "gc", "ga"];
    palette.forEach(function (cls, i) {
      var a0 = (i / palette.length) * Math.PI * 2 - Math.PI / 2;
      var a1 = a0 + (Math.PI * 2) / palette.length - 0.14;
      el("path", {
        class: "arc " + cls,
        d: "M" + (cx + r * Math.cos(a0)) + " " + (cy + r * Math.sin(a0)) + " A " + r + " " + r + " 0 0 1 " + (cx + r * Math.cos(a1)) + " " + (cy + r * Math.sin(a1))
      }, svg);
    });
    var stopAt = sc.id === "none" ? 360 : sc.id === "paralog" ? 125 : 228;
    var rad = ((stopAt - 90) * Math.PI) / 180;
    var spot = { x: cx + r * Math.cos(rad), y: cy + r * Math.sin(rad) };

    var lens = el("g", { class: "lens" }, svg);
    el("line", { x1: cx + 18, y1: cy - r + 18, x2: cx + 40, y2: cy - r + 40 }, lens);
    el("circle", { cx: cx, cy: cy - r, r: 24 }, lens);
    lens.style.transformOrigin = cx + "px " + cy + "px";
    lens.style.transformBox = "view-box";
    animate(lens, [{ transform: "rotate(0deg)" }, { transform: "rotate(" + stopAt + "deg)" }], { duration: 2600, delay: 300, easing: "cubic-bezier(.45,.05,.3,1)" });

    var panel = el("g", { class: "panel" }, svg);
    el("rect", { x: 562, y: 110, width: 150, height: 150, rx: 14 }, panel);
    text(panel, 580, 140, "identity", { class: "small-label" });
    var idText = text(panel, 580, 170, "0.0%", { class: "big-num" });
    text(panel, 580, 204, "coverage", { class: "small-label" });
    var covText = text(panel, 580, 234, "0.0%", { class: "big-num" });
    enter(panel, "fade", 0.2);

    later(2950, function () {
      if (sc.id === "none") {
        idText.textContent = "< 90%";
        covText.textContent = "< 80%";
        enter(badge(svg, 330, 352, "no hit above the thresholds", "is-muted"), "pop", 0);
        return;
      }
      var glow = el("circle", { class: "glow", cx: spot.x, cy: spot.y, r: 18 }, svg);
      enter(glow, "pulse", 0);
      countUp(idText, sc.identity, "%", 900);
      countUp(covText, sc.coverage, "%", 900);
      enter(badge(svg, 330, 352, sc.id === "paralog" ? "a match, elsewhere" : "a match", "is-good"), "pop", 0.4);
    });
  }

  function sceneLook(svg) {
    var sc = scenario();
    var y = 222;
    el("line", { class: "contig", x1: 30, y1: y + 15, x2: 690, y2: y + 15 }, svg);
    var paralog = sc.id === "paralog";
    var left = paralog ? [["ge", "E"], ["gf", "F"]] : [["ga", "A"], ["gb", "B"]];
    var right = paralog ? [["gg", "G"], ["gh", "H"]] : [["gc", "C"], ["gd", "D"]];
    enter(gene(svg, 50, y, 96, left[0][0], left[0][1]), "fade", 0);
    enter(gene(svg, 158, y, 96, left[1][0], left[1][1]), "fade", 0.05);
    enter(gene(svg, 462, y, 96, right[0][0], right[0][1]), "fade", 0.1);
    enter(gene(svg, 570, y, 96, right[1][0], right[1][1]), "fade", 0.15);

    var lx = 280;
    var lw = 156;
    var label = "";
    if (sc.id === "annotated") {
      enter(gene(svg, lx, y, lw, "gp", "CDS, no cluster"), "fade", 0.2);
      label = "a CDS the pangenome tool left out";
    } else if (sc.id === "stop") {
      enter(gene(svg, lx, y, lw, "gp", ""), "fade", 0.2);
      label = "the gene is there, with a stop codon inside";
    } else if (sc.id === "unannotated") {
      var bracket = el("path", { class: "bracket", d: "M" + lx + " " + (y + 34) + "v10h" + lw + "v-10" }, svg);
      enter(bracket, "fade", 0.2);
      enter(text(svg, lx + lw / 2, y + 66, "no CDS annotated here", { class: "small-label", "text-anchor": "middle" }), "fade", 0.3);
      label = "the gene is there, the GFF missed it";
    } else if (sc.id === "split" || paralog) {
      enter(gene(svg, lx, y, lw, "go", paralog ? "cluster Z" : "cluster Y"), "fade", 0.2);
      label = paralog ? "a related gene of another family" : "the same gene, filed under another cluster";
    } else if (sc.id === "pieces") {
      enter(gene(svg, lx, y, 62, "gp", "1"), "fade", 0.2);
      var is = el("g", { class: "is-el" }, svg);
      el("rect", { x: lx + 64, y: y - 4, width: 28, height: 38, rx: 5 }, is);
      text(is, lx + 78, y + 21, "IS", { "text-anchor": "middle" });
      enter(is, "pop", 0.3);
      enter(gene(svg, lx + 94, y, 62, "gp", "2"), "fade", 0.25);
      label = "one gene, broken by an insertion";
    } else {
      label = "nothing there reaches 90% identity";
    }

    // The query drops onto what it matched.
    var q = el("g", {}, svg);
    gene(q, lx, 92, lw, "gx", "gene X");
    if (sc.id === "none") {
      animate(q, [{ transform: "translate(0,-60px)", opacity: 0 }, { transform: "translate(0,0)", opacity: 0.45 }], { duration: 700 });
      var cross = el("g", { class: "cross" }, svg);
      el("path", { d: "M" + (lx + lw / 2 - 14) + " 150 l28 28 M" + (lx + lw / 2 + 14) + " 150 l-28 28" }, cross);
      enter(cross, "pop", 0.8);
    } else {
      animate(q, [{ transform: "translate(0,-60px)", opacity: 0 }, { transform: "translate(0,0)", opacity: 1 }], { duration: 700 });
      // Each guide joins a point of the query to where it lands. The pieces
      // split the query in the middle, around the insertion.
      var mid = lx + lw / 2;
      var guides = sc.id === "pieces" ? [[lx, lx], [mid, lx + 62], [mid, lx + 94], [lx + lw, lx + lw]] : [[lx, lx], [lx + lw, lx + lw]];
      guides.forEach(function (pair, k) {
        enter(el("line", { class: "guide", x1: pair[0], y1: 126, x2: pair[1], y2: y - 4, pathLength: 1 }, svg), "draw", 0.6 + k * 0.08);
      });
      enter(text(svg, lx + lw + 14, 112, sc.identity.toFixed(1) + "% identical", { class: "small-label" }), "fade", 0.9);
      if (sc.id === "stop") {
        var stopSign = el("g", { class: "stop-sign" }, svg);
        var sx = lx + 112;
        el("polygon", { points: [sx - 6, y + 3, sx + 6, y + 3, sx + 12, y + 9, sx + 12, y + 21, sx + 6, y + 27, sx - 6, y + 27, sx - 12, y + 21, sx - 12, y + 9].join(" ") }, stopSign);
        text(stopSign, sx, y + 18, "stop", { "text-anchor": "middle" });
        enter(stopSign, "pop", 1.1);
        enter(el("rect", { class: "stop-tick", x: sx - 2, y: 92, width: 4, height: 30 }, svg), "pop", 1.1);
      }
    }
    enter(badge(svg, 360, 330, sc.id === "none" ? "no rescue" : sc.type, sc.id === "none" ? "is-muted" : ""), "pop", 1.2);
    enter(text(svg, 360, 372, label, { class: "small-label", "text-anchor": "middle" }), "fade", 1.3);
  }

  // The fading lives on an inner group: the entrance animation sets the
  // outer group's opacity and would undo it.
  function track(svg, y, genes, faint, dx) {
    var g = el("g", { class: "track" }, svg);
    var body = el("g", faint ? { class: "faint" } : {}, g);
    var shift = dx || 0;
    el("line", { class: "contig", x1: 60 + shift, y1: y + 15, x2: 660 + shift, y2: y + 15 }, body);
    var xs = [80, 190, 300, 450, 560];
    genes.forEach(function (spec, i) {
      gene(body, xs[i] + shift, y, i === 2 ? 140 : 96, spec[0], spec[1]);
    });
    return { g: g, xs: xs };
  }

  function sceneNeighbours(svg) {
    var sc = scenario();
    var none = sc.id === "none";
    var differ = state.synteny === "differ";
    var off = state.synteny === "off";
    var hitCls = sc.id === "split" || sc.id === "paralog" ? "go" : sc.id === "unannotated" ? "gx-ghost" : "gp";
    var top = differ ? [["ge", "E"], ["gf", "F"], [hitCls, "hit"], ["gg", "G"], ["gh", "H"]] : [["ga", "A"], ["gb", "B"], [hitCls, "hit"], ["gc", "C"], ["gd", "D"]];
    text(svg, 60, 326, "the genomes that have gene X", { class: "axis" });
    var faintAll = off || none;
    track(svg, 266, [["ga", ""], ["gb", ""], ["gx", ""], ["gc", ""], ["gd", ""]], true, 10);
    var b = track(svg, 256, [["ga", "A"], ["gb", "B"], ["gx", "gene X"], ["gc", "C"], ["gd", "D"]], faintAll);
    enter(b.g, "fade", 0.15);
    if (none) {
      // No hit in G3, so there is no place to look around.
      enter(badge(svg, 360, 150, "nothing found, nothing to compare", "is-muted"), "pop", 0.4);
      return;
    }
    text(svg, 60, 66, "G3, around the hit", { class: "axis" });
    var t = track(svg, 84, top, faintAll);
    enter(t.g, "fade", 0);
    if (off) {
      enter(badge(svg, 360, 190, "not compared without --synteny", "is-muted"), "pop", 0.4);
      return;
    }
    [0, 1, 3, 4].forEach(function (i, k) {
      var x = t.xs[i] + 44;
      var path = el("path", {
        class: differ ? "link link-bad" : "link link-good",
        d: "M" + x + " 120 C " + x + " 170, " + x + " 200, " + x + " 250",
        pathLength: 1
      }, svg);
      enter(path, "draw", 0.4 + k * 0.18);
    });
    enter(badge(svg, 360, 190, differ ? "different neighbours: another locus" : "same neighbours: same locus", differ ? "is-bad" : "is-good"), "pop", 1.2);
  }

  function sceneDecide(svg) {
    var s = evidenceState();
    var out = outcome();
    var chip = el("g", { class: "chip-svg" }, svg);
    el("rect", { x: 22, y: 168, width: 196, height: 44, rx: 12 }, chip);
    text(chip, 120, 196, s, { "text-anchor": "middle", class: s.length > 18 ? "long" : "" });
    enter(chip, "pop", 0);
    var gates = { strict: 300, balanced: 420, sensitive: 540 };
    Object.keys(gates).forEach(function (p) {
      var gx = gates[p];
      var on = p === state.policy;
      var allowed = CORRECTS[p].indexOf(s) !== -1;
      var g = el("g", { class: "gate" + (on ? " is-on" : "") + (allowed ? " allows" : "") }, svg);
      text(g, gx, 96, p, { "text-anchor": "middle" });
      el("rect", { class: "frame", x: gx - 40, y: 112, width: 80, height: 150, rx: 8 }, g);
      var door = el("rect", { class: "door", x: gx - 34, y: 118, width: 68, height: 138, rx: 5 }, g);
      enter(g, "fade", 0.1);
      if (on && out === "putback") {
        // The door rolls up inside its frame.
        door.style.transformBox = "fill-box";
        door.style.transformOrigin = "50% 0%";
        later(1300, function () {
          animate(door, [{ transform: "scaleY(1)" }, { transform: "scaleY(0.05)" }], { duration: 600 });
        });
      }
    });
    var target = gates[state.policy];
    var token = el("circle", { class: "token", cx: 218, cy: 190, r: 11 }, svg);
    token.style.opacity = 0;
    if (out === "absent") {
      enter(badge(svg, 420, 330, "nothing found: the cell stays absent", "is-muted"), "pop", 0.6);
      return;
    }
    var put = out === "putback";
    var result = el("g", { class: "result" }, svg);
    el("path", { class: "flow-arrow", d: "M" + target + " 268 V 294 M" + (target - 6) + " 287 L " + target + " 294 L " + (target + 6) + " 287", pathLength: 1 }, result);
    el("rect", { class: "r-cell " + (put ? "r-put" : "r-kept"), x: target - 24, y: 300, width: 48, height: 48, rx: 8 }, result);
    // The label goes on whichever side of the cell has room.
    var room = target < 500;
    var lines = put ? ["present", "in the corrected matrix"] : ["kept out,", "with the reason"];
    var anchor = { "text-anchor": room ? "start" : "end" };
    text(result, room ? target + 40 : target - 40, 318, lines[0], Object.assign({ class: "small-label result-head" }, anchor));
    text(result, room ? target + 40 : target - 40, 342, lines[1], Object.assign({ class: "small-label" }, anchor));
    result.style.opacity = 0;
    var toGate = target - 218;
    animate(token, [{ transform: "translate(0,0)", opacity: 1 }, { transform: "translate(" + toGate + "px, 0)", opacity: 1 }], { duration: 1000, delay: 300 });
    if (put) {
      later(1900, function () {
        animate(token, [{ transform: "translate(" + toGate + "px, 0)", opacity: 1 }, { transform: "translate(" + toGate + "px, 134px)", opacity: 0 }], { duration: 700 });
        animate(result, [{ opacity: 0 }, { opacity: 1 }], { duration: 500, delay: 400 });
      });
    } else {
      later(1400, function () {
        var gate = svg.querySelector(".gate.is-on .door");
        if (gate) animate(gate, [{ transform: "translateX(0)" }, { transform: "translateX(-5px)" }, { transform: "translateX(5px)" }, { transform: "translateX(0)" }], { duration: 400, fill: "none" });
        animate(token, [{ transform: "translate(" + toGate + "px, 0)", opacity: 1 }, { transform: "translate(" + toGate + "px, 134px)", opacity: 0 }], { duration: 700, delay: 350 });
        animate(result, [{ opacity: 0 }, { opacity: 1 }], { duration: 500, delay: 800 });
      });
    }
  }

  function doc(svg, x, y, w, h, title) {
    var g = el("g", { class: "doc" }, svg);
    el("rect", { class: "doc-bg", x: x, y: y, width: w, height: h, rx: 12 }, g);
    text(g, x + 14, y + 26, title, { class: "doc-title" });
    return g;
  }

  function sceneWrite(svg) {
    var s = evidenceState();
    var out = outcome();
    var sc = scenario();
    var d1 = doc(svg, 20, 60, 214, 220, "corrected matrix");
    var d2 = doc(svg, 252, 60, 214, 220, "states matrix");
    var d3 = doc(svg, 484, 60, 216, 220, "orthology_evidence.tsv");
    [d1, d2].forEach(function (d, k) {
      var x0 = (k ? 252 : 20) + 28;
      for (var r = 0; r < 4; r++) {
        for (var c = 0; c < 5; c++) {
          var isTarget = r === 2 && c === 2;
          var cls = "m-cell m-present";
          if (isTarget) cls = k === 0 ? (out === "putback" ? "m-cell m-put" : "m-cell m-empty") : "m-cell m-state";
          if (!isTarget && (r * 5 + c) % 7 === 3) cls = "m-cell m-empty";
          el("rect", { class: cls, x: x0 + c * 32, y: 92 + r * 32, width: 26, height: 26, rx: 5 }, d);
        }
      }
    });
    text(d1, 127, 252, out === "putback" ? "G3, gene X: present" : "G3, gene X: absent", { class: "small-label", "text-anchor": "middle" });
    text(d2, 359, 252, s, { class: "small-label state-label" + (s.length > 18 ? " long" : ""), "text-anchor": "middle" });
    if (s === "ABSENT") {
      text(d3, 592, 140, "no row:", { class: "small-label", "text-anchor": "middle" });
      text(d3, 592, 164, "nothing was found", { class: "small-label", "text-anchor": "middle" });
    } else {
      var rows = [
        ["gene_cluster", "gene_X"],
        ["genome", "G3"],
        ["rescue_type", sc.written],
        ["evidence_state", s],
        ["matrix_action", matrixAction()]
      ];
      rows.forEach(function (row, k) {
        text(d3, 498, 104 + k * 30, row[0], { class: "tsv-key" });
        text(d3, 498, 119 + k * 30, row[1], { class: "tsv-val" });
      });
    }
    enter(d1, "rise", 0);
    enter(d2, "rise", 0.25);
    enter(d3, "rise", 0.5);
    var hl = el("rect", { class: "hl", x: 112, y: 152, width: 30, height: 30, rx: 7 }, svg);
    enter(hl, "pulse", 1.0);
    enter(badge(svg, 360, 330, "nothing is thrown away", "is-good"), "pop", 1.2);
  }

  var STEPS = [
    {
      title: "The matrix says absent",
      kicker: "Step 1",
      scene: sceneMatrix,
      text: function () {
        return "Your pangenome tool marked gene X absent in genome G3. Every empty cell like this one is a question Recore asks the assembly.";
      }
    },
    {
      title: "Take the gene's sequence",
      kicker: "Step 2",
      scene: sceneSequence,
      text: function () {
        return "Recore takes gene X's representative sequence from the pan-genome reference. With <code>--multi-allele</code> it takes every allele, so a diverged copy is not missed.";
      }
    },
    {
      title: "Search the assembly",
      kicker: "Step 3",
      scene: sceneSearch,
      text: function () {
        var base = "minimap2 searches G3's assembly. A hit counts if it is at least 90% identical and covers 80% of the gene (<code>--min-identity</code>, <code>--min-coverage</code>).";
        return scenario().id === "none" ? base + " Here nothing reaches them: the cell was searched, and the gene is not there." : base;
      }
    },
    {
      title: "Look at what is there",
      kicker: "Step 4",
      scene: sceneLook,
      text: function () {
        return {
          annotated: "The hit sits on a CDS the annotation has, but that the pangenome tool left out of every cluster. Recore calls it <b>DROPPED</b>.",
          unannotated: "The whole gene is in the assembly, but the GFF has no CDS there. Recore calls it <b>MISSING_CDS</b>.",
          split: "The hit sits on a gene of another cluster, at 99% identity: the same gene, split between two clusters. Recore calls it <b>SPLIT</b>.",
          stop: "The gene is there, but a stop codon interrupts it: a pseudogene. Recore reads the CDS in its own frame and counts the stops.",
          pieces: "The gene is in two pieces, here broken by an insertion sequence, but the pieces cover it in order. Recore calls it <b>FRAGMENTED</b>.",
          paralog: "The hit is a related gene of another family, 92% identical, at another place in the genome. Recore calls it <b>PARALOG</b>: a lookalike, not the missing gene.",
          none: "Nothing in G3 is close enough to gene X. That is a real absence, and the states file says it was searched."
        }[scenario().id];
      }
    },
    {
      title: "Check the neighbours",
      kicker: "Step 5",
      scene: sceneNeighbours,
      control: "synteny",
      text: function () {
        var base = "With <code>--synteny</code>, Recore compares the genes around the hit with the genes around gene X in the genomes that have it.";
        if (scenario().id === "none") return base + " Here there is no hit to compare.";
        if (state.synteny === "off") return base + " Without it, the neighbours are not compared, and the evidence has to stand on the sequence alone.";
        if (state.synteny === "differ") return base + " Different neighbours mean a copy somewhere else: the candidate is rejected and kept as audit evidence.";
        return base + " The same neighbours mean the same place in the genome: strong evidence it is gene X.";
      }
    },
    {
      title: "Decide",
      kicker: "Step 6",
      scene: sceneDecide,
      control: "policy",
      text: function () {
        var s = evidenceState();
        var out = outcome();
        var p = "<code>" + state.policy + "</code>";
        var sentence = "The evidence is <b>" + s + "</b>. <code>--matrix-policy</code> decides which evidence becomes presence in the corrected matrix.";
        if (out === "absent") return sentence + " Nothing was found, so the cell stays absent.";
        if (out === "putback") {
          var guard = state.policy === "sensitive" && (s === "PSEUDOGENE" || s === "PARTIAL_HOMOLOGY") ? ", if it also passes the quality guard" : "";
          return sentence + " Under " + p + " this one goes in" + guard + ": the cell becomes present.";
        }
        var anyPolicy = Object.keys(CORRECTS).some(function (k) {
          return CORRECTS[k].indexOf(s) !== -1;
        });
        if (!anyPolicy) return sentence + " No policy takes <b>" + s + "</b>: under all three it stays out, written down with the reason.";
        return sentence + " Under " + p + " it stays out, written down with the reason.";
      }
    },
    {
      title: "Write it down",
      kicker: "Step 7",
      scene: sceneWrite,
      text: function () {
        return "Every decision is written with its evidence: the corrected matrix for your next tool, a state for every cell, and one row per decision in <code>orthology_evidence.tsv</code>. Kept-out evidence is not lost.";
      }
    }
  ];

  /* ── Page ── */

  var stepsList = document.querySelector("[data-steps]");
  var kicker = document.querySelector("[data-kicker]");
  var title = document.querySelector("[data-step-title]");
  var body = document.querySelector("[data-step-text]");
  var control = document.querySelector("[data-control]");
  var progress = document.querySelector("[data-progress]");
  var prevButton = document.querySelector("[data-prev]");
  var nextButton = document.querySelector("[data-next]");
  var playButton = document.querySelector("[data-play]");
  var replayButton = document.querySelector("[data-replay]");
  var scenarioBox = document.querySelector("[data-scenarios]");
  var svgTitle = document.getElementById("stage-title");
  var svgDesc = document.getElementById("stage-desc");

  SCENARIOS.forEach(function (sc) {
    var b = document.createElement("button");
    b.type = "button";
    b.textContent = sc.label;
    b.setAttribute("aria-pressed", sc.id === state.scenario ? "true" : "false");
    b.addEventListener("click", function () {
      stopPlay();
      state.scenario = sc.id;
      // A lookalike from elsewhere has, by definition, other neighbours.
      state.synteny = sc.id === "paralog" ? "differ" : "agree";
      scenarioBox.querySelectorAll("button").forEach(function (other) {
        other.setAttribute("aria-pressed", other === b ? "true" : "false");
      });
      render();
    });
    scenarioBox.appendChild(b);
  });

  STEPS.forEach(function (step, i) {
    var li = document.createElement("li");
    var b = document.createElement("button");
    b.type = "button";
    b.innerHTML = '<span class="n">' + (i + 1) + "</span><span>" + step.title + "</span>";
    b.addEventListener("click", function () {
      stopPlay();
      go(i);
    });
    li.appendChild(b);
    stepsList.appendChild(li);
  });

  function segmented(options, current, onPick, label) {
    var wrap = document.createElement("div");
    wrap.className = "seg";
    wrap.setAttribute("role", "group");
    wrap.setAttribute("aria-label", label);
    options.forEach(function (opt) {
      var b = document.createElement("button");
      b.type = "button";
      b.innerHTML = opt[1];
      b.setAttribute("aria-pressed", opt[0] === current ? "true" : "false");
      b.addEventListener("click", function () {
        stopPlay();
        onPick(opt[0]);
        render();
      });
      wrap.appendChild(b);
    });
    return wrap;
  }

  function render() {
    var step = STEPS[state.step];
    Array.prototype.forEach.call(stage.querySelectorAll(":scope > :not(title):not(desc)"), function (n) {
      n.remove();
    });
    timers.forEach(clearTimeout);
    timers = [];
    var scene = el("g", { class: "scene" }, stage);
    step.scene(scene);

    kicker.textContent = step.kicker + " of " + STEPS.length;
    title.textContent = step.title;
    body.innerHTML = step.text();
    svgTitle.textContent = step.kicker + ": " + step.title;
    svgDesc.textContent = body.textContent;
    control.innerHTML = "";
    if (step.control === "synteny" && scenario().id !== "none") {
      control.appendChild(
        segmented(
          [["agree", "Neighbours agree"], ["differ", "Neighbours differ"], ["off", "Not checked (no <code>--synteny</code>)"]],
          state.synteny,
          function (v) {
            state.synteny = v;
          },
          "Neighbours"
        )
      );
    }
    if (step.control === "policy" && outcome() !== "absent") {
      control.appendChild(
        segmented(
          [["strict", "strict"], ["balanced", "balanced"], ["sensitive", "sensitive"]],
          state.policy,
          function (v) {
            state.policy = v;
          },
          "Matrix policy"
        )
      );
    }
    stepsList.querySelectorAll("button").forEach(function (b, i) {
      if (i === state.step) b.setAttribute("aria-current", "step");
      else b.removeAttribute("aria-current");
      b.classList.toggle("done", i < state.step);
    });
    progress.style.width = ((state.step + 1) / STEPS.length) * 100 + "%";
    prevButton.disabled = state.step === 0;
    nextButton.textContent = state.step === STEPS.length - 1 ? "Start again" : "Next";
  }

  function go(i) {
    state.step = (i + STEPS.length) % STEPS.length;
    render();
  }

  function stopPlay() {
    if (playTimer) clearInterval(playTimer);
    playTimer = null;
    playButton.textContent = "Play all";
  }

  prevButton.addEventListener("click", function () {
    stopPlay();
    go(state.step - 1);
  });
  nextButton.addEventListener("click", function () {
    stopPlay();
    go(state.step + 1);
  });
  replayButton.addEventListener("click", function () {
    stopPlay();
    render();
  });
  playButton.addEventListener("click", function () {
    if (playTimer) {
      stopPlay();
      return;
    }
    playButton.textContent = "Pause";
    go(0);
    playTimer = setInterval(function () {
      if (state.step === STEPS.length - 1) {
        stopPlay();
        return;
      }
      go(state.step + 1);
    }, 5200);
  });
  document.addEventListener("keydown", function (event) {
    if (/input|textarea|select/i.test(document.activeElement.tagName)) return;
    if (event.key === "ArrowRight") {
      stopPlay();
      go(state.step + 1);
    } else if (event.key === "ArrowLeft") {
      stopPlay();
      go(state.step - 1);
    }
  });

  // Labels and pills are sized for the font of the moment, so redraw when
  // the phone sizes switch on or off.
  if (window.matchMedia) {
    var phone = window.matchMedia("(max-width: 600px)");
    if (phone.addEventListener) phone.addEventListener("change", render);
  }

  render();
})();
