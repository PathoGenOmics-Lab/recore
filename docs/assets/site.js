/* Recore website: theme, tabs, copy buttons, the command builder, the
   animated matrix and the reference filter. Plain JavaScript, no build step. */
(function () {
  "use strict";

  var reduceMotion = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  /* ── Theme ── */

  function systemDark() {
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches;
  }

  function currentTheme() {
    var set = document.documentElement.dataset.theme;
    if (set === "light" || set === "dark") return set;
    return systemDark() ? "dark" : "light";
  }

  document.querySelectorAll("[data-theme-toggle]").forEach(function (button) {
    button.addEventListener("click", function () {
      var next = currentTheme() === "dark" ? "light" : "dark";
      document.documentElement.dataset.theme = next;
      try {
        localStorage.setItem("recore-theme", next);
      } catch (e) {
        /* private mode: the choice lasts for this page only */
      }
    });
  });

  /* ── Tabs ── */

  document.querySelectorAll('[role="tablist"]').forEach(function (list) {
    var tabs = Array.prototype.slice.call(list.querySelectorAll('[role="tab"]'));
    function select(tab) {
      tabs.forEach(function (t) {
        var on = t === tab;
        t.setAttribute("aria-selected", on ? "true" : "false");
        t.tabIndex = on ? 0 : -1;
        var panel = document.getElementById(t.getAttribute("aria-controls"));
        if (panel) panel.hidden = !on;
      });
    }
    tabs.forEach(function (tab, i) {
      tab.tabIndex = tab.getAttribute("aria-selected") === "true" ? 0 : -1;
      tab.addEventListener("click", function () {
        select(tab);
      });
      tab.addEventListener("keydown", function (event) {
        var step = event.key === "ArrowRight" ? 1 : event.key === "ArrowLeft" ? -1 : 0;
        if (!step) return;
        event.preventDefault();
        var next = tabs[(i + step + tabs.length) % tabs.length];
        select(next);
        next.focus();
      });
    });
  });

  /* ── Copy buttons ── */

  function plainText(node) {
    var clone = node.cloneNode(true);
    clone.querySelectorAll(".p").forEach(function (prompt) {
      prompt.remove();
    });
    return clone.textContent.replace(/\s+$/, "") + "\n";
  }

  function copyText(text) {
    if (navigator.clipboard && window.isSecureContext) {
      return navigator.clipboard.writeText(text);
    }
    return new Promise(function (resolve, reject) {
      var area = document.createElement("textarea");
      area.value = text;
      area.setAttribute("readonly", "");
      area.style.position = "fixed";
      area.style.opacity = "0";
      document.body.appendChild(area);
      area.select();
      try {
        document.execCommand("copy") ? resolve() : reject();
      } catch (e) {
        reject(e);
      }
      area.remove();
    });
  }

  function addCopy(button, source) {
    button.addEventListener("click", function () {
      copyText(plainText(source)).then(
        function () {
          button.textContent = "Copied";
          button.classList.add("done");
          setTimeout(function () {
            button.textContent = "Copy";
            button.classList.remove("done");
          }, 1600);
        },
        function () {
          button.textContent = "Select and copy";
        }
      );
    });
  }

  document.querySelectorAll("[data-copy-from]").forEach(function (button) {
    var source = document.querySelector(button.getAttribute("data-copy-from"));
    if (source) addCopy(button, source);
  });

  // Every other code block gets a button too.
  document.querySelectorAll(".code").forEach(function (block) {
    if (block.querySelector(".copy")) return;
    var pre = block.querySelector("pre");
    if (!pre) return;
    var button = document.createElement("button");
    button.type = "button";
    button.className = "copy";
    button.textContent = "Copy";
    block.appendChild(button);
    addCopy(button, pre);
  });

  /* ── Command builder ── */

  var TOOLS = [
    { id: "panaroo", label: "Panaroo", p: "gene_presence_absence.csv", r: "pan_genome_reference.fa", note: "" },
    { id: "roary", label: "Roary", p: "gene_presence_absence.csv", r: "pan_genome_reference.fa", note: "" },
    {
      id: "ppanggolin",
      label: "PPanGGOLiN",
      p: "matrix.csv",
      r: "all_nucleotide_families.fasta",
      note: "matrix.csv comes from ppanggolin all (or write_pangenome --csv), the FASTA from ppanggolin fasta --gene_families all."
    },
    {
      id: "panta",
      label: "PanTA",
      p: "gene_presence_absence.csv",
      r: "representative_clusters_nucl.fasta",
      note: ""
    },
    {
      id: "pirate",
      label: "PIRATE",
      p: "PIRATE.gene_families.tsv",
      r: "representative_sequences.ffn",
      note: "Give Recore PIRATE's own table, not the output of PIRATE_to_roary.pl."
    }
  ];

  var builder = document.querySelector("[data-command]");
  if (builder) {
    var toolBox = document.querySelector("[data-tools]");
    var toolNote = document.querySelector("[data-tool-note]");
    var explain = document.querySelector("[data-explain]");
    var layout = document.querySelector("[data-layout]");
    var state = { tool: "panaroo" };

    TOOLS.forEach(function (tool) {
      var button = document.createElement("button");
      button.type = "button";
      button.textContent = tool.label;
      button.setAttribute("aria-pressed", tool.id === state.tool ? "true" : "false");
      button.addEventListener("click", function () {
        state.tool = tool.id;
        toolBox.querySelectorAll("button").forEach(function (b) {
          b.setAttribute("aria-pressed", b === button ? "true" : "false");
        });
        render();
      });
      toolBox.appendChild(button);
    });

    document.querySelectorAll("[data-opt], [data-opt-value]").forEach(function (input) {
      input.addEventListener("input", render);
      input.addEventListener("change", render);
    });

    function opt(name) {
      var box = document.querySelector('[data-opt="' + name + '"]');
      return !!(box && box.checked);
    }

    function value(name, fallback) {
      var input = document.querySelector('[data-opt-value="' + name + '"]');
      var v = input ? input.value.trim() : "";
      return v || fallback;
    }

    function esc(text) {
      return text.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
    }

    function shellWord(text) {
      return /^[A-Za-z0-9_.\/:@%+=,-]+$/.test(text) ? text : "'" + text.replace(/'/g, "'\\''") + "'";
    }

    function render() {
      var tool = TOOLS.filter(function (t) {
        return t.id === state.tool;
      })[0];
      var threads = String(Math.max(1, Math.min(256, parseInt(value("threads", "8"), 10) || 8)));
      var ref = shellWord(value("ref", "H37Rv"));
      var lines = [];
      var notes = [];
      function line(parts) {
        lines.push(
          parts
            .map(function (part) {
              return part.charAt(0) === "-"
                ? '<span class="f">' + esc(part) + "</span>"
                : '<span class="v">' + esc(part) + "</span>";
            })
            .join(" ")
        );
      }
      line(["-p", tool.p]);
      notes.push(["-p  -r", "the matrix and the representative sequences " + tool.label + " wrote"]);
      line(["-r", tool.r]);
      line(["-g", "genomes/"]);
      notes.push(["-g genomes/", "one assembly per genome, named like the matrix columns"]);
      if (opt("gff")) {
        line(["-a", "annotations/", "--synteny"]);
        notes.push(["-a  --synteny", "the GFFs " + tool.label + " was run on; each gene's neighbours are checked"]);
      }
      if (opt("tree")) {
        line(["--tree", "species.nwk"]);
        notes.push(["--tree", "a rooted Newick tree, tips named like the matrix columns"]);
      }
      if (opt("ref")) {
        line(["--reference-genome", ref]);
        notes.push(["--reference-genome", "a column of the matrix, or its AS identifier"]);
      }
      if (opt("report")) {
        line(["--write-html-report"]);
        notes.push(["--write-html-report", "writes recore_report.html"]);
      }
      line(["-t", threads, "-o", "recore_output/"]);
      notes.push(["-t  -o", "threads, and the folder the results go to"]);

      builder.innerHTML =
        "recore \\\n" +
        lines
          .map(function (l, i) {
            return "  " + l + (i < lines.length - 1 ? " \\" : "");
          })
          .join("\n");

      explain.innerHTML = notes
        .map(function (n) {
          return "<li><code>" + esc(n[0]) + "</code><span>" + esc(n[1]) + "</span></li>";
        })
        .join("");

      toolNote.textContent = tool.note;

      var rows = [
        ["your_project/", ""],
        ["├── " + tool.p, "-p"],
        ["├── " + tool.r, "-r"],
        ["├── genomes/", "-g"],
        ["│   ├── L1_AS00101.fasta", "one per matrix column"],
        ["│   └── …", ""]
      ];
      if (opt("gff")) {
        rows.push(["├── annotations/", "-a"]);
        rows.push(["│   ├── L1_AS00101.gff", "same names"]);
        rows.push(["│   └── …", ""]);
      }
      if (opt("tree")) rows.push(["├── species.nwk", "--tree"]);
      var last = rows.length - 1;
      rows[last][0] = rows[last][0].replace(/^├/, "└");
      var width = rows.reduce(function (w, r) {
        return Math.max(w, r[0].length);
      }, 0);
      layout.innerHTML = rows
        .map(function (r) {
          var pad = new Array(width - r[0].length + 4).join(" ");
          var name = esc(r[0]).replace(/(L1_AS00101\.\w+)/, '<span class="hl">$1</span>');
          return name + (r[1] ? pad + '<span class="cm"># ' + esc(r[1]) + "</span>" : "");
        })
        .join("\n");
    }

    render();
  }

  /* ── The animated matrix ── */

  var grid = document.querySelector("[data-demo-grid]");
  if (grid) {
    var COLS = 16;
    var ROWS = 10;
    // Shaped like the example run: two genes falsely absent from three
    // genomes and one from two (put back), a pseudogene across one lineage
    // (kept out), a gene really missing from eight genomes and a few
    // scattered true absences.
    var OUTCOME = {};
    function mark(kind, col, rows) {
      rows.forEach(function (r) {
        OUTCOME[r * COLS + col] = kind;
      });
    }
    mark("rescued", 3, [1, 4, 7]);
    mark("rescued", 8, [2, 5, 6]);
    mark("rescued", 12, [0, 9]);
    mark("kept", 5, [3, 4, 5]);
    mark("absent", 14, [0, 1, 2, 3, 4, 5, 6, 7]);
    mark("absent", 10, [2]);
    mark("absent", 1, [6]);
    mark("absent", 6, [9]);
    mark("absent", 15, [4]);

    grid.style.setProperty("--cols", COLS);
    var cells = [];
    for (var i = 0; i < ROWS * COLS; i++) {
      var cell = document.createElement("span");
      cell.className = "cell";
      grid.appendChild(cell);
      cells.push(cell);
    }

    var stepLabel = document.querySelector("[data-demo-step]");
    var tallies = {
      rescued: document.querySelector('[data-tally="rescued"]'),
      kept: document.querySelector('[data-tally="kept"]'),
      absent: document.querySelector('[data-tally="absent"]')
    };
    var counts;
    var generation = 0;
    var visible = true;

    function setTally() {
      Object.keys(tallies).forEach(function (k) {
        if (tallies[k]) tallies[k].textContent = String(counts[k]);
      });
    }

    function reset() {
      counts = { rescued: 0, kept: 0, absent: 0 };
      cells.forEach(function (c, idx) {
        c.className = OUTCOME[idx] ? "cell gap" : "cell";
      });
      setTally();
      if (stepLabel) stepLabel.textContent = "Input matrix";
    }

    function finalState() {
      counts = { rescued: 0, kept: 0, absent: 0 };
      cells.forEach(function (c, idx) {
        var kind = OUTCOME[idx];
        c.className = !kind ? "cell" : kind === "absent" ? "cell gap" : "cell " + kind;
        if (kind) counts[kind] += 1;
      });
      setTally();
      if (stepLabel) stepLabel.textContent = "After Recore";
    }

    function wait(ms, gen) {
      return new Promise(function (resolve) {
        function tick() {
          if (gen !== generation) return;
          if (!visible || document.hidden) {
            setTimeout(tick, 250);
            return;
          }
          resolve();
        }
        setTimeout(tick, ms);
      });
    }

    function run() {
      var gen = ++generation;
      reset();
      var chain = wait(1400, gen).then(function () {
        if (stepLabel) stepLabel.textContent = "Searching the assemblies";
      });
      for (var col = 0; col < COLS; col++) {
        (function (col) {
          chain = chain.then(function () {
            var gaps = [];
            for (var r = 0; r < ROWS; r++) {
              var idx = r * COLS + col;
              if (OUTCOME[idx]) gaps.push(idx);
            }
            gaps.forEach(function (idx) {
              cells[idx].classList.add("scan");
            });
            return wait(gaps.length ? 360 : 90, gen).then(function () {
              gaps.forEach(function (idx) {
                var kind = OUTCOME[idx];
                cells[idx].classList.remove("scan");
                if (kind !== "absent") {
                  cells[idx].classList.remove("gap");
                  cells[idx].classList.add(kind);
                }
                counts[kind] += 1;
              });
              setTally();
            });
          });
        })(col);
      }
      chain
        .then(function () {
          if (stepLabel) stepLabel.textContent = "Done";
          return wait(4200, gen);
        })
        .then(function () {
          if (gen === generation) run();
        });
    }

    if (reduceMotion) {
      finalState();
    } else {
      if ("IntersectionObserver" in window) {
        new IntersectionObserver(function (entries) {
          visible = entries[0].isIntersecting;
        }).observe(grid);
      }
      run();
    }
  }

  /* ── Embedded report ── */

  // The report takes the mouse wheel for zooming, so it only takes the mouse
  // after a click, and gives it back when the pointer leaves: otherwise the
  // page stops scrolling wherever the report is under the pointer.
  document.querySelectorAll("[data-report-stage]").forEach(function (stage) {
    var cover = stage.querySelector(".report-cover");
    var frame = stage.querySelector("iframe");
    if (!cover || !frame) return;
    cover.addEventListener("click", function () {
      stage.classList.add("live");
      try {
        frame.focus();
      } catch (e) {
        /* nothing to focus yet */
      }
    });
    stage.addEventListener("mouseleave", function () {
      stage.classList.remove("live");
    });
  });

  /* ── Reference filter ── */

  document.querySelectorAll("[data-filter]").forEach(function (input) {
    var table = document.querySelector(input.getAttribute("data-filter"));
    if (!table) return;
    var empty = document.querySelector(input.getAttribute("data-filter-empty"));
    input.addEventListener("input", function () {
      var query = input.value.trim().toLowerCase();
      var rows = Array.prototype.slice.call(table.tBodies[0].rows);
      var shown = 0;
      var group = null;
      var groupHasRows = false;
      function closeGroup() {
        if (group) group.hidden = !groupHasRows;
      }
      rows.forEach(function (row) {
        if (row.classList.contains("group-row")) {
          closeGroup();
          group = row;
          groupHasRows = false;
          return;
        }
        var match = !query || row.textContent.toLowerCase().indexOf(query) !== -1;
        row.hidden = !match;
        if (match) {
          shown += 1;
          groupHasRows = true;
        }
      });
      closeGroup();
      if (empty) empty.style.display = shown ? "none" : "block";
    });
  });
})();
