(() => {
  const root = document.documentElement;
  // The head bootstrap already resolved the theme before first paint; only sync
  // the toggle label and notify listeners if it somehow did not run.
  applyTheme(root.getAttribute("data-theme") === "dark" ? "dark" : "light");

  function applyTheme(mode) {
    root.setAttribute("data-theme", mode);
    const button = document.getElementById("themeToggle");
    if (button) {
      const isDark = mode === "dark";
      button.setAttribute("aria-pressed", String(isDark));
      button.querySelector(".label").textContent = isDark ? "Light mode" : "Dark mode";
    }
    window.dispatchEvent(new CustomEvent("recore:theme-change", {detail: {mode}}));
  }

  document.addEventListener("DOMContentLoaded", () => {
    const button = document.getElementById("themeToggle");
    if (button) {
      button.addEventListener("click", () => {
        const next = root.getAttribute("data-theme") === "dark" ? "light" : "dark";
        localStorage.setItem("recore.theme", next);
        applyTheme(next);
      });
    }

    setupNavObserver();
    setupCandidatesTable();
    setupSortableTables();
  });

  function setupNavObserver() {
    const links = Array.from(document.querySelectorAll(".nav-bar a[data-target]"));
    if (!links.length) return;
    const targets = new Map();
    links.forEach(link => {
      const id = link.dataset.target;
      const el = document.getElementById(id);
      if (el) targets.set(id, link);
    });
    if (!targets.size) return;

    links.forEach(link => {
      link.addEventListener("click", event => {
        const el = document.getElementById(link.dataset.target);
        if (!el) return;
        event.preventDefault();
        el.scrollIntoView({behavior: "smooth", block: "start"});
        history.replaceState(null, "", `#${link.dataset.target}`);
      });
    });

    const observer = new IntersectionObserver(entries => {
      entries.forEach(entry => {
        const link = targets.get(entry.target.id);
        if (!link) return;
        if (entry.isIntersecting) {
          links.forEach(other => other.classList.remove("active"));
          link.classList.add("active");
        }
      });
    }, {rootMargin: "-40% 0px -50% 0px", threshold: 0.01});

    targets.forEach((_link, id) => {
      const el = document.getElementById(id);
      if (el) observer.observe(el);
    });
  }

  function setupCandidatesTable() {
    const table = document.getElementById("candidatesTable");
    if (!table) return;
    const search = document.getElementById("candidatesSearch");
    const priority = document.getElementById("candidatesPriority");
    const count = document.getElementById("candidatesCount");
    const tbody = table.tBodies[0];
    if (!tbody) return;
    const rows = Array.from(tbody.rows);
    const total = rows.length;

    const update = () => {
      const q = (search?.value || "").trim().toLowerCase();
      const prio = (priority?.value || "").trim();
      let visible = 0;
      rows.forEach(row => {
        const haystack = row.dataset.searchText || row.textContent.toLowerCase();
        const rowPrio = row.dataset.priority || "";
        const matchText = !q || haystack.includes(q);
        const matchPrio = !prio || rowPrio === prio;
        const show = matchText && matchPrio;
        row.classList.toggle("hidden", !show);
        if (show) visible += 1;
      });
      if (count) {
        count.textContent = visible === total
          ? `${total} candidates`
          : `${visible} / ${total} candidates`;
      }
    };
    rows.forEach(row => {
      row.dataset.searchText = row.textContent.toLowerCase();
    });
    if (search) search.addEventListener("input", update);
    if (priority) priority.addEventListener("change", update);
    update();
  }

  function setupSortableTables() {
    document.querySelectorAll("table[data-sortable]").forEach(table => {
      const tbody = table.tBodies[0];
      if (!tbody) return;
      const headers = Array.from(table.tHead?.rows[0]?.cells || []);
      headers.forEach((header, idx) => {
        if (!header.classList.contains("sortable")) return;
        header.addEventListener("click", () => {
          const isAsc = header.classList.contains("sort-asc");
          headers.forEach(h => h.classList.remove("sort-asc", "sort-desc"));
          header.classList.add(isAsc ? "sort-desc" : "sort-asc");
          const dir = isAsc ? -1 : 1;
          const type = header.dataset.sortType || "string";
          const rows = Array.from(tbody.rows);
          rows.sort((a, b) => {
            const av = sortValue(a.cells[idx], type);
            const bv = sortValue(b.cells[idx], type);
            if (av < bv) return -1 * dir;
            if (av > bv) return 1 * dir;
            return 0;
          });
          rows.forEach(row => tbody.appendChild(row));
        });
      });
    });
  }

  function sortValue(cell, type) {
    if (!cell) return type === "number" ? -Infinity : "";
    const raw = (cell.dataset.sortValue ?? cell.textContent ?? "").trim();
    if (type === "number") {
      const n = parseFloat(raw.replace(/[,_]/g, ""));
      return Number.isFinite(n) ? n : -Infinity;
    }
    return raw.toLowerCase();
  }
})();
