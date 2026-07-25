(() => {
const data = window.recoreExplorerData;
if (!data || !data.clusters.length) return;
function readVar(name, fallback) {
  const value = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  return value || fallback;
}
function refreshThemeColors() {
  colors.P = readVar("--present", "#155eef");
  colors.R = readVar("--rescue", "#12b76a");
  colors.N = readVar("--nsc", "#7a5af8");
  colors.A = readVar("--absent", "#ffffff");
  strokes.P = colors.P;
  strokes.R = colors.R;
  strokes.N = colors.N;
  strokes.A = readVar("--absent-stroke", "#d0d5dd");
  panelBg = readVar("--panel", "#ffffff");
  headerBg = readVar("--bg-subtle", "#f8fafc");
  hoverColTint = readVar("--warning-soft", "#fff7ed");
  selectColTint = readVar("--accent-soft", "#eef4ff");
  selectStroke = readVar("--accent-ink", "#1d4ed8");
  gridStroke = readVar("--faint", "#98a2b3");
  labelInk = readVar("--subtle", "#475467");
}
const clusterByName = new Map(data.clusters.map(cluster => [cluster.name, cluster]));
const colors = {P:"#155eef", R:"#12b76a", N:"#7a5af8", A:"#ffffff"};
const strokes = {P:"#155eef", R:"#067647", N:"#6941c6", A:"#d0d5dd"};
let panelBg = "#ffffff";
let headerBg = "#f8fafc";
let hoverColTint = "#fff7ed";
let selectColTint = "#eef4ff";
let selectStroke = "#1d4ed8";
let gridStroke = "#98a2b3";
let labelInk = "#475467";
const labels = {P:"present in input", R:"rescued by Recore", N:"supported by NSC", A:"absent"};
const processColors = {
  clade_associated_gain:"#079455",
  ancestral_retained_no_disruption:"#155eef",
  duplication_or_paralog_expansion:"#7a5af8",
  deletion_loss_candidate:"#b42318",
  reference_or_lineage_loss_candidate:"#d92d20",
  clean_insertion_or_acquisition_candidate:"#079455",
  IS_interruption_supported:"#f79009",
  pseudogenization:"#b54708",
  gene_fragmentation_or_structural_break:"#dc6803",
  partial_reference_homology:"#0e9384",
  reference_gene_allelic_derivative:"#2970ff",
  reference_gene_disruption_or_truncation:"#c11574",
  translocation_or_context_shift_candidate:"#9e77ed",
  annotation_or_clustering_artifact:"#667085",
  unstable_or_homoplastic_accessory_signal:"#dd2590",
  unresolved_accessory_variation:"#475467",
};
const processDescriptions = {
  clade_associated_gain:"Legacy label for a clade-associated accessory gain; inspect current mechanism fields for final wording.",
  ancestral_retained_no_disruption:"Conserved ancestral presence across the assessed tree; not an obvious gain/loss event.",
  duplication_or_paralog_expansion:"Signal compatible with copy-number expansion or paralogous locus complexity.",
  deletion_loss_candidate:"Cluster is absent in a pattern compatible with deletion or lineage loss.",
  reference_or_lineage_loss_candidate:"Reference or a lineage appears to have lost a locus present elsewhere.",
  clean_insertion_or_acquisition_candidate:"Presence pattern is consistent with a gain/acquisition/insertion candidate.",
  IS_interruption_supported:"Loss/disruption is supported by nearby mobile element or insertion-sequence evidence.",
  pseudogenization:"Evidence suggests gene disruption, frameshift, internal stop, or degraded ORF.",
  gene_fragmentation_or_structural_break:"Partial homology or split alignment suggests fragmentation or a local structural break.",
  partial_reference_homology:"Cluster has partial similarity to a reference locus rather than being fully novel.",
  reference_gene_allelic_derivative:"Cluster likely represents an allelic or event-derived version of a reference gene.",
  reference_gene_disruption_or_truncation:"Reference-like sequence is present but length/frame evidence suggests disruption or truncation.",
  translocation_or_context_shift_candidate:"Sequence appears real but its genomic context differs from the reference locus.",
  annotation_or_clustering_artifact:"Pattern is more likely driven by annotation or clustering behavior than biology.",
  unstable_or_homoplastic_accessory_signal:"Patchy pattern with repeated appearances/disappearances; treat cautiously.",
  unresolved_accessory_variation:"Accessory signal remains unresolved from the available matrix/tree/evidence.",
};
const lineageColors = {
  L4:"#C94434",
  L2:"#3B5BA9",
  L3:"#8B6DAF",
  L1:"#D44D6E",
  L7:"#D4C847",
  L5:"#A07040",
  L6:"#4A9B68",
  L9:"#5CB85C",
  L10:"#7DA094",
  A1:"#C4A546",
  A2:"#6DB89B",
  A3:"#7AAED4",
  A4:"#D88EAF",
  L8:"#E08B38",
  "M. canettii":"#8C8C8C",
  unknown:"#667085",
};
const lineageDescriptions = {
  L4:"Euro-American",
  L2:"East Asian (Beijing)",
  L3:"East African-Indian (CAS)",
  L1:"Indo-Oceanic",
  L7:"Ethiopian",
  L5:"West African 1 (M. africanum)",
  L6:"West African 2 (M. africanum)",
  L9:"East African",
  L10:"Horn of Africa",
  A1:"Animal-adapted",
  A2:"Animal-adapted",
  A3:"Animal-adapted",
  A4:"Animal-adapted",
  L8:"African Great Lakes",
  "M. canettii":"Outgroup",
  unknown:"Lineage not parsed from genome name",
};
const lineageOrder = ["L4","L2","L3","L1","L7","L5","L6","L9","L10","A1","A2","A3","A4","L8","M. canettii","unknown"];
const search = document.getElementById("clusterSearch");
const omnibox = document.getElementById("clusterOmnibox");
const omniboxName = document.getElementById("omniboxName");
const omniboxNote = document.getElementById("omniboxNote");
const omniboxKey = document.getElementById("omniboxKey");
const clusterPosition = document.getElementById("clusterPosition");
const prevClusterButton = document.getElementById("prevCluster");
const nextClusterButton = document.getElementById("nextCluster");
const palette = document.getElementById("palette");
const paletteList = document.getElementById("paletteList");
const paletteCount = document.getElementById("paletteCount");
const paletteScrim = document.getElementById("paletteScrim");
const shortcutSheet = document.getElementById("shortcutSheet");
const shortcutsButton = document.getElementById("shortcutsButton");
const treeMenu = document.getElementById("treeMenu");
const treeMenuBtn = document.getElementById("treeMenuBtn");
const matrixMenu = document.getElementById("matrixMenu");
const matrixMenuBtn = document.getElementById("matrixMenuBtn");
const clusterTab = document.getElementById("clusterTab");
const genomeTab = document.getElementById("genomeTab");
const genomeDetails = document.getElementById("genomeDetails");
const overviewStrip = document.getElementById("overviewStrip");
const overviewCanvas = document.getElementById("overviewCanvas");
const overviewWindow = document.getElementById("overviewWindow");
const treeScaleBar = document.getElementById("treeScaleBar");
const closeSyntenyButton = document.getElementById("closeSynteny");
const sizeControl = document.getElementById("cellSize");
const rowHeightControl = document.getElementById("rowHeight");
const colorMode = document.getElementById("tipColorMode");
const matrixColorMode = document.getElementById("matrixColorMode");
const columnOrder = document.getElementById("columnOrder");
const rootMode = document.getElementById("rootMode");
const labelToggle = document.getElementById("tipLabels");
const referenceToggle = document.getElementById("referenceOnly");
const syntenyToggle = document.getElementById("toggleSynteny");
const syntenyZoomControl = document.getElementById("syntenyZoom");
const fullscreenButton = document.getElementById("explorerFullscreen");
const phylogramButton = document.getElementById("phylogramMode");
const cladogramButton = document.getElementById("cladogramMode");
const explorerSection = document.getElementById("explorer");
const inspectorPane = document.getElementById("inspectorPane");
const exportButton = document.getElementById("exportView");
const details = document.getElementById("clusterDetails");
const processLegend = document.getElementById("processLegend");
const lineageLegend = document.getElementById("lineageLegend");
const status = document.getElementById("matrixStatus");
const headerCanvas = document.getElementById("columnHeaderCanvas");
const headerCtx = headerCanvas.getContext("2d");
const canvas = document.getElementById("matrixCanvas");
const matrixViewport = canvas ? canvas.parentElement : null;
const matrixSpacer = document.getElementById("matrixSpacer");
const headerSpacer = document.getElementById("headerSpacer");
const treeCanvas = document.getElementById("treeCanvas");
const treeCtx = treeCanvas ? treeCanvas.getContext("2d") : null;
const treeViewport = treeCanvas ? treeCanvas.parentElement : null;
const treeSpacer = document.getElementById("treeSpacer");
const treeHeadSpacer = document.getElementById("treeHeadSpacer");
const ctx = canvas.getContext("2d");
const matrixPane = document.getElementById("matrixPane");
const treePane = document.getElementById("treePane");
const treeScroll = document.getElementById("treeScroll");
const lineageStripCanvas = document.getElementById("lineageStrip");
const lineageStripCtx = lineageStripCanvas ? lineageStripCanvas.getContext("2d") : null;
const headerScroller = document.getElementById("headerScroller");
const viewerBody = document.getElementById("viewerBody");
const viewerSplit = document.getElementById("viewerSplit");
const treeHandle = document.getElementById("treeHandle");
const syntenyHandle = document.getElementById("syntenyHandle");
const topTreeMeta = document.getElementById("topTreeMeta");
const topMatrixMeta = document.getElementById("topMatrixMeta");
const topSyntenyMeta = document.getElementById("topSyntenyMeta");
const tooltip = document.getElementById("matrixTooltip");
const syntenyPanel = document.getElementById("syntenyPanel");
const syntenySummary = document.getElementById("syntenySummary");
const syntenySvg = document.getElementById("syntenySvg");
const LAYOUT_KEY = "recore.explorer.layout.v2";
const DEFAULT_TREE_W = 320;
const DEFAULT_SYNTENY_H = 300;
const MIN_TREE_W = 180;
const MIN_MATRIX_W = 220;
const MIN_SYNTENY_H = 170;
const MAX_SYNTENY_FRACTION = 0.72;
const DEFAULT_INSPECTOR_W = 320;
const MIN_INSPECTOR_W = 220;
const LINEAGE_STRIP_W = 6;
// Single source for the zoom bounds. The sliders, currentCellSize/currentRowHeight
// and the wheel handler must agree, or the wheel snaps the value into its own range.
const MIN_COL_ZOOM = 0.25;
const MAX_COL_ZOOM = 90;
const MIN_ROW_ZOOM = 1;
const MAX_ROW_ZOOM = 40;
// Target pixel span used to pick the opening zoom from the run size.
const FIT_TARGET_H = 560;
const FIT_TARGET_W = 720;
// Column names are only worth drawing once a cell is wide enough to carry a
// rotated label without the text of neighbouring columns touching.
const LABEL_MIN_CELL_W = 9;
const LABEL_BAND_H = 82;
// Above this a name fits across the column, and turning your head to read a
// dozen columns is a cost the layout no longer has to impose.
const LABEL_HORIZONTAL_CELL_W = 54;
const LABEL_BAND_H_FLAT = 26;
const LABEL_FONT = "10px Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif";
const defaultStatus = "Hover a cell for genome, cluster and evidence state. Wheel over the column header zooms columns, over the tree zooms rows. Drag the seams between panes to resize.";
let selected = 0;
let hoverRow = -1;
let hoverCol = -1;
let leafOrder = [];
/* Matrix rows follow tree-tip order while cell sets are keyed by matrix genome
   index, so drawing needs the inverse mapping. -1 marks a genome with no tip. */
let rowOfGenome = new Int32Array(0);
let visibleClusters = [];
let genomeStats = [];
let maxGenomeStats = {R:1, N:1, A:1, correction:1};
let maxDepth = 0;
let maxDistance = 0;
let treeHeight = 0;
let treeWidth = DEFAULT_TREE_W;
const headerH = 24;
// Chrome refuses a 2D canvas wider or taller than 16384 device pixels, and does
// it silently: allocation succeeds and every draw is discarded, so the matrix
// comes out blank rather than erroring. The old 32000 sat above that line.
const maxCanvasDim = 16384;
const maxCanvasPixels = 90000000;
let cellH = 12;
let cellW = 8;
let virtualW = 1;
let virtualH = 1;
// Cells may be fractional, but never thinner than this or they vanish.
const MIN_CELL_PX = 0.35;
let treeDirty = true;
let renderQueued = false;
let detailsKey = "";
let syntenyKey = "";
let treeMode = "phylogram";
let referenceOnly = false;
let syntenyHidden = true;
let selectedSyntenyBlock = 0;
let selectedSyntenyReference = false;
let selectedSyntenyCarriers = new Set();
let mechanismFilters = new Set();
/* Rows matching the genome search, marked on the tree and the lineage strip. */
let highlightedRows = new Set();
let overviewKey = "";
let treeRoot = data.tree;

/* ---- Cell decoding -------------------------------------------------------
   Cells arrive as one base64 blob per cluster holding the genome indices that
   are not absent, split by state (see cellcodec.rs).

   Rows stay in set form rather than being expanded into a dense array. That is
   what keeps this viable at scale: a decoded cluster costs about what its
   encoding costs, so every cluster can be held at once (a few MB for 30,000
   clusters) instead of 10 KB each, and drawing a column costs the number of
   cells that actually carry information rather than the number of genomes.

   A set is one of:
     {kind:"all"}                 every genome, drawn as one filled rect
     {kind:"list", idx}           sorted genome indices
     {kind:"complement", idx}     every genome except these, drawn as a filled
                                  rect with the exclusions punched back out */
const STATE_CODES = ["P", "R", "N", "A"];
const EMPTY_SET = {kind: "list", idx: new Int32Array(0)};
const decodedRows = new Map();

function decodeBase64(text) {
  const binary = atob(text);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i += 1) bytes[i] = binary.charCodeAt(i);
  return bytes;
}
function decodeClusterSets(encoded, nGenomes) {
  const sets = [EMPTY_SET, EMPTY_SET, EMPTY_SET];
  if (!encoded) return sets;
  const bytes = decodeBase64(encoded);
  let pos = 0;
  const readVarint = () => {
    let value = 0;
    let shift = 0;
    for (;;) {
      const byte = bytes[pos];
      pos += 1;
      value += (byte & 0x7f) * Math.pow(2, shift);
      if ((byte & 0x80) === 0) return value;
      shift += 7;
    }
  };
  const readIndices = () => {
    const count = readVarint();
    const idx = new Int32Array(count);
    let previous = 0;
    for (let i = 0; i < count; i += 1) {
      previous += readVarint();
      idx[i] = previous;
    }
    return idx;
  };
  for (let state = 0; state < 3; state += 1) {
    const mode = bytes[pos];
    pos += 1;
    if (mode === 0) {
      sets[state] = {kind: "list", idx: readIndices()};
    } else if (mode === 1) {
      const idx = readIndices();
      sets[state] = idx.length === 0 ? {kind: "all"} : {kind: "complement", idx};
    } else if (mode === 2) {
      // Bitset mode only wins near 50% prevalence, where a plain list is the
      // cheapest thing to draw from anyway.
      const len = Math.ceil(nGenomes / 8);
      const idx = [];
      for (let g = 0; g < nGenomes; g += 1) {
        if (bytes[pos + (g >> 3)] & (1 << (g & 7))) idx.push(g);
      }
      pos += len;
      sets[state] = {kind: "list", idx: Int32Array.from(idx)};
    } else {
      return sets;
    }
  }
  return sets;
}
function buildRowOfGenome() {
  rowOfGenome = new Int32Array(data.genomes.length).fill(-1);
  for (let row = 0; row < leafOrder.length; row += 1) {
    const genomeIdx = leafOrder[row];
    if (genomeIdx >= 0 && genomeIdx < rowOfGenome.length) rowOfGenome[genomeIdx] = row;
  }
}
function clusterSets(clusterIdx) {
  let sets = decodedRows.get(clusterIdx);
  if (!sets) {
    sets = decodeClusterSets(data.cells[clusterIdx], data.genomes.length);
    decodedRows.set(clusterIdx, sets);
  }
  return sets;
}
function setSize(set, nGenomes) {
  if (set.kind === "all") return nGenomes;
  if (set.kind === "complement") return nGenomes - set.idx.length;
  return set.idx.length;
}
/** Index of the first entry >= value, or idx.length. */
function lowerBound(idx, value) {
  let lo = 0;
  let hi = idx.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (idx[mid] < value) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}
function setContains(set, value) {
  if (set.kind === "all") return true;
  const at = lowerBound(set.idx, value);
  const found = at < set.idx.length && set.idx[at] === value;
  return set.kind === "complement" ? !found : found;
}
function stateAt(clusterIdx, genomeIdx) {
  const sets = clusterSets(clusterIdx);
  // Later states win, matching the writer's ordering.
  for (let state = 2; state >= 0; state -= 1) {
    if (setContains(sets[state], genomeIdx)) return STATE_CODES[state];
  }
  return "A";
}
function countsForCluster(clusterIdx) {
  const n = data.genomes.length;
  const sets = clusterSets(clusterIdx);
  const counts = {
    P: setSize(sets[0], n),
    R: setSize(sets[1], n),
    N: setSize(sets[2], n),
    A: 0,
  };
  counts.A = Math.max(0, n - counts.P - counts.R - counts.N);
  return counts;
}
function passesActiveFilters(idx, includeMechanism = true) {
  const cluster = data.clusters[idx];
  if (referenceOnly && !cluster.referenceAbsent) return false;
  if (includeMechanism && mechanismFilters.size && !mechanismFilters.has(processKind(cluster))) return false;
  return true;
}
function mechanismFilterKey() {
  return [...mechanismFilters].sort().join("|");
}
function mechanismFilterSummary() {
  const filters = [...mechanismFilters].sort();
  if (!filters.length) return "";
  const labels = filters.map(processLabel);
  if (labels.length <= 2) return labels.join(" + ");
  return `${labels.length} mechanisms`;
}
/* Column order is a real analysis control, not decoration. Panaroo emits
   clusters in its own order, which over 30,000 columns is arbitrary: sorting by
   prevalence or mechanism is what turns the matrix from a colour wash into
   something where clade-associated blocks group together. */
function sortVisibleClusters(indices) {
  const mode = columnOrder ? columnOrder.value : "matrix";
  if (mode === "matrix") return indices;
  const key = idx => {
    const cluster = data.clusters[idx] || {};
    if (mode === "prevalence") return -(cluster.prevalence || 0);
    if (mode === "name") return String(cluster.name || "");
    if (mode === "mechanism") return processLabel(processKind(cluster));
    if (mode === "corrected") {
      const counts = countsForCluster(idx);
      return -(counts.R + counts.N);
    }
    return idx;
  };
  return indices
    .map(idx => ({idx, key: key(idx)}))
    .sort((a, b) => {
      if (a.key < b.key) return -1;
      if (a.key > b.key) return 1;
      // Ties keep matrix order, so the view is stable between renders.
      return a.idx - b.idx;
    })
    .map(entry => entry.idx);
}
function refreshVisibleClusters() {
  const previousSelected = selected;
  visibleClusters = sortVisibleClusters(data.clusters
    .map((_, idx) => idx)
    .filter(idx => passesActiveFilters(idx, true)));
  if (!visibleClusters.length) {
    if (referenceOnly) {
      referenceOnly = false;
      visibleClusters = sortVisibleClusters(data.clusters
        .map((_, idx) => idx)
        .filter(idx => passesActiveFilters(idx, true)));
    }
    if (!visibleClusters.length && mechanismFilters.size) {
      mechanismFilters.clear();
      visibleClusters = sortVisibleClusters(data.clusters.map((_, idx) => idx));
    }
  }
  if (!visibleClusters.includes(selected)) {
    selected = visibleClusters.length ? visibleClusters[0] : 0;
  }
  if (selected !== previousSelected) {
    selectedSyntenyBlock = 0;
    selectedSyntenyReference = false;
    selectedSyntenyCarriers = new Set();
    syntenyKey = "";
  }
  if (referenceToggle) {
    referenceToggle.checked = referenceOnly;
    referenceToggle.setAttribute(
      "data-tip",
      data.referenceGenome
        ? `Show only clusters absent from ${data.referenceGenome.name}`
        : "No reference genome was supplied for this run",
    );
    referenceToggle.disabled = !data.referenceGenome;
  }
}
function selectedVisibleColumn() {
  return visibleClusters.indexOf(selected);
}
function originalClusterAtVisibleColumn(col) {
  return visibleClusters[col] ?? -1;
}
function computeGenomeStats() {
  // Totals come precomputed from the writer: deriving them here would mean
  // decoding every cluster before the first paint.
  genomeStats = data.genomes.map((_, gi) => {
    const totals = (data.genomeStats && data.genomeStats[gi]) || [0, 0, 0, 0];
    return {P: totals[0], R: totals[1], N: totals[2], A: totals[3]};
  });
  maxGenomeStats = genomeStats.reduce((acc, s) => ({
    R: Math.max(acc.R, s.R),
    N: Math.max(acc.N, s.N),
    A: Math.max(acc.A, s.A),
    correction: Math.max(acc.correction, s.R + s.N),
  }), {R:1, N:1, A:1, correction:1});
}
function gradientColor(value, maxValue, low, high) {
  if (!maxValue) return low;
  const t = Math.max(0, Math.min(1, value / maxValue));
  const parse = hex => [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16));
  const a = parse(low);
  const b = parse(high);
  const rgb = a.map((v, i) => Math.round(v + (b[i] - v) * t));
  return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
}
function processColor(mechanism) {
  if (processColors[mechanism]) return processColors[mechanism];
  let hash = 0;
  const text = mechanism || "unresolved_accessory_variation";
  for (let i = 0; i < text.length; i += 1) hash = (hash * 31 + text.charCodeAt(i)) >>> 0;
  const hue = hash % 360;
  return `hsl(${hue} 58% 43%)`;
}
function processKind(cluster) {
  return (cluster && cluster.mechanism) || "unresolved_accessory_variation";
}
function processLabel(mechanism) {
  return (mechanism || "unresolved_accessory_variation").replace(/_/g, " ");
}
function processDescription(mechanism) {
  return processDescriptions[mechanism] || "Evolutionary process not explicitly classified; inspect TSV evidence before interpretation.";
}
function lineageOfGenomeName(name) {
  const text = String(name || "");
  if (text.startsWith("M. canettii")) return "M. canettii";
  const match = text.match(/^(L10|L[1-9]|A[1-4])(?=[_-])/);
  return match ? match[1] : "unknown";
}
function lineageColor(lineage) {
  return lineageColors[lineage] || lineageColors.unknown;
}
function lineageLabel(lineage) {
  const description = lineageDescriptions[lineage] || "";
  return description ? `${lineage} - ${description}` : lineage;
}
function renderLineageLegend() {
  if (!lineageLegend) return;
  if (colorMode.value !== "lineage") {
    lineageLegend.innerHTML = "";
    lineageLegend.style.display = "none";
    return;
  }
  const counts = new Map();
  leafOrder.forEach(genomeIdx => {
    const lineage = lineageOfGenomeName(data.genomes[genomeIdx]);
    counts.set(lineage, (counts.get(lineage) || 0) + 1);
  });
  const items = lineageOrder
    .filter(lineage => counts.has(lineage))
    .map(lineage => `<span class="lineage-chip"><i style="background:${lineageColor(lineage)}"></i><strong>${escapeHtml(lineage)}</strong><small>${counts.get(lineage)} tips</small></span>`)
    .join("");
  lineageLegend.style.display = "flex";
  lineageLegend.innerHTML = `<span class="lineage-title">MTBC lineage colors</span>${items}`;
}
function visibleProcessCounts() {
  const counts = new Map();
  data.clusters
    .map((_, idx) => idx)
    .filter(idx => passesActiveFilters(idx, false))
    .forEach(idx => {
    const mechanism = processKind(data.clusters[idx]);
    const row = counts.get(mechanism) || {clusters:0, present:0};
    row.clusters += 1;
    row.present += Number(data.clusters[idx]?.present || 0);
    counts.set(mechanism, row);
  });
  return [...counts.entries()].sort((a, b) =>
    b[1].clusters - a[1].clusters || processLabel(a[0]).localeCompare(processLabel(b[0]))
  );
}
function renderProcessLegend() {
  if (!processLegend) return;
  const cluster = data.clusters[selected] || {};
  const selectedMechanism = processKind(cluster);
  const filterLabel = mechanismFilters.size
    ? `<button type="button" class="process-clear" data-clear-mechanism="1">Showing ${escapeHtml(mechanismFilterSummary())} (${visibleClusters.length} clusters). Clear filters</button>`
    : `<span>Click mechanism cards to add/remove matrix-column filters.</span>`;
  const items = visibleProcessCounts().map(([mechanism, row]) => {
    const active = mechanism === selectedMechanism ? " active" : "";
    const filtered = mechanismFilters.has(mechanism) ? " filtered" : "";
    return `<button type="button" class="process-item${active}${filtered}" data-mechanism="${escapeHtml(mechanism)}"><span class="process-swatch" style="background:${processColor(mechanism)}"></span><div><strong>${escapeHtml(processLabel(mechanism))}</strong><small>${plural(row.clusters, "cluster")}; ${plural(row.present, "present cell")}. ${escapeHtml(processDescription(mechanism))}</small></div></button>`;
  }).join("");
  processLegend.innerHTML = `<div class="process-selected"><span class="process-chip" style="--process-color:${processColor(selectedMechanism)}"><i></i>${escapeHtml(processLabel(selectedMechanism))}</span><span>${escapeHtml(processDescription(selectedMechanism))}</span>${filterLabel}</div><div class="process-list">${items}</div>`;
}
function tipStyle(genomeIdx) {
  const mode = colorMode.value;
  if (mode === "cluster") {
    const state = stateAt(selected, genomeIdx);
    return {fill: colors[state], stroke: strokes[state], title: labels[state]};
  }
  if (mode === "process") {
    const cluster = data.clusters[selected] || {};
    const state = stateAt(selected, genomeIdx);
    const mechanism = processKind(cluster);
    const color = processColor(mechanism);
    return state === "A"
      ? {fill: colors.A, stroke: strokes.A, title: `${processLabel(mechanism)} absent`}
      : {fill: color, stroke: color, title: `${processLabel(mechanism)}; ${labels[state]}`};
  }
  if (mode === "lineage") {
    const genome = data.genomes[genomeIdx] || "";
    const lineage = lineageOfGenomeName(genome);
    const color = lineageColor(lineage);
    return {fill: color, stroke: color, title: lineageLabel(lineage)};
  }
  const stats = genomeStats[genomeIdx] || {R:0, N:0, A:0};
  if (mode === "rescue") {
    return {fill: gradientColor(stats.R, maxGenomeStats.R, "#e0f2fe", "#0369a1"), stroke: "#0369a1", title: `rescued cells ${stats.R}`};
  }
  if (mode === "absence") {
    return {fill: gradientColor(stats.A, maxGenomeStats.A, "#f8fafc", "#b42318"), stroke: "#b42318", title: `absent cells ${stats.A}`};
  }
  const correction = stats.R + stats.N;
  return {fill: gradientColor(correction, maxGenomeStats.correction, "#ecfdf3", "#067647"), stroke: "#067647", title: `corrected/NSC cells ${correction}`};
}
function subtreeLineage(node) {
  if (node._subtreeLineage !== undefined) return node._subtreeLineage;
  if (!node.children || !node.children.length) {
    node._subtreeLineage = lineageOfGenomeName(tipName(node));
    return node._subtreeLineage;
  }
  const childLineages = (node.children || []).map(subtreeLineage);
  const first = childLineages[0] || "unknown";
  node._subtreeLineage = childLineages.every(lineage => lineage === first) ? first : "mixed";
  return node._subtreeLineage;
}
function treeEdgeStyle(child) {
  if (colorMode.value !== "lineage") return "";
  const lineage = subtreeLineage(child);
  if (!lineage || lineage === "mixed") return "";
  return ` style="stroke:${lineageColor(lineage)};stroke-width:1.35;stroke-opacity:.82"`;
}
function selectedSyntenyBlockData() {
  const cluster = data.clusters[selected];
  const blocks = cluster && data.synteny ? (data.synteny[cluster.name] || []) : [];
  if (!blocks.length) return null;
  selectedSyntenyBlock = Math.max(0, Math.min(selectedSyntenyBlock, blocks.length - 1));
  return blocks[selectedSyntenyBlock] || null;
}
function updateSelectedSyntenyCarriers(blocks) {
  selectedSyntenyBlock = Math.max(0, Math.min(selectedSyntenyBlock, Math.max(0, blocks.length - 1)));
  if (selectedSyntenyReference) {
    const refIdx = data.referenceGenome && Number.isInteger(data.referenceGenome.index)
      ? data.referenceGenome.index
      : null;
    selectedSyntenyCarriers = refIdx === null ? new Set() : new Set([refIdx]);
    return;
  }
  const block = blocks[selectedSyntenyBlock];
  selectedSyntenyCarriers = new Set(block ? syntenyCarrierIndices(block) : []);
}
function isSelectedSyntenyCarrier(genomeIdx) {
  return selectedSyntenyCarriers.has(genomeIdx);
}
function syntenyCarrierIndices(block) {
  if (Array.isArray(block.genomes)) return block.genomes;
  const bits = String(block.genomeBits || "");
  const out = [];
  for (let i = 0; i < bits.length; i += 1) {
    const value = parseInt(bits[i], 16);
    if (!Number.isFinite(value)) continue;
    for (let bit = 0; bit < 4; bit += 1) {
      if (value & (1 << bit)) out.push(i * 4 + bit);
    }
  }
  return out;
}
function cloneNode(node) {
  return {
    label: node.label || "",
    genome: node.genome,
    length: Number(node.length || 0),
    children: (node.children || []).map(cloneNode),
  };
}
function tipName(node) {
  return node.genome === null || node.genome === undefined
    ? (node.label || "")
    : (data.genomes[node.genome] || node.label || "");
}
function isL8Tip(node) {
  return (!node.children || !node.children.length) && /^L8[_-]/.test(tipName(node));
}
function assignTreeIds(node, parent, nodes, edges) {
  const id = nodes.length;
  node._id = id;
  node._parent = parent;
  nodes.push(node);
  (node.children || []).forEach(child => {
    const childId = assignTreeIds(child, node, nodes, edges);
    const length = Number(child.length || 0);
    edges.push([id, childId, length]);
  });
  return id;
}
function collectL8Leaves(node, out) {
  if (isL8Tip(node)) out.push(node);
  (node.children || []).forEach(child => collectL8Leaves(child, out));
}
function commonAncestor(nodes) {
  if (!nodes.length) return null;
  const paths = nodes.map(node => {
    const path = [];
    let cursor = node;
    while (cursor) {
      path.push(cursor);
      cursor = cursor._parent;
    }
    return path.reverse();
  });
  let mrca = paths[0][0];
  for (let i = 0; ; i += 1) {
    const candidate = paths[0][i];
    if (!candidate || paths.some(path => path[i] !== candidate)) break;
    mrca = candidate;
  }
  return mrca;
}
function rootByL8(node) {
  const root = cloneNode(node);
  const nodes = [];
  const edges = [];
  assignTreeIds(root, null, nodes, edges);
  const leaves = [];
  collectL8Leaves(root, leaves);
  if (!leaves.length) return root;
  const rootNode = leaves.length === 1 && leaves[0]._parent ? leaves[0]._parent : commonAncestor(leaves);
  const adj = nodes.map(() => []);
  edges.forEach(([a, b, length]) => {
    adj[a].push({id:b, length});
    adj[b].push({id:a, length});
  });
  function build(id, parentId, length) {
    const source = nodes[id];
    const children = adj[id]
      .filter(edge => edge.id !== parentId)
      .map(edge => build(edge.id, id, edge.length));
    children.sort((a, b) => Number(hasL8Tip(b)) - Number(hasL8Tip(a)) || countTips(b) - countTips(a));
    return {label: source.label || "", genome: source.genome, length, children};
  }
  return build(rootNode ? rootNode._id : root._id, null, 0);
}
function hasL8Tip(node) {
  return isL8Tip(node) || (node.children || []).some(hasL8Tip);
}
function countTips(node) {
  if (!node.children || !node.children.length) return 1;
  return node.children.reduce((sum, child) => sum + countTips(child), 0);
}
function findLeafByGenome(node, genomeIdx) {
  if ((!node.children || !node.children.length) && node.genome === genomeIdx) return node;
  for (const child of node.children || []) {
    const found = findLeafByGenome(child, genomeIdx);
    if (found) return found;
  }
  return null;
}
function layout(node, depth, distance = 0) {
  node.depth = depth;
  node.distance = distance;
  maxDepth = Math.max(maxDepth, depth);
  maxDistance = Math.max(maxDistance, distance);
  if (!node.children || node.children.length === 0) {
    node.row = leafOrder.length;
    leafOrder.push(node.genome);
    return node.row;
  }
  let sum = 0;
  let count = 0;
  node.children.forEach(child => { sum += layout(child, depth + 1, distance + Number(child.length || 0)); count += 1; });
  node.row = count ? sum / count : leafOrder.length;
  return node.row;
}
/* Keeps the selection inside the visible set and repaints the find field. This
   used to build one <option> per cluster on every keystroke, which at 30,000
   clusters meant 30,000 DOM nodes for a list nobody could scroll anyway. */
function populateOptions() {
  const previousSelected = selected;
  if (visibleClusters.length && !visibleClusters.includes(selected)) {
    selected = visibleClusters[0];
  }
  if (selected !== previousSelected) {
    selectedSyntenyBlock = 0;
    selectedSyntenyReference = false;
    selectedSyntenyCarriers = new Set();
    syntenyKey = "";
  }
  const cluster = data.clusters[selected];
  if (omniboxName) omniboxName.textContent = cluster ? cluster.name : "No cluster";
  if (omniboxNote) {
    omniboxNote.textContent = cluster
      ? (cluster.annotation || "No annotation provided")
      : "nothing matches the current filters";
  }
  if (clusterPosition) {
    const at = visibleClusters.indexOf(selected);
    clusterPosition.textContent = at < 0
      ? `0 / ${visibleClusters.length}`
      : `${at + 1} / ${visibleClusters.length}`;
  }
  if (omnibox) {
    omnibox.setAttribute(
      "data-tip",
      cluster
        ? `${cluster.name}: ${cluster.annotation || "no annotation"}. Search any cluster or genome.`
        : "Search any cluster or genome.",
    );
  }
  if (prevClusterButton) prevClusterButton.disabled = visibleClusters.indexOf(selected) <= 0;
  if (nextClusterButton) {
    const at = visibleClusters.indexOf(selected);
    nextClusterButton.disabled = at < 0 || at >= visibleClusters.length - 1;
  }
}

/** Move the selection by `delta` columns within the visible set. */
function stepCluster(delta) {
  if (!visibleClusters.length) return;
  const at = visibleClusters.indexOf(selected);
  const next = Math.max(0, Math.min(visibleClusters.length - 1, (at < 0 ? 0 : at) + delta));
  if (visibleClusters[next] === selected) return;
  selected = visibleClusters[next];
  selectedSyntenyBlock = 0;
  selectedSyntenyReference = false;
  syntenyKey = "";
  populateOptions();
  treeDirty = true;
  render(true);
}
function currentCellSize() {
  const value = Number(sizeControl.value || 8);
  return Math.max(MIN_COL_ZOOM, Math.min(MAX_COL_ZOOM, value));
}
function currentRowHeight() {
  const value = Number(rowHeightControl.value || 8);
  return Math.max(MIN_ROW_ZOOM, Math.min(MAX_ROW_ZOOM, value));
}
function clusterLabel(cluster) {
  return cluster.annotation ? `${cluster.name} - ${cluster.annotation}` : cluster.name;
}
function referenceInterpretation(cluster) {
  const referenceCell = (cluster.referenceCell || "").trim();
  const targetGene = (cluster.targetReferenceGene || "").trim();
  const targetCluster = (cluster.targetReferenceCluster || "").trim();
  const eventClass = (cluster.referenceEventClass || "").trim();
  const parts = [];
  if (referenceCell) {
    parts.push(`<span>reference gene(s): ${escapeHtml(referenceCell)}</span>`);
  }
  if (targetGene || targetCluster || eventClass) {
    const sourceBits = [];
    if (targetGene) sourceBits.push(`source gene ${targetGene}`);
    if (targetCluster) sourceBits.push(`source cluster ${targetCluster}`);
    if (eventClass) sourceBits.push(eventClass);
    parts.push(`<span>reference-derived/event source: ${escapeHtml(sourceBits.join("; "))}</span>`);
  }
  if (!parts.length) {
    const mechanism = processKind(cluster);
    const text = mechanism === "clean_insertion_or_acquisition_candidate"
      ? "no reference equivalent detected; candidate de novo/reference-absent locus"
      : "no reference equivalent detected by Recore";
    parts.push(`<span>${escapeHtml(text)}</span>`);
  }
  const state = cluster.referenceState ? `reference state ${cluster.referenceState}` : "reference state not assessed";
  return `<div class="reference-meta"><strong>Reference equivalent</strong>${parts.join("")}<span>${escapeHtml(state)}</span></div>`;
}
/** Horizontal room the branches actually get, once labels are accounted for. */
function treeUsableWidth() {
  const labelPad = labelToggle.checked && leafOrder.length <= 1500 ? 150 : 24;
  return Math.max(80, treeWidth - 22 - labelPad);
}
function xForNode(node) {
  const usable = treeUsableWidth();
  if (treeMode === "phylogram" && maxDistance > 0) {
    return 14 + node.distance / maxDistance * usable;
  }
  const tipX = 14 + usable;
  if (!node.children || !node.children.length) {
    return tipX;
  }
  return 14 + node.depth * (usable / Math.max(1, maxDepth));
}
function drawNode(node, parts) {
  const x = xForNode(node);
  const y = node.row * cellH + cellH / 2;
  if (node.children && node.children.length) {
    node.children.forEach(child => {
      const cx = xForNode(child);
      const cy = child.row * cellH + cellH / 2;
      const branchTitle = treeMode === "phylogram" ? `<title>branch length ${formatBranchLength(child.length)}</title>` : "";
      const edgeStyle = treeEdgeStyle(child);
      parts.push(`<line x1="${x}" y1="${cy}" x2="${cx}" y2="${cy}" class="tree-edge"${edgeStyle}>${branchTitle}</line>`);
      parts.push(`<line x1="${x}" y1="${y}" x2="${x}" y2="${cy}" class="tree-edge"${edgeStyle}/>`);
      drawNode(child, parts);
    });
  } else {
    const style = tipStyle(node.genome);
    const label = data.genomes[node.genome] || node.label || "";
    const radius = Math.max(0.18, Math.min(5, cellH / 2 - 0.08));
    const strokeWidth = Math.max(0.2, Math.min(1.2, cellH * 0.45));
    parts.push(`<circle cx="${x}" cy="${y}" r="${radius}" fill="${style.fill}" stroke="${style.stroke}" stroke-width="${strokeWidth}"><title>${escapeHtml(label)}: ${escapeHtml(style.title)}</title></circle>`);
    if (highlightedRows.has(node.row)) {
    treeCtx.strokeStyle = "#0ea5e9";
    treeCtx.lineWidth = 1.6;
    treeCtx.beginPath();
    treeCtx.arc(x, y, radius + Math.max(2, Math.min(4, cellH * 0.4)), 0, Math.PI * 2);
    treeCtx.stroke();
  }
  if (isSelectedSyntenyCarrier(node.genome)) {
      const ringRadius = radius + Math.max(1.2, Math.min(3, cellH * 0.30));
      const ringTitle = selectedSyntenyReference ? "selected H37Rv reference context" : `selected synteny #${selectedSyntenyBlock + 1}`;
      parts.push(`<circle cx="${x}" cy="${y}" r="${ringRadius}" class="synteny-tip-ring"><title>${escapeHtml(label)}: ${escapeHtml(ringTitle)}</title></circle>`);
    }
    if (labelToggle.checked && leafOrder.length <= 1500 && cellH >= 4) {
      parts.push(`<text x="${x + 9}" y="${y + 3}" class="tip-label">${escapeHtml(label)}</text>`);
    }
  }
}
function formatBranchLength(value) {
  const n = Number(value || 0);
  if (!Number.isFinite(n)) return "0";
  if (n === 0) return "0";
  if (Math.abs(n) < 0.001) return n.toExponential(2);
  return n.toPrecision(4);
}
function plural(count, noun) {
  return `${count} ${noun}${count === 1 ? "" : "s"}`;
}
function escapeHtml(text) {
  return String(text).replace(/[&<>"']/g, ch => ({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;","'":"&#39;"}[ch]));
}
function clipText(text, maxLen) {
  const value = String(text || "");
  return value.length > maxLen ? `${value.slice(0, Math.max(1, maxLen - 3))}...` : value;
}
function compactAnnotation(text) {
  return String(text || "").split(/[;|]/)[0].trim();
}
function referenceLocusForCluster(name) {
  const cluster = clusterByName.get(name);
  if (!cluster) return "";
  const referenceCell = String(cluster.referenceCell || "")
    .split(";")
    .map(item => item.trim())
    .find(Boolean);
  if (referenceCell) return displayReferenceLocus(referenceCell);
  const targetGene = String(cluster.targetReferenceGene || "").trim();
  if (targetGene) return displayReferenceLocus(targetGene);
  return "";
}
function displayReferenceLocus(value) {
  const text = String(value || "").replace(/(_pseudo|_dup_\d+)$/g, "");
  const parts = text.split(/[;:,|\s]+/);
  const rv = parts.find(part => /^Rv\d{4}c?$/.test(part));
  if (rv) return rv;
  const embedded = text.match(/(?:^|_)(Rv\d{4}c?)(?:_|$)/);
  return embedded ? embedded[1] : text;
}
function syntenyGeneLabel(name) {
  const cluster = clusterByName.get(name);
  const annotation = compactAnnotation(cluster && cluster.annotation);
  const locus = referenceLocusForCluster(name);
  if (locus && annotation && annotation !== name) return `${locus} (${name}) ${annotation}`;
  if (locus) return `${locus} (${name})`;
  if (annotation && annotation !== name) return `${name} ${annotation}`;
  return name || "NA";
}
function firstExampleLocus(block) {
  const first = String(block.examples || "").split(";").map(item => item.trim()).find(Boolean) || "";
  if (!first) return "example locus unavailable";
  const parts = first.split(":");
  return parts.length > 1 ? parts.slice(1).join(":") : first;
}
function fallbackSyntenyGenes(block, clusterName) {
  const upstream = [...(block.upstream || [])].reverse().map(name => ({name, cluster:name, geneId:"", strand:"+", target:false}));
  const downstream = (block.downstream || []).map(name => ({name, cluster:name, geneId:"", strand:"+", target:false}));
  return [...upstream, {name:clusterName, cluster:clusterName, geneId:"", strand:"+", target:true}, ...downstream];
}
function normalizeSyntenyGene(gene, clusterName, reference = false) {
  const source = Array.isArray(gene)
    ? {cluster:gene[0], geneId:gene[1], strand:gene[2], contig:gene[3], start:gene[4], end:gene[5], length:gene[6]}
    : gene;
  const cluster = source.cluster || source.name || "";
  const gap = cluster === "__GAP__" || source.gap;
  return {
    name: cluster,
    cluster,
    geneId: source.geneId || source.gene_id || "",
    strand: source.strand || "+",
    contig: source.contig || "",
    start: source.start || "",
    end: source.end || "",
    length: source.length || "",
    target: cluster === clusterName,
    reference,
    gap,
  };
}
function blockGenes(block, clusterName) {
  if (Array.isArray(block.genes) && block.genes.length) {
    return block.genes.map(gene => normalizeSyntenyGene(gene, clusterName, false));
  }
  return fallbackSyntenyGenes(block, clusterName);
}
function referenceGenesForBlock(block, clusterName) {
  if (Array.isArray(block.referenceGenes) && block.referenceGenes.length) {
    return block.referenceGenes.map(gene => normalizeSyntenyGene(gene, clusterName, true));
  }
  return [{name:"__GAP__", cluster:"__GAP__", geneId:"H37Rv context not resolved", strand:".", contig:"", start:"", end:"", length:"", target:false, reference:true, gap:true}];
}
function syntenyDisplayLabel(gene) {
  if (gene.gap) return gene.geneId || "H37Rv gap";
  const direct = gene.geneId ? displayReferenceLocus(gene.geneId) : "";
  if (direct && /^Rv\d{4}c?$/.test(direct)) return direct;
  if (gene.reference && direct) return direct;
  return referenceLocusForCluster(gene.cluster) || direct || gene.cluster || "NA";
}
function syntenyIdentitySignature(genes) {
  return (genes || [])
    .map(gene => gene.gap ? "__GAP__" : syntenyDisplayLabel(gene))
    .filter(label => label && label !== "NA")
    .join("|");
}
function syntenySecondaryLabel(gene) {
  if (gene.gap) return gene.length ? `${gene.length} bp gap` : "context gap";
  const cluster = gene.cluster || "";
  const annotation = compactAnnotation((clusterByName.get(cluster) || {}).annotation);
  return annotation || cluster;
}
function syntenyCoordLabel(gene) {
  const coords = gene.contig && gene.start && gene.end
    ? `${gene.contig}:${gene.start}-${gene.end}`
    : "";
  const length = gene.length ? `${gene.length} bp` : "";
  const strand = gene.strand && gene.strand !== "." ? gene.strand : "";
  return [coords, strand, length].filter(Boolean).join(" ");
}
function syntenyLinkKey(gene) {
  if (!gene || gene.gap || !gene.cluster || gene.cluster === "NA") return "";
  return gene.cluster;
}
function clusterIndexForName(name) {
  const idx = data.clusters.findIndex(cluster => cluster.name === name);
  return idx >= 0 ? idx : -1;
}
function selectMatrixClusterByName(name) {
  const idx = clusterIndexForName(name);
  if (idx < 0) return false;
  if (!passesActiveFilters(idx, true)) {
    referenceOnly = false;
    mechanismFilters.clear();
    refreshVisibleClusters();
  }
  selected = idx;
  selectedSyntenyBlock = 0;
  selectedSyntenyReference = false;
  selectedSyntenyCarriers = new Set();
  if (referenceToggle) referenceToggle.checked = referenceOnly;
  populateOptions();
  resetMatrixCache();
  render(true);
  return true;
}
function syntenyZoom() {
  const value = Number(syntenyZoomControl?.value || 1);
  return Math.max(0.45, Math.min(2.4, value));
}
function syntenyLength(gene) {
  const value = Number(gene && gene.length);
  return Number.isFinite(value) && value > 0 ? value : 600;
}
/* Arrow width is a real base-pair scale, not a rank. The previous formula was
   34px + 132px * (len / p90), so a 720bp gene and a 1200bp gene came out 113px
   and 166px: a 1.67x difference drawn as 1.47x, and everything above the 90th
   percentile drawn identically. A gene you cannot measure by eye is a box with
   a label on it, not a gene. */
const SYNTENY_TARGET_PX = 120;
const SYNTENY_MIN_GENE_PX = 26;
function syntenyLengthScale(rows) {
  const lengths = rows
    .flatMap(row => row.genes || [])
    .map(syntenyLength)
    .filter(value => value > 0)
    .sort((a, b) => a - b);
  if (!lengths.length) return 1200;
  // Median, so one unusually long gene does not squash all the others.
  return Math.max(200, lengths[Math.floor(lengths.length / 2)]);
}
/** Pixels per base pair for the current scale and zoom. */
function syntenyPxPerBp(scaleLen) {
  return (SYNTENY_TARGET_PX * syntenyZoom()) / Math.max(1, scaleLen);
}
function syntenyGeneWidth(gene, scaleLen) {
  const len = syntenyLength(gene);
  if (!len) return SYNTENY_MIN_GENE_PX * syntenyZoom();
  return Math.max(SYNTENY_MIN_GENE_PX, len * syntenyPxPerBp(scaleLen));
}
/* Rows are anchored on the focal gene so the comparison between alternative
   blocks is vertical. Laying every row out from the left edge only lined up
   while all blocks happened to carry the same number of upstream neighbours. */
function syntenyRowLayout(genes, left, gap, scaleLen) {
  const widths = genes.map(gene => syntenyGeneWidth(gene, scaleLen));
  const targetIdx = genes.findIndex(gene => gene.target);
  let offset = 0;
  if (targetIdx >= 0) {
    let before = 0;
    for (let i = 0; i < targetIdx; i += 1) before += widths[i] + gap;
    offset = before;
  }
  let x = left;
  const items = genes.map((gene, i) => {
    const item = {gene, x, width: widths[i]};
    x += widths[i] + gap;
    return item;
  });
  // Shift so every row's focal gene sits at the same distance from the left.
  const shift = offset;
  items.forEach(item => { item.x += (syntenyAnchor(scaleLen) - shift); });
  return items;
}
/** Distance from the row start to the focal gene, shared by every row. */
let syntenyAnchorPx = 0;
function syntenyAnchor() {
  return syntenyAnchorPx;
}
function arrowPath(x, y, w, h, strand = "+") {
  const head = Math.min(12, Math.max(7, w * 0.24));
  if (strand === "-") {
    return `M ${x + w} ${y} H ${x + head} L ${x} ${y + h / 2} L ${x + head} ${y + h} H ${x + w} Z`;
  }
  return `M ${x} ${y} H ${x + w - head} L ${x + w} ${y + h / 2} L ${x + w - head} ${y + h} H ${x} Z`;
}
/** Round a base-pair span to a readable 1/2/5 x 10^n step for the scale bar. */
function niceScaleBp(raw) {
  const magnitude = Math.pow(10, Math.floor(Math.log10(Math.max(1, raw))));
  const scaled = raw / magnitude;
  const step = scaled >= 5 ? 5 : scaled >= 2 ? 2 : 1;
  return step * magnitude;
}
function formatBp(bp) {
  if (bp >= 1000) return `${(bp / 1000).toLocaleString(undefined, {maximumFractionDigits: 1})} kb`;
  return `${Math.round(bp)} bp`;
}
function renderSyntenyPanel() {
  if (syntenyHidden) return;
  if (!syntenyPanel || !syntenySvg || !syntenySummary) return;
  const cluster = data.clusters[selected];
  const blocks = cluster && data.synteny ? (data.synteny[cluster.name] || []) : [];
  updateSelectedSyntenyCarriers(blocks);
  const scrollContainer = syntenyPanel.querySelector(".synteny-scroll") || syntenyPanel;
  const paneWidth = Math.max(360, scrollContainer.clientWidth || syntenyPanel.clientWidth || 600);
  const nextKey = cluster
    ? `${cluster.name}|${blocks.length}|${selectedSyntenyBlock}|${selectedSyntenyReference}|${syntenyZoom().toFixed(2)}|${paneWidth}`
    : "none";
  if (nextKey === syntenyKey) return;
  syntenyKey = nextKey;
  treeDirty = true;

  if (!cluster) {
    selectedSyntenyCarriers = new Set();
    syntenySummary.textContent = "No cluster selected.";
    syntenySvg.setAttribute("width", String(paneWidth));
    syntenySvg.setAttribute("height", "80");
    syntenySvg.innerHTML = `<text x="14" y="44" class="synteny-empty-text">Select a cluster to inspect flanking loci.</text>`;
    return;
  }
  if (!blocks.length) {
    selectedSyntenyCarriers = new Set();
    syntenySummary.textContent = `${cluster.name}: no 5-neighbour synteny blocks in this run.`;
    syntenySvg.setAttribute("width", String(paneWidth));
    syntenySvg.setAttribute("height", "80");
    syntenySvg.innerHTML = `<text x="14" y="44" class="synteny-empty-text">No unique synteny block data for ${escapeHtml(cluster.name)}.</text>`;
    return;
  }

  const rows = [];
  blocks.forEach((block, idx) => {
    const representativeGenes = blockGenes(block, cluster.name);
    rows.push({kind:"block", block, blockIndex:idx, genes:representativeGenes});
    if (idx === 0) {
      const referenceGenes = referenceGenesForBlock(block, cluster.name);
      if (syntenyIdentitySignature(referenceGenes) !== syntenyIdentitySignature(representativeGenes)) {
        rows.push({kind:"reference", block, blockIndex:idx, genes:referenceGenes});
      }
    }
  });
  const hasReferenceRow = rows.some(row => row.kind === "reference");
  if (selectedSyntenyReference && !hasReferenceRow) {
    selectedSyntenyReference = false;
    updateSelectedSyntenyCarriers(blocks);
  }

  // Layout constants
  const META_W = 168;
  const ROW_H = 78;
  const GENE_H = 22;
  const GAP = 12 * syntenyZoom();
  const TOP = 18;
  const PADDING_X = 12;
  const innerLeft = PADDING_X + META_W + 8;
  const lengthScale = syntenyLengthScale(rows);
  // Widest run of upstream neighbours across all rows: every row is shifted so
  // its focal gene lands here, which is what makes the rows comparable.
  syntenyAnchorPx = 0;
  rows.forEach(row => {
    const genes = row.genes || [];
    const targetIdx = genes.findIndex(gene => gene.target);
    if (targetIdx <= 0) return;
    let before = 0;
    for (let i = 0; i < targetIdx; i += 1) {
      before += syntenyGeneWidth(genes[i], lengthScale) + GAP;
    }
    syntenyAnchorPx = Math.max(syntenyAnchorPx, before);
  });
  const rowLayouts = rows.map(row => syntenyRowLayout(row.genes, innerLeft, GAP, lengthScale));
  const measuredWidth = Math.max(
    paneWidth,
    ...rowLayouts.map(layout => {
      const last = layout[layout.length - 1];
      return last ? last.x + last.width + 16 : innerLeft + 80;
    }),
  );
  const width = measuredWidth;
  const height = TOP + rows.length * ROW_H + 12;

  const parts = [];
  const selectedBlock = selectedSyntenyBlockData();
  const selectedCount = selectedSyntenyReference
    ? selectedSyntenyCarriers.size
    : selectedBlock ? syntenyCarrierIndices(selectedBlock).length : 0;
  const selectedLabel = selectedSyntenyReference
    ? `${data.referenceGenome ? data.referenceGenome.name : "reference"} context`
    : `#${selectedSyntenyBlock + 1}`;
  const referenceSummary = hasReferenceRow
    ? " Reference row inferred from shared neighbours."
    : "";
  const pxPerBp = syntenyPxPerBp(lengthScale);
  syntenySummary.textContent = `${cluster.name}: ${blocks.length} unique 5-neighbour block${blocks.length === 1 ? "" : "s"}; selected ${selectedLabel} marks ${selectedCount} tip${selectedCount === 1 ? "" : "s"}.${referenceSummary}`;

  // Collect per-row positions BEFORE drawing so we can draw connection lines on
  // a single layer that sits above row backgrounds but below the arrows.
  const rowPositions = [];
  rows.forEach((row, rowIdx) => {
    const layout = rowLayouts[rowIdx];
    const y = TOP + rowIdx * ROW_H + 14;
    const positions = new Map();
    let targetCenter = null;
    layout.forEach(item => {
      const gene = item.gene;
      const center = {x: item.x + item.width / 2, y: y + GENE_H / 2};
      const key = syntenyLinkKey(gene);
      if (key && !positions.has(key)) positions.set(key, center);
      if (gene.target) targetCenter = center;
    });
    rowPositions.push({positions, targetCenter, y});
  });

  // Scale bar: arrows are drawn to base-pair scale, so say what the scale is.
  const barBp = niceScaleBp(200 / pxPerBp);
  const barPx = barBp * pxPerBp;
  parts.push(
    `<line x1="${PADDING_X + 6}" y1="10" x2="${PADDING_X + 6 + barPx}" y2="10" class="synteny-scale-bar"/>` +
    `<line x1="${PADDING_X + 6}" y1="6" x2="${PADDING_X + 6}" y2="14" class="synteny-scale-bar"/>` +
    `<line x1="${PADDING_X + 6 + barPx}" y1="6" x2="${PADDING_X + 6 + barPx}" y2="14" class="synteny-scale-bar"/>` +
    `<text x="${PADDING_X + 12 + barPx}" y="13" class="synteny-scale-text">${formatBp(barBp)}</text>`
  );

  // Layer 1: row backgrounds.
  rows.forEach((row, rowIdx) => {
    const rowY = TOP + rowIdx * ROW_H;
    const isSelected = row.kind === "reference"
      ? selectedSyntenyReference
      : row.blockIndex === selectedSyntenyBlock && !selectedSyntenyReference;
    const classes = ["synteny-row-bg"];
    if (row.kind === "reference") classes.push("synteny-reference-bg");
    if (isSelected) classes.push("synteny-row-selected");
    parts.push(
      `<rect x="${PADDING_X}" y="${rowY}" width="${width - PADDING_X * 2}" height="${ROW_H - 4}" rx="10" class="${classes.join(" ")}"/>`
    );
  });

  // Layer 2: connection lines between adjacent rows for shared neighbour clusters.
  rowPositions.forEach((rp, rowIdx) => {
    if (rowIdx === 0) return;
    const prev = rowPositions[rowIdx - 1];
    rp.positions.forEach((pos, name) => {
      const previous = prev.positions.get(name);
      if (!previous) return;
      const midY = (previous.y + pos.y) / 2;
      const label = escapeHtml(syntenyGeneLabel(name));
      parts.push(
        `<path d="M ${previous.x} ${previous.y + GENE_H / 2 - 2} C ${previous.x} ${midY} ${pos.x} ${midY} ${pos.x} ${pos.y - GENE_H / 2 + 2}" class="synteny-link" data-link-cluster="${escapeHtml(name)}"><title>shared neighbour: ${label}</title></path>`
      );
    });
  });

  // Layer 3: meta column dividers + content per row.
  rows.forEach((row, rowIdx) => {
    const rowY = TOP + rowIdx * ROW_H;
    const baseY = rowY + 14;
    const block = row.block;
    const isSelected = row.kind === "reference"
      ? selectedSyntenyReference
      : row.blockIndex === selectedSyntenyBlock && !selectedSyntenyReference;
    const rowAttrs = row.kind === "block"
      ? ` class="synteny-row" data-block-index="${row.blockIndex}"`
      : ` class="synteny-reference-row" data-reference-context="1"`;
    parts.push(`<g${rowAttrs}>`);
    // Hit area (transparent rectangle the user can click to select the row)
    parts.push(
      `<rect x="${PADDING_X}" y="${rowY}" width="${META_W}" height="${ROW_H - 4}" class="synteny-meta-hit" rx="10"/>`
    );

    // Meta column
    if (row.kind === "reference") {
      parts.push(`<text x="${PADDING_X + 12}" y="${rowY + 18}" class="synteny-rank">REF</text>`);
      parts.push(`<text x="${PADDING_X + 12}" y="${rowY + 36}" class="synteny-meta-title">H37Rv context</text>`);
      parts.push(
        `<text x="${PADDING_X + 12}" y="${rowY + 52}" class="synteny-meta-sub">inferred from shared neighbours</text>`
      );
      const refText = selectedSyntenyReference ? "tip marked" : "click row to mark";
      parts.push(`<text x="${PADDING_X + 12}" y="${rowY + 68}" class="synteny-meta-sub">${escapeHtml(refText)}</text>`);
    } else {
      parts.push(`<text x="${PADDING_X + 12}" y="${rowY + 18}" class="synteny-rank">#${block.rank || row.blockIndex + 1}</text>`);
      const fraction = block.fraction === null || block.fraction === undefined || block.fraction === "NA"
        ? ""
        : ` · ${Math.round(Number(block.fraction) * 100)}%`;
      parts.push(
        `<text x="${PADDING_X + 12 + 26}" y="${rowY + 18}" class="synteny-meta-rank-sub">n=${block.nGenomes ?? "?"}${fraction}</text>`
      );
      const titleText = block.representativeGenome
        ? `${block.representativeGenome}`
        : (firstExampleLocus(block) || cluster.name);
      parts.push(`<text x="${PADDING_X + 12}" y="${rowY + 36}" class="synteny-meta-title">${escapeHtml(clipText(titleText, 22))}</text>`);
      const coordText = block.representativeContig && block.representativeStart && block.representativeEnd
        ? `${block.representativeContig}:${formatCoord(block.representativeStart)}-${formatCoord(block.representativeEnd)}`
        : firstExampleLocus(block) || "";
      const lenSuffix = block.representativeLength ? ` · ${block.representativeLength} bp` : "";
      parts.push(
        `<text x="${PADDING_X + 12}" y="${rowY + 52}" class="synteny-meta-sub">${escapeHtml(clipText(coordText + lenSuffix, 30))}</text>`
      );
      const interp = (block.interpretation || "").replace(/_/g, " ");
      const selectedTag = isSelected ? "selected" : (interp || "click to mark tips");
      parts.push(`<text x="${PADDING_X + 12}" y="${rowY + 68}" class="synteny-meta-sub">${escapeHtml(clipText(selectedTag, 26))}</text>`);
    }
    // Vertical divider after meta column
    parts.push(
      `<line x1="${PADDING_X + META_W + 1}" y1="${rowY + 6}" x2="${PADDING_X + META_W + 1}" y2="${rowY + ROW_H - 10}" class="synteny-divider"/>`
    );

    // Gene arrows
    const layout = rowLayouts[rowIdx];
    layout.forEach(item => {
      const gene = item.gene;
      const x = item.x;
      const geneW = item.width;
      const cls = ["synteny-gene"];
      if (gene.gap) cls.push("synteny-gap");
      else if (gene.target) cls.push("synteny-target");
      else if (row.kind === "reference") cls.push("synteny-reference-gene");
      else cls.push("synteny-neighbor");
      const label = syntenyDisplayLabel(gene);
      const coord = syntenyCoordLabel(gene);
      const title = [syntenyGeneLabel(gene.cluster || gene.name), gene.geneId, coord].filter(Boolean).join(" | ");
      const geneAttrs = gene.gap || !gene.cluster || gene.cluster === "NA"
        ? ""
        : ` data-gene-cluster="${escapeHtml(gene.cluster)}"`;
      parts.push(
        `<path d="${arrowPath(x, baseY, geneW, GENE_H, gene.strand)}" class="${cls.join(" ")}"${geneAttrs}><title>${escapeHtml(title)}</title></path>`
      );
      const labelClass = gene.target ? "synteny-target-label" : "synteny-arrow-label";
      parts.push(
        `<text x="${x + geneW / 2}" y="${baseY + 14}" text-anchor="middle" class="${labelClass}">${escapeHtml(clipText(label, 11))}</text>`
      );
      const subLabel = clipText(syntenySecondaryLabel(gene), 12);
      if (subLabel) {
        parts.push(
          `<text x="${x + geneW / 2}" y="${baseY + GENE_H + 12}" text-anchor="middle" class="synteny-sub-label">${escapeHtml(subLabel)}</text>`
        );
      }
      if (gene.length) {
        parts.push(
          `<text x="${x + geneW / 2}" y="${baseY + GENE_H + 24}" text-anchor="middle" class="synteny-len-label">${escapeHtml(gene.length)} bp</text>`
        );
      }
    });
    parts.push("</g>");
  });

  syntenySvg.setAttribute("width", String(width));
  syntenySvg.setAttribute("height", String(height));
  syntenySvg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  syntenySvg.innerHTML = parts.join("");
  fitSyntenyHeight(height);
  // Each row opens on its meta column: rank, carriers, representative genome and
  // coordinates. Centring the focal gene, which is what a narrow side panel had
  // to do, pushes that column off the left edge now the band is full width, so
  // only scroll when the focal gene would otherwise be out of view.
  if (scrollContainer) {
    scrollContainer.scrollLeft = 0;
    scrollContainer.scrollTop = 0;
    const firstTarget = (rowLayouts[0] || []).find(item => item.gene.target);
    if (firstTarget) {
      const targetRight = firstTarget.x + firstTarget.width;
      const visible = scrollContainer.clientWidth;
      if (visible > 0 && targetRight > visible) {
        scrollContainer.scrollLeft = Math.max(0, targetRight - visible + 80);
      }
    }
  }
}
function formatCoord(value) {
  const n = Number(value);
  if (!Number.isFinite(n)) return String(value);
  return n.toLocaleString();
}
function loadLayout() {
  try {
    const raw = localStorage.getItem(LAYOUT_KEY);
    if (!raw) return null;
    return JSON.parse(raw);
  } catch (_e) {
    return null;
  }
}
function persistLayout() {
  try {
    const payload = {
      treeW: Math.round(layoutState.treeW),
      syntenyH: Math.round(layoutState.syntenyH),
      inspectorW: Math.round(layoutState.inspectorW),
      syntenyUserSized: layoutState.syntenyUserSized,
      syntenyHidden,
    };
    localStorage.setItem(LAYOUT_KEY, JSON.stringify(payload));
  } catch (_e) {
    /* localStorage unavailable */
  }
}
const layoutState = (() => {
  const stored = loadLayout() || {};
  return {
    treeW: Number.isFinite(stored.treeW) ? Math.max(MIN_TREE_W, stored.treeW) : DEFAULT_TREE_W,
    syntenyH: Number.isFinite(stored.syntenyH)
      ? Math.max(MIN_SYNTENY_H, stored.syntenyH)
      : DEFAULT_SYNTENY_H,
    inspectorW: Number.isFinite(stored.inspectorW)
      ? Math.max(MIN_INSPECTOR_W, stored.inspectorW)
      : DEFAULT_INSPECTOR_W,
    syntenyUserSized: stored.syntenyUserSized === true,
    syntenyHiddenStored: stored.syntenyHidden,
  };
})();
function applyLayout() {
  if (!viewerSplit || !viewerBody || !explorerSection) return;
  // Below this the shell has not been laid out yet, and reacting would shrink
  // the panes to their minimums and persist that.
  const bodyW = viewerBody.clientWidth || 0;
  const measured = bodyW >= MIN_TREE_W + MIN_MATRIX_W;
  let treeW = layoutState.treeW;
  let inspectorW = layoutState.inspectorW;
  const inspectorShown = !!inspectorPane && getComputedStyle(inspectorPane).display !== "none";
  if (measured) {
    // The matrix keeps its floor; the tree and the inspector give way first.
    const used = treeW + (inspectorShown ? inspectorW : 0) + 12;
    const availableForMatrix = bodyW - used;
    if (availableForMatrix < MIN_MATRIX_W) {
      let overshoot = MIN_MATRIX_W - availableForMatrix;
      const shrinkTree = Math.min(treeW - MIN_TREE_W, overshoot * 0.5);
      treeW = Math.max(MIN_TREE_W, treeW - shrinkTree);
      overshoot -= shrinkTree;
      if (inspectorShown && overshoot > 0) {
        inspectorW = Math.max(MIN_INSPECTOR_W, inspectorW - overshoot);
      }
    }
    layoutState.treeW = treeW;
    layoutState.inspectorW = inspectorW;
  }
  // Synteny is a horizontal band now, so its size is a height bounded by the
  // workspace: gene arrows need the full width, not a slice of it.
  const workspaceH = explorerSection.clientHeight || 0;
  let syntenyH = layoutState.syntenyH;
  if (workspaceH > 0) {
    syntenyH = Math.max(
      MIN_SYNTENY_H,
      Math.min(syntenyH, Math.round(workspaceH * MAX_SYNTENY_FRACTION)),
    );
    layoutState.syntenyH = syntenyH;
  }
  explorerSection.style.setProperty("--synteny-h", `${Math.round(syntenyH)}px`);
  viewerBody.style.setProperty("--tree-w", `${treeW}px`);
  viewerBody.style.setProperty("--inspector-w", `${inspectorW}px`);
  viewerSplit.style.setProperty("--tree-w", `${treeW}px`);
  viewerSplit.style.setProperty("--inspector-w", `${inspectorW}px`);
  treeWidth = Math.max(160, Math.round(treeW) - LINEAGE_STRIP_W - 4);
}
/* Two synteny rows are the point of the panel, so open at the height that shows
   them rather than at a number picked in advance. Once the reader drags the
   seam, that is their decision and the content stops overriding it. */
function fitSyntenyHeight(contentHeight) {
  if (layoutState.syntenyUserSized || syntenyHidden || !explorerSection) return;
  const headH = syntenyPanel ? (syntenyPanel.querySelector(".pane-head")?.offsetHeight || 34) : 34;
  const workspaceH = explorerSection.clientHeight || 0;
  const ceiling = workspaceH > 0 ? Math.round(workspaceH * MAX_SYNTENY_FRACTION) : 520;
  const want = Math.max(MIN_SYNTENY_H, Math.min(ceiling, Math.round(contentHeight) + headH + 10));
  if (Math.abs(want - layoutState.syntenyH) < 2) return;
  layoutState.syntenyH = want;
  applyLayout();
}
function applyViewerSplit() {
  applyLayout();
}
function applySyntenyVisibility() {
  if (!explorerSection || !syntenyToggle) return;
  explorerSection.classList.toggle("synteny-hidden", syntenyHidden);
  syntenyToggle.classList.toggle("active", !syntenyHidden);
  syntenyToggle.setAttribute(
    "data-tip",
    syntenyHidden
      ? "Show the flanking-gene context for the selected cluster (S)"
      : "Hide the synteny panel (S)",
  );
  syntenyToggle.setAttribute("aria-pressed", syntenyHidden ? "false" : "true");
  applyLayout();
  persistLayout();
  if (!syntenyHidden) {
    syntenyKey = "";
    renderSyntenyPanel();
  } else {
    selectedSyntenyCarriers = new Set();
    treeDirty = true;
  }
}
/* Three resizable seams: tree|matrix and matrix|inspector move horizontally,
   matrix|synteny moves vertically. */
function setupDragHandle(handle, side) {
  if (!handle) return;
  const vertical = side === "synteny";
  let start = 0;
  let startTreeW = 0;
  let startInspectorW = 0;
  let startSyntenyH = 0;
  let bodyW = 0;
  let active = false;
  const applyDelta = delta => {
    if (side === "tree") {
      const maxTree = Math.max(MIN_TREE_W, bodyW - MIN_MATRIX_W - startInspectorW - 12);
      layoutState.treeW = Math.max(MIN_TREE_W, Math.min(maxTree, startTreeW + delta));
    } else if (side === "inspector") {
      const maxInspector = Math.max(MIN_INSPECTOR_W, bodyW - MIN_MATRIX_W - layoutState.treeW - 12);
      layoutState.inspectorW = Math.max(
        MIN_INSPECTOR_W,
        Math.min(maxInspector, startInspectorW - delta),
      );
    } else {
      // Dragging up grows the synteny band, so the delta is inverted. Once the
      // reader has picked a height, stop second-guessing it from the content.
      layoutState.syntenyH = Math.max(MIN_SYNTENY_H, startSyntenyH - delta);
      layoutState.syntenyUserSized = true;
    }
    applyLayout();
  };
  const capture = () => {
    bodyW = viewerBody ? viewerBody.clientWidth || 0 : 0;
    startTreeW = layoutState.treeW;
    startInspectorW = layoutState.inspectorW;
    startSyntenyH = layoutState.syntenyH;
  };
  const onMove = event => {
    if (!active) return;
    applyDelta((vertical ? event.clientY : event.clientX) - start);
    scheduleRender(false);
  };
  const onUp = event => {
    if (!active) return;
    active = false;
    handle.classList.remove("dragging");
    document.body.classList.remove("recore-dragging", "recore-dragging-row", "recore-dragging-col");
    persistLayout();
    if (handle.hasPointerCapture && handle.hasPointerCapture(event.pointerId)) {
      handle.releasePointerCapture(event.pointerId);
    }
    treeDirty = true;
    resetMatrixCache();
    syntenyKey = "";
    render(false);
  };
  handle.addEventListener("pointerdown", event => {
    if (event.button !== 0) return;
    capture();
    active = true;
    start = vertical ? event.clientY : event.clientX;
    handle.classList.add("dragging");
    document.body.classList.add("recore-dragging");
    document.body.classList.add(vertical ? "recore-dragging-row" : "recore-dragging-col");
    try { handle.setPointerCapture(event.pointerId); } catch (_e) {}
    event.preventDefault();
  });
  handle.addEventListener("pointermove", onMove);
  handle.addEventListener("pointerup", onUp);
  handle.addEventListener("pointercancel", onUp);
  handle.addEventListener("keydown", event => {
    const keys = vertical ? ["ArrowUp", "ArrowDown"] : ["ArrowLeft", "ArrowRight"];
    if (!keys.includes(event.key)) return;
    const step = event.shiftKey ? 32 : 12;
    capture();
    applyDelta((event.key === keys[1] ? 1 : -1) * step);
    persistLayout();
    treeDirty = true;
    resetMatrixCache();
    syntenyKey = "";
    render(false);
    event.preventDefault();
  });
}
function treeCaption() {
  const rootLabel = hasL8Tip(treeRoot) ? "L8-rooted" : "input-rooted";
  return treeMode === "phylogram" && maxDistance > 0
    ? `phylogram, ${rootLabel}, longest root-to-tip path ${formatBranchLength(maxDistance)}`
    : `cladogram, ${rootLabel}, branch lengths ignored`;
}
function updateTopBarMeta() {
  if (topTreeMeta) {
    // Short in the header, full text on hover: the caption is the first thing
    // that should give way when the pane narrows.
    topTreeMeta.textContent = `${leafOrder.length} tips`;
    topTreeMeta.setAttribute("data-tip", `${leafOrder.length} matched tips · ${treeCaption()}`);
  }
  if (topMatrixMeta) topMatrixMeta.textContent = `${visibleClusters.length} / ${data.clusters.length} clusters`;
  if (topSyntenyMeta) {
    const cluster = data.clusters[selected];
    const blocks = cluster && data.synteny ? (data.synteny[cluster.name] || []) : [];
    topSyntenyMeta.textContent = cluster
      ? (blocks.length ? `${cluster.name} – ${blocks.length} unique block${blocks.length === 1 ? "" : "s"}` : `${cluster.name} – no synteny data`)
      : "select a cluster";
  }
}
function referenceMechanismRow(cluster) {
  const referenceCell = (cluster.referenceCell || "").trim();
  const targetGene = (cluster.targetReferenceGene || "").trim();
  const targetCluster = (cluster.targetReferenceCluster || "").trim();
  const eventClass = (cluster.referenceEventClass || "").trim();
  const refState = (cluster.referenceState || "").trim();
  const bits = [];
  if (referenceCell) bits.push(`Reference gene(s): <strong>${escapeHtml(referenceCell)}</strong>`);
  if (targetGene || targetCluster || eventClass) {
    const source = [];
    if (targetGene) source.push(`source gene <strong>${escapeHtml(targetGene)}</strong>`);
    if (targetCluster) source.push(`source cluster <strong>${escapeHtml(targetCluster)}</strong>`);
    if (eventClass) source.push(escapeHtml(eventClass));
    bits.push(`Event source: ${source.join("; ")}`);
  }
  if (refState) {
    // referenceState is the single-letter cell code; spell it out rather than
    // making the reader look up what "P" means.
    const spelled = labels[refState] || refState;
    bits.push(`State in reference: <strong>${escapeHtml(spelled)}</strong>`);
  }
  if (!bits.length) {
    const mechanism = processKind(cluster);
    const text = mechanism === "clean_insertion_or_acquisition_candidate"
      ? "No reference equivalent detected; candidate de novo/reference-absent locus."
      : "No reference equivalent detected by Recore.";
    bits.push(text);
  }
  return `<dt>Reference (audit)</dt><dd>${bits.join("<br>")}</dd>`;
}
/* Canvases are laid out in CSS pixels but must be allocated at the device ratio,
   or a HiDPI screen upscales a 1x bitmap and the cells come out soft. The cap
   still applies, though the viewport-sized canvases are nowhere near it now. */
function canvasScale(cssW, cssH) {
  const dpr = Math.max(1, Math.min(3, window.devicePixelRatio || 1));
  const byDim = maxCanvasDim / Math.max(1, cssW, cssH);
  const byArea = Math.sqrt(maxCanvasPixels / Math.max(1, cssW * cssH));
  return Math.max(1, Math.min(dpr, byDim, byArea));
}
function sizeCanvas(target, context, cssW, cssH, scale) {
  target.width = Math.max(1, Math.round(cssW * scale));
  target.height = Math.max(1, Math.round(cssH * scale));
  target.style.width = `${cssW}px`;
  target.style.height = `${cssH}px`;
  context.setTransform(scale, 0, 0, scale, 0, 0);
}

/* ---- Virtual geometry ----------------------------------------------------
   The matrix is as large as clusters x genomes, which no canvas can hold: the
   browser caps a 2D canvas at 16384px per dimension, so 30,000 columns cannot
   fit even at one pixel each. Every pane therefore draws only its visible
   window onto a viewport-sized canvas, and a spacer carries the virtual size so
   the browser keeps owning the scrollbars. Cell sizes may be fractional; below
   one pixel per cell the drawing collapses onto shared pixels, which is what
   lets a 300 million cell matrix render at all. */
function updateCellGeometry() {
  const rows = leafOrder.length;
  const cols = visibleClusters.length;
  cellH = currentRowHeight();
  cellW = Math.max(MIN_CELL_PX, currentCellSize() * 0.82);
  virtualW = Math.max(1, Math.ceil(cols * cellW));
  virtualH = Math.max(1, Math.ceil(rows * cellH));
  treeHeight = virtualH;
  if (matrixSpacer) {
    matrixSpacer.style.width = `${virtualW}px`;
    matrixSpacer.style.height = `${virtualH}px`;
  }
  if (headerSpacer) {
    headerSpacer.style.width = `${virtualW}px`;
    headerSpacer.style.height = "1px";
  }
  if (treeSpacer) {
    treeSpacer.style.width = "1px";
    treeSpacer.style.height = `${virtualH}px`;
  }
  return {rows, cols};
}
function matrixScrollLeft() {
  return matrixPane ? matrixPane.scrollLeft : 0;
}
function matrixScrollTop() {
  return matrixPane ? matrixPane.scrollTop : 0;
}
/** Half-open range of columns touching the visible window. */
function visibleColumnRange(width) {
  const left = matrixScrollLeft();
  const first = Math.max(0, Math.floor(left / cellW));
  const last = Math.min(visibleClusters.length, Math.ceil((left + width) / cellW) + 1);
  return [first, last];
}
/** Half-open range of rows touching the visible window. */
function visibleRowRange(height) {
  const top = matrixScrollTop();
  const first = Math.max(0, Math.floor(top / cellH));
  const last = Math.min(leafOrder.length, Math.ceil((top + height) / cellH) + 1);
  return [first, last];
}
function viewportSize(element) {
  if (!element) return {w: 0, h: 0};
  return {w: Math.max(1, element.clientWidth), h: Math.max(1, element.clientHeight)};
}

function renderMatrix(scrollToSelected) {
  const {rows, cols} = updateCellGeometry();
  if (scrollToSelected) {
    const selectedCol = selectedVisibleColumn();
    if (selectedCol >= 0 && matrixPane) {
      const view = viewportSize(matrixViewport);
      const x = selectedCol * cellW;
      if (x < matrixPane.scrollLeft || x > matrixPane.scrollLeft + view.w - cellW) {
        matrixPane.scrollLeft = Math.max(0, x - view.w / 2);
      }
    }
  }
  const view = viewportSize(matrixViewport);
  const scale = canvasScale(view.w, view.h);
  sizeCanvas(canvas, ctx, view.w, view.h, scale);
  ctx.fillStyle = panelBg;
  ctx.fillRect(0, 0, view.w, view.h);

  const left = matrixScrollLeft();
  const top = matrixScrollTop();
  const [firstCol, lastCol] = visibleColumnRange(view.w);
  const [firstRow, lastRow] = visibleRowRange(view.h);
  // Cells never shrink below a hairline, so a dense matrix stays legible instead
  // of fading out when several clusters share a pixel column.
  const drawW = Math.max(1, cellW > 2 ? cellW - 1 : cellW);
  const drawH = Math.max(1, cellH > 2 ? cellH - 1 : cellH);
  const gutterRows = cellH >= 2;
  const byProcess = matrixColorMode.value === "process";
  const spanTop = firstRow * cellH - top;
  const spanHeight = (lastRow - firstRow) * cellH;

  for (let c = firstCol; c < lastCol; c += 1) {
    const clusterIdx = originalClusterAtVisibleColumn(c);
    if (clusterIdx < 0) continue;
    const sets = clusterSets(clusterIdx);
    const x = c * cellW - left;
    const processFill = byProcess
      ? processColor(processKind(data.clusters[clusterIdx]))
      : null;
    const fillWholeColumn = () => {
      if (gutterRows) {
        for (let r = firstRow; r < lastRow; r += 1) {
          ctx.fillRect(x, r * cellH - top, drawW, drawH);
        }
      } else {
        ctx.fillRect(x, spanTop, drawW, spanHeight);
      }
    };
    for (let state = 0; state < 3; state += 1) {
      const set = sets[state];
      ctx.fillStyle = processFill || colors[STATE_CODES[state]] || colors.A;
      if (set.kind === "all") {
        fillWholeColumn();
        continue;
      }
      if (set.kind === "complement") {
        fillWholeColumn();
        ctx.fillStyle = panelBg;
        for (let i = 0; i < set.idx.length; i += 1) {
          const row = rowOfGenome[set.idx[i]];
          if (row >= firstRow && row < lastRow) {
            ctx.fillRect(x, row * cellH - top, drawW, Math.max(1, cellH));
          }
        }
        continue;
      }
      // Only the genomes inside the visible band are worth walking.
      const idx = set.idx;
      for (let i = 0; i < idx.length; i += 1) {
        const row = rowOfGenome[idx[i]];
        if (row >= firstRow && row < lastRow) {
          ctx.fillRect(x, row * cellH - top, drawW, drawH);
        }
      }
    }
  }

  const bandH = Math.max(0, Math.min(view.h, virtualH - top));
  const selectedCol = selectedVisibleColumn();
  if (selectedCol >= firstCol && selectedCol < lastCol && bandH > 0) {
    const x = selectedCol * cellW - left;
    ctx.fillStyle = "rgba(37,99,235,.13)";
    ctx.fillRect(x, 0, Math.max(1, cellW), bandH);
    ctx.strokeStyle = selectStroke;
    ctx.lineWidth = 1;
    ctx.strokeRect(x + 0.5, 0.5, Math.max(1, cellW - 1), bandH - 1);
  }
  if (hoverRow >= firstRow && hoverRow < lastRow) {
    ctx.fillStyle = "rgba(245,158,11,.14)";
    ctx.fillRect(0, hoverRow * cellH - top, Math.min(view.w, virtualW - left), Math.max(1, cellH));
  }
  if (hoverCol >= firstCol && hoverCol < lastCol && bandH > 0) {
    ctx.fillStyle = "rgba(245,158,11,.10)";
    ctx.fillRect(hoverCol * cellW - left, 0, Math.max(1, cellW), bandH);
  }
  if (rows === 0 || cols === 0) {
    ctx.fillStyle = labelInk;
    ctx.font = LABEL_FONT;
    ctx.fillText("No clusters match the current filters.", 14, 22);
  }
}

function columnLabelsFlat() {
  return cellW >= LABEL_HORIZONTAL_CELL_W;
}
function columnLabelBand() {
  if (cellW < LABEL_MIN_CELL_W || visibleClusters.length > 4000) return 0;
  return columnLabelsFlat() ? LABEL_BAND_H_FLAT : LABEL_BAND_H;
}
function renderColumnHeader() {
  const cols = visibleClusters.length;
  const band = columnLabelBand();
  const cssH = headerH + band;
  const view = viewportSize(headerScroller);
  const width = Math.max(1, view.w);
  sizeCanvas(headerCanvas, headerCtx, width, cssH, canvasScale(width, cssH));
  const headerWrap = headerCanvas.parentElement;
  if (headerWrap) headerWrap.style.height = `${cssH}px`;
  // Keep the tree pane's first row level with the matrix's first row.
  if (treeHeadSpacer) treeHeadSpacer.style.height = `${cssH}px`;
  headerCtx.fillStyle = headerBg;
  headerCtx.fillRect(0, 0, width, cssH);

  const left = matrixScrollLeft();
  const [firstCol, lastCol] = visibleColumnRange(width);
  const selectedCol = selectedVisibleColumn();
  const byProcess = matrixColorMode.value === "process";
  const tickW = Math.max(1, cellW > 2 ? cellW - 1 : cellW);

  for (let c = firstCol; c < lastCol; c += 1) {
    const clusterIdx = originalClusterAtVisibleColumn(c);
    if (clusterIdx < 0) continue;
    const x = c * cellW - left;
    headerCtx.fillStyle = byProcess
      ? processColor(processKind(data.clusters[clusterIdx]))
      : gridStroke;
    headerCtx.fillRect(x, cssH - 9, tickW, 6);
    if (byProcess && band) {
      headerCtx.fillRect(x, band, tickW, 6);
    }
  }
  if (hoverCol >= firstCol && hoverCol < lastCol) {
    headerCtx.fillStyle = hoverColTint;
    headerCtx.fillRect(hoverCol * cellW - left, 0, Math.max(1, cellW), cssH);
  }
  if (selectedCol >= firstCol && selectedCol < lastCol) {
    const x = selectedCol * cellW - left;
    headerCtx.fillStyle = selectColTint;
    headerCtx.fillRect(x, 0, Math.max(1, cellW), cssH);
    headerCtx.strokeStyle = selectStroke;
    headerCtx.strokeRect(x + 0.5, 0.5, Math.max(1, cellW - 1), cssH - 1);
  }
  if (!band) return;
  const flat = columnLabelsFlat();
  headerCtx.save();
  headerCtx.font = LABEL_FONT;
  headerCtx.textBaseline = "middle";
  // Rotated: after rotate(-90deg) the local +x axis points up the screen, so
  // left-aligned text grows upward from the tick row and reads bottom to top.
  headerCtx.textAlign = flat ? "center" : "left";
  for (let c = firstCol; c < lastCol; c += 1) {
    const cluster = data.clusters[originalClusterAtVisibleColumn(c)];
    if (!cluster) continue;
    headerCtx.fillStyle = c === selectedCol ? selectStroke : labelInk;
    if (flat) {
      headerCtx.fillText(
        fitLabel(cluster.name, cellW - 8),
        c * cellW - left + cellW / 2,
        band / 2,
      );
      continue;
    }
    headerCtx.save();
    headerCtx.translate(c * cellW - left + cellW / 2, band - 6);
    headerCtx.rotate(-Math.PI / 2);
    headerCtx.fillText(fitLabel(cluster.name, band - 12), 0, 0);
    headerCtx.restore();
  }
  headerCtx.restore();
}
function fitLabel(text, maxWidth) {
  const value = String(text || "");
  if (headerCtx.measureText(value).width <= maxWidth) return value;
  for (let len = value.length - 1; len > 1; len -= 1) {
    const candidate = `${value.slice(0, len)}…`;
    if (headerCtx.measureText(candidate).width <= maxWidth) return candidate;
  }
  return "…";
}

/* The tree is a canvas rather than SVG: 10,000 tips is roughly 40,000 SVG nodes,
   which the browser will not lay out at an interactive frame rate. Only the tips
   inside the visible band are drawn, and internal edges are clipped to it. */
/* ---- Rerooting ------------------------------------------------------------
   The displayed root drives every ancestral reading in the panel, so it has to
   be the reader's choice rather than whatever the Newick happened to carry.
   assignTreeIds already yields an adjacency list, so rooting anywhere is a
   matter of picking the node and rebuilding downwards from it. */
function treeGraph(node) {
  const root = cloneNode(node);
  const nodes = [];
  const edges = [];
  assignTreeIds(root, null, nodes, edges);
  const adj = nodes.map(() => []);
  edges.forEach(([a, b, length]) => {
    adj[a].push({id: b, length: Number(length) || 0});
    adj[b].push({id: a, length: Number(length) || 0});
  });
  return {root, nodes, adj};
}
function buildFrom(nodes, adj, id, parentId, length) {
  const source = nodes[id];
  const children = adj[id]
    .filter(edge => edge.id !== parentId)
    .map(edge => buildFrom(nodes, adj, edge.id, id, edge.length));
  children.sort((a, b) => countTips(b) - countTips(a));
  return {label: source.label || "", genome: source.genome, length, children};
}
/** Farthest node from a start, by summed branch length. */
function farthestFrom(adj, start) {
  const dist = new Array(adj.length).fill(-1);
  const prev = new Array(adj.length).fill(-1);
  const stack = [start];
  dist[start] = 0;
  let best = start;
  while (stack.length) {
    const id = stack.pop();
    for (const edge of adj[id]) {
      if (dist[edge.id] >= 0) continue;
      dist[edge.id] = dist[id] + edge.length;
      prev[edge.id] = id;
      if (dist[edge.id] > dist[best]) best = edge.id;
      stack.push(edge.id);
    }
  }
  return {best, dist, prev};
}
/* Midpoint rooting: the two most distant tips define the tree's diameter, and
   the root goes halfway along the path between them. */
function midpointRoot(nodes, adj) {
  if (!adj.length) return 0;
  const first = farthestFrom(adj, 0).best;
  const {best: second, dist, prev} = farthestFrom(adj, first);
  const half = dist[second] / 2;
  let cursor = second;
  while (prev[cursor] >= 0 && dist[cursor] > half) {
    const parent = prev[cursor];
    // Stop at whichever end of the straddling edge is nearer the midpoint.
    if (dist[parent] <= half) {
      return half - dist[parent] <= dist[cursor] - half ? parent : cursor;
    }
    cursor = parent;
  }
  return cursor;
}
function rerootTree(mode) {
  if (mode === "input") return cloneNode(data.tree);
  if (mode === "l8") return rootByL8(data.tree);
  const {nodes, adj} = treeGraph(data.tree);
  if (mode === "midpoint") return buildFrom(nodes, adj, midpointRoot(nodes, adj), null, 0);
  if (mode.startsWith("genome:")) {
    const genomeIdx = Number(mode.slice(7));
    const id = nodes.findIndex(n => n.genome === genomeIdx);
    if (id >= 0) {
      // Root on the tip's parent, so the chosen genome hangs off the root as an
      // outgroup rather than becoming a degenerate root itself.
      const parent = adj[id][0];
      return buildFrom(nodes, adj, parent ? parent.id : id, null, 0);
    }
  }
  return rootByL8(data.tree);
}
/** Rebuild layout state after the displayed root changes. */
function applyRoot(mode) {
  treeRoot = rerootTree(mode);
  leafOrder = [];
  maxDepth = 0;
  maxDistance = 0;
  layout(treeRoot, 0, 0);
  buildRowOfGenome();
  treeDirty = true;
  resetMatrixCache();
  syntenyKey = "";
}
function renderTree() {
  if (!treeCanvas || !treeCtx) return;
  const view = viewportSize(treeViewport);
  const width = Math.max(1, view.w - LINEAGE_STRIP_W);
  treeWidth = width;
  sizeCanvas(treeCanvas, treeCtx, width, view.h, canvasScale(width, view.h));
  treeCtx.fillStyle = panelBg;
  treeCtx.fillRect(0, 0, width, view.h);

  const top = matrixScrollTop();
  const [firstRow, lastRow] = visibleRowRange(view.h);
  const showLabels = labelToggle.checked && cellH >= 5 && leafOrder.length <= 4000;

  treeCtx.lineWidth = 1;
  treeCtx.strokeStyle = gridStroke;
  drawTreeNode(treeRoot, top, firstRow, lastRow, view.h, showLabels);

  if (hoverRow >= firstRow && hoverRow < lastRow) {
    const y = hoverRow * cellH + cellH / 2 - top;
    treeCtx.strokeStyle = "rgba(249,115,22,.55)";
    treeCtx.beginPath();
    treeCtx.moveTo(0, y);
    treeCtx.lineTo(width, y);
    treeCtx.stroke();
  }
  renderLineageStrip(view.h, top, firstRow, lastRow);
}
function drawTreeNode(node, top, firstRow, lastRow, viewH, showLabels) {
  const rowSpan = nodeRowSpan(node);
  // Whole subtrees outside the band cost nothing.
  if (rowSpan.last < firstRow || rowSpan.first > lastRow) return;
  const x = xForNode(node);
  const y = node.row * cellH + cellH / 2 - top;
  if (node.children && node.children.length) {
    for (const child of node.children) {
      const childSpan = nodeRowSpan(child);
      if (childSpan.last < firstRow || childSpan.first > lastRow) continue;
      const cx = xForNode(child);
      const cy = child.row * cellH + cellH / 2 - top;
      treeCtx.strokeStyle = treeEdgeColor(child);
      treeCtx.beginPath();
      treeCtx.moveTo(x, cy);
      treeCtx.lineTo(cx, cy);
      treeCtx.moveTo(x, y);
      treeCtx.lineTo(x, cy);
      treeCtx.stroke();
      drawTreeNode(child, top, firstRow, lastRow, viewH, showLabels);
    }
    return;
  }
  if (node.row < firstRow || node.row >= lastRow) return;
  const style = tipStyle(node.genome);
  const radius = Math.max(0.5, Math.min(4.5, cellH / 2 - 0.5));
  treeCtx.fillStyle = style.fill;
  treeCtx.beginPath();
  treeCtx.arc(x, y, radius, 0, Math.PI * 2);
  treeCtx.fill();
  if (radius >= 1.5) {
    treeCtx.strokeStyle = style.stroke;
    treeCtx.lineWidth = 1;
    treeCtx.stroke();
  }
  if (highlightedRows.has(node.row)) {
    treeCtx.strokeStyle = "#0ea5e9";
    treeCtx.lineWidth = 1.6;
    treeCtx.beginPath();
    treeCtx.arc(x, y, radius + Math.max(2, Math.min(4, cellH * 0.4)), 0, Math.PI * 2);
    treeCtx.stroke();
  }
  if (isSelectedSyntenyCarrier(node.genome)) {
    treeCtx.strokeStyle = "#f97316";
    treeCtx.lineWidth = 1.4;
    treeCtx.beginPath();
    treeCtx.arc(x, y, radius + Math.max(1.4, Math.min(3, cellH * 0.3)), 0, Math.PI * 2);
    treeCtx.stroke();
  }
  if (showLabels) {
    treeCtx.fillStyle = labelInk;
    treeCtx.font = LABEL_FONT;
    treeCtx.textBaseline = "middle";
    treeCtx.textAlign = "left";
    treeCtx.fillText(fitTipLabel(data.genomes[node.genome] || ""), x + radius + 4, y);
  }
}
function fitTipLabel(text) {
  const budget = Math.max(24, treeWidth - xForNode({depth: maxDepth, distance: maxDistance, children: []}) - 10);
  const value = String(text);
  if (treeCtx.measureText(value).width <= budget) return value;
  for (let len = value.length - 1; len > 1; len -= 1) {
    const candidate = `${value.slice(0, len)}…`;
    if (treeCtx.measureText(candidate).width <= budget) return candidate;
  }
  return "…";
}
/** Row span covered by a subtree, memoised: used to skip offscreen branches. */
function nodeRowSpan(node) {
  if (node._span) return node._span;
  if (!node.children || !node.children.length) {
    node._span = {first: node.row, last: node.row};
    return node._span;
  }
  let first = Infinity;
  let last = -Infinity;
  for (const child of node.children) {
    const span = nodeRowSpan(child);
    if (span.first < first) first = span.first;
    if (span.last > last) last = span.last;
  }
  node._span = {first, last};
  return node._span;
}
function treeEdgeColor(child) {
  if (colorMode.value !== "lineage") return gridStroke;
  const lineage = subtreeLineage(child);
  if (!lineage || lineage === "mixed") return gridStroke;
  return lineageColor(lineage);
}
function renderLineageStrip(viewH, top, firstRow, lastRow) {
  if (!lineageStripCanvas || !lineageStripCtx) return;
  sizeCanvas(
    lineageStripCanvas,
    lineageStripCtx,
    LINEAGE_STRIP_W,
    viewH,
    canvasScale(LINEAGE_STRIP_W, viewH),
  );
  lineageStripCtx.fillStyle = panelBg;
  lineageStripCtx.fillRect(0, 0, LINEAGE_STRIP_W, viewH);
  for (let r = firstRow; r < lastRow; r += 1) {
    const genomeIdx = leafOrder[r];
    const lineage = lineageOfGenomeName(data.genomes[genomeIdx] || "");
    lineageStripCtx.fillStyle = lineageColor(lineage);
    lineageStripCtx.fillRect(0, r * cellH - top, LINEAGE_STRIP_W, Math.max(1, cellH));
  }
}

/** One labelled count with its share, used by both inspector tabs. */
function countRow(label, value, total, key) {
  const share = total ? `${(value / total * 100).toFixed(1)}%` : "0%";
  return `<div class="count-row"><span class="sw ${key}"></span>` +
    `<span class="count-name">${escapeHtml(label)}</span>` +
    `<b class="num">${value.toLocaleString()}</b><small class="num">${share}</small></div>`;
}

/** Prevalence as a bar, because "16 / 16" and "3 / 16" should not look alike. */
function prevalenceBar(present, total, pct) {
  const width = total ? Math.max(1, Math.round(present / total * 100)) : 0;
  return `<div class="prevalence-bar"><i style="width:${width}%"></i></div>` +
    `<p class="prevalence num">${present.toLocaleString()} / ${total.toLocaleString()} genomes &middot; ${pct}</p>`;
}

/* The inspector reflects the selected cluster; keyed so it is not rebuilt on
   every scroll frame. */
function renderClusterDetails() {
  const cluster = data.clusters[selected];
  const nextDetailsKey = cluster
    ? `${selected}|${referenceOnly}|${mechanismFilterKey()}|${visibleClusters.length}`
    : "none";
  if (nextDetailsKey !== detailsKey && cluster) {
    const counts = countsForCluster(selected);
    const mechanism = processKind(cluster);
    const mechanismLabel = processLabel(mechanism);
    const mechanismDescription = processDescription(mechanism);
    const mechanismColor = processColor(mechanism);
    const prevalencePct = `${Math.round(cluster.prevalence * 1000) / 10}%`;
    const refRow = referenceMechanismRow(cluster);
    const supportRow = cluster.support
      ? `<dt>Support tier</dt><dd>${escapeHtml(cluster.support)}</dd>`
      : "";
    const phyloRow = cluster.phylo
      ? `<dt>Phylogenetic pattern</dt><dd>${escapeHtml(cluster.phylo)}</dd>`
      : "";
    const detailRow = cluster.mechanismDetail
      ? `<dt>Mechanism detail</dt><dd>${escapeHtml(cluster.mechanismDetail)}</dd>`
      : "";
    const wordingRow = cluster.wording
      ? `<dt>Recommended wording</dt><dd><span class="process-wording">${escapeHtml(cluster.wording)}</span></dd>`
      : "";
    const total = Math.max(1, data.genomes.length);
    details.innerHTML = `
      <div class="cluster-head">
        <h3>${escapeHtml(cluster.name)}</h3>
        <p class="annotation">${escapeHtml(cluster.annotation || "No annotation provided")}</p>
      </div>
      ${prevalenceBar(cluster.present, total, prevalencePct)}
      <div class="card-title">Evidence in ${total.toLocaleString()} genomes</div>
      <div class="count-rows">
        ${countRow("present in input", counts.P, total, "p")}
        ${countRow("rescued by Recore", counts.R, total, "r")}
        ${countRow("supported by NSC", counts.N, total, "n")}
        ${countRow("absent", counts.A, total, "a")}
      </div>
      <div class="card-title">Mechanism</div>
      <span class="process-chip" style="--process-color:${mechanismColor}"><i></i>${escapeHtml(mechanismLabel)}</span>
      <p class="process-wording">${escapeHtml(mechanismDescription)}</p>
      <dl class="cluster-dl">
        ${supportRow}
        ${phyloRow}
        ${detailRow}
        ${wordingRow}
      </dl>
      ${refRow}
    `;
    renderProcessLegend();
    detailsKey = nextDetailsKey;
  } else if (nextDetailsKey !== detailsKey) {
    details.innerHTML = "";
    processLegend.innerHTML = "";
    detailsKey = nextDetailsKey;
  }
}

function render(scrollToSelected = false) {
  renderMatrix(scrollToSelected);
  renderClusterDetails();
  renderColumnHeader();
  renderSyntenyPanel();
  renderTree();
  renderLineageStrip();
  renderLineageLegend();
  renderTreeScaleBar();
  renderOverview();
  updateTopBarMeta();
}
function scheduleRender(scrollToSelected = false) {
  if (scrollToSelected) {
    render(true);
    return;
  }
  if (renderQueued) return;
  renderQueued = true;
  requestAnimationFrame(() => {
    renderQueued = false;
    render(false);
  });
}
function setHover(row, col, text) {
  const changed = hoverRow !== row || hoverCol !== col;
  hoverRow = row;
  hoverCol = col;
  if (text) status.textContent = text;
  if (changed) scheduleRender(false);
}
function setTreeMode(mode) {
  treeMode = mode;
  phylogramButton.classList.toggle("active", mode === "phylogram");
  cladogramButton.classList.toggle("active", mode === "cladogram");
  treeDirty = true;
  render(false);
}
function resetMatrixCache() {
  hoverCol = -1;
  treeDirty = true;
}
refreshThemeColors();
const defaultRootMode = hasL8Tip(data.tree) ? "l8" : "input";
applyRoot(defaultRootMode);
if (rootMode) rootMode.value = defaultRootMode;
computeGenomeStats();
refreshVisibleClusters();
// Open at a zoom that suits the run. A 16-genome, 12-cluster matrix at the old
// fixed sliders was a thumbnail in a half-empty pane; a 3000-tip tree needed the
// opposite. Column width is fitted to the space the matrix pane actually has.
function fitZoomToData() {
  const tips = Math.max(1, leafOrder.length);
  const cols = Math.max(1, visibleClusters.length);
  const bodyW = viewerBody ? viewerBody.clientWidth : 0;
  const matrixW = bodyW > 0
    ? Math.max(220, bodyW - layoutState.treeW - layoutState.inspectorW - 24)
    : FIT_TARGET_W;
  const viewH = matrixViewport ? matrixViewport.clientHeight : 0;
  const targetH = viewH > 120 ? viewH - 4 : FIT_TARGET_H;
  const rowFit = Math.max(MIN_ROW_ZOOM, Math.min(MAX_ROW_ZOOM, Math.floor(targetH / tips)));
  // cellW is round(slider * 0.82), so undo that factor when fitting to a width.
  const colFit = Math.max(MIN_COL_ZOOM, Math.min(MAX_COL_ZOOM, Math.floor((matrixW - 2) / cols / 0.82)));
  rowHeightControl.value = String(rowFit);
  sizeControl.value = String(colFit);
  labelToggle.checked = tips <= 500 && rowFit >= 4;
}
populateOptions();
if (layoutState.syntenyHiddenStored === false) {
  syntenyHidden = false;
}
applyLayout();
applySyntenyVisibility();
fitZoomToData();
setupDragHandle(treeHandle, "tree");
setupDragHandle(document.getElementById("inspectorHandle"), "inspector");
setupDragHandle(syntenyHandle, "synteny");
if (typeof ResizeObserver !== "undefined") {
  let lastWidth = 0;
  const ro = new ResizeObserver(entries => {
    const width = Math.round(entries[0]?.contentRect.width || 0);
    if (width === lastWidth) return;
    lastWidth = width;
    applyLayout();
    resetMatrixCache();
    render(false);
  });
  ro.observe(viewerBody);
}
let lastPixelRatio = window.devicePixelRatio || 1;
window.addEventListener("resize", () => {
  const ratio = window.devicePixelRatio || 1;
  if (ratio !== lastPixelRatio) {
    // Moved to a screen with a different density, or the page was zoomed.
    lastPixelRatio = ratio;
  }
  applyLayout();
  resetMatrixCache();
  syntenyKey = "";
  render(false);
});
status.textContent = defaultStatus;
render(true);
// The first fit runs before the column header knows how tall it needs to be, and
// the header can halve once labels turn horizontal. One more pass, after the
// panes have their real size, so the data actually fills them.
requestAnimationFrame(() => {
  fitZoomToData();
  resetMatrixCache();
  render(true);
});
window.addEventListener("recore:theme-change", () => {
  refreshThemeColors();
  resetMatrixCache();
  render(false);
});
/* ---- Pane menus ----------------------------------------------------------
   One menu per pane, holding everything that shapes that pane. Only one is open
   at a time, and a click anywhere else closes it. */
const paneMenus = [
  {menu: treeMenu, button: treeMenuBtn},
  {menu: matrixMenu, button: matrixMenuBtn},
].filter(entry => entry.menu && entry.button);

function setPaneMenu(target, open) {
  paneMenus.forEach(({menu, button}) => {
    const show = menu === target && open;
    menu.hidden = !show;
    button.setAttribute("aria-expanded", String(show));
    button.classList.toggle("active", show);
  });
  syncMenuIndicators();
}

function anyPaneMenuOpen() {
  return paneMenus.some(({menu}) => !menu.hidden);
}

/* The buttons have to say something when their menu is shut, or a filter left on
   two minutes ago silently explains a matrix that looks wrong. */
function syncMenuIndicators() {
  if (matrixMenuBtn && matrixMenu.hidden) {
    matrixMenuBtn.classList.toggle("active", mechanismFilters.size > 0 || referenceOnly);
  }
  if (treeMenuBtn && treeMenu.hidden) {
    treeMenuBtn.classList.toggle("active", rootMode ? rootMode.value !== defaultRootMode : false);
  }
}

paneMenus.forEach(({menu, button}) => {
  button.addEventListener("click", event => {
    event.stopPropagation();
    setPaneMenu(menu, menu.hidden);
  });
  menu.addEventListener("click", event => event.stopPropagation());
});
document.addEventListener("pointerdown", event => {
  if (!anyPaneMenuOpen()) return;
  if (event.target.closest(".pane-menu, #treeMenuBtn, #matrixMenuBtn")) return;
  setPaneMenu(null, false);
});

colorMode.addEventListener("change", () => { treeDirty = true; render(false); });
labelToggle.addEventListener("change", () => { treeDirty = true; render(false); });
phylogramButton.addEventListener("click", () => setTreeMode("phylogram"));
cladogramButton.addEventListener("click", () => setTreeMode("cladogram"));
referenceToggle.addEventListener("change", () => {
  if (!data.referenceGenome) {
    referenceToggle.checked = false;
    return;
  }
  referenceOnly = referenceToggle.checked;
  refreshVisibleClusters();
  populateOptions();
  resetMatrixCache();
  render(true);
});
if (prevClusterButton) prevClusterButton.addEventListener("click", () => stepCluster(-1));
if (nextClusterButton) nextClusterButton.addEventListener("click", () => stepCluster(1));
function setSynteny(open) {
  syntenyHidden = !open;
  applySyntenyVisibility();
  render(false);
}
syntenyToggle.addEventListener("click", () => setSynteny(syntenyHidden));
if (closeSyntenyButton) closeSyntenyButton.addEventListener("click", () => setSynteny(false));
syntenySvg.addEventListener("click", event => {
  const gene = event.target.closest(".synteny-gene[data-gene-cluster]");
  if (gene) {
    const clusterName = gene.getAttribute("data-gene-cluster") || "";
    if (selectMatrixClusterByName(clusterName)) {
      event.preventDefault();
      event.stopPropagation();
    }
    return;
  }
  const referenceRow = event.target.closest(".synteny-reference-row[data-reference-context]");
  if (referenceRow) {
    selectedSyntenyReference = true;
    syntenyKey = "";
    treeDirty = true;
    render(false);
    return;
  }
  const row = event.target.closest(".synteny-row[data-block-index]");
  if (!row) return;
  selectedSyntenyBlock = Number(row.getAttribute("data-block-index") || 0);
  selectedSyntenyReference = false;
  syntenyKey = "";
  treeDirty = true;
  render(false);
});
syntenySvg.addEventListener("mouseover", event => {
  const gene = event.target.closest(".synteny-gene[data-gene-cluster]");
  if (gene) {
    highlightSyntenyCluster(gene.getAttribute("data-gene-cluster"));
    return;
  }
  const link = event.target.closest(".synteny-link[data-link-cluster]");
  if (link) {
    highlightSyntenyCluster(link.getAttribute("data-link-cluster"));
  }
});
syntenySvg.addEventListener("mouseout", event => {
  const target = event.target;
  if (!target) return;
  if (target.closest(".synteny-gene[data-gene-cluster]") || target.closest(".synteny-link[data-link-cluster]")) {
    clearSyntenyHighlight();
  }
});
function highlightSyntenyCluster(clusterName) {
  if (!clusterName) return;
  syntenySvg.querySelectorAll(".synteny-gene[data-gene-cluster]").forEach(node => {
    const match = node.getAttribute("data-gene-cluster") === clusterName;
    node.classList.toggle("cluster-highlight", match);
    node.classList.toggle("dimmed", !match && !node.classList.contains("synteny-target"));
  });
  syntenySvg.querySelectorAll(".synteny-link[data-link-cluster]").forEach(link => {
    const match = link.getAttribute("data-link-cluster") === clusterName;
    link.classList.toggle("highlight", match);
    link.classList.toggle("dimmed", !match);
  });
}
function clearSyntenyHighlight() {
  syntenySvg.querySelectorAll(".synteny-gene").forEach(node => {
    node.classList.remove("cluster-highlight", "dimmed");
  });
  syntenySvg.querySelectorAll(".synteny-link").forEach(link => {
    link.classList.remove("highlight", "dimmed");
  });
}
processLegend.addEventListener("click", event => {
  const clear = event.target.closest("[data-clear-mechanism]");
  const card = event.target.closest("[data-mechanism]");
  if (!clear && !card) return;
  if (clear) {
    mechanismFilters.clear();
  } else {
    const mechanism = card.getAttribute("data-mechanism") || "";
    if (mechanismFilters.has(mechanism)) {
      mechanismFilters.delete(mechanism);
    } else {
      mechanismFilters.add(mechanism);
    }
  }
  refreshVisibleClusters();
  populateOptions();
  resetMatrixCache();
  syncMenuIndicators();
  render(true);
});
async function toggleFullscreen() {
  if (!document.fullscreenElement && explorerSection.requestFullscreen) {
    await explorerSection.requestFullscreen();
  } else if (document.exitFullscreen) {
    await document.exitFullscreen();
  }
}
fullscreenButton.addEventListener("click", toggleFullscreen);
document.addEventListener("fullscreenchange", () => {
  // Setting textContent here used to wipe the inline SVG and leave a blank box.
  const on = document.fullscreenElement === explorerSection;
  fullscreenButton.classList.toggle("active", on);
  fullscreenButton.setAttribute("aria-pressed", String(on));
  fullscreenButton.setAttribute("data-tip", on ? "Leave fullscreen (F)" : "Fullscreen (F)");
});
sizeControl.addEventListener("input", () => { resetMatrixCache(); render(false); });
rowHeightControl.addEventListener("input", () => { resetMatrixCache(); render(false); });
syntenyZoomControl.addEventListener("input", () => { syntenyKey = ""; render(false); });
matrixColorMode.addEventListener("change", () => { resetMatrixCache(); render(false); });
if (columnOrder) {
  columnOrder.addEventListener("change", () => {
    refreshVisibleClusters();
    populateOptions();
    resetMatrixCache();
    render(true);
  });
}
matrixPane.addEventListener("scroll", () => { headerScroller.scrollLeft = matrixPane.scrollLeft; });
headerScroller.addEventListener("scroll", () => { matrixPane.scrollLeft = headerScroller.scrollLeft; });
function adjustRange(control, delta, min, max) {
  const value = Number(control.value || min);
  const next = Math.max(min, Math.min(max, value + delta));
  if (next === value) return false;
  control.value = String(next);
  return true;
}
function wheelZoom(event, forcedAxis = "") {
  const columnZoom = forcedAxis === "column" || event.shiftKey;
  const rowZoom = forcedAxis === "row" || event.ctrlKey || event.metaKey || event.altKey;
  if (!columnZoom && !rowZoom) return;
  event.preventDefault();
  event.stopPropagation();
  const delta = event.deltaY < 0 ? 1 : -1;
  const changed = columnZoom
    ? adjustRange(sizeControl, delta, MIN_COL_ZOOM, MAX_COL_ZOOM)
    : adjustRange(rowHeightControl, delta, MIN_ROW_ZOOM, MAX_ROW_ZOOM);
  if (changed) {
    resetMatrixCache();
    render(false);
  }
}
viewerBody.addEventListener("wheel", wheelZoom, {passive:false});
if (treePane) treePane.addEventListener("wheel", event => wheelZoom(event, "row"), {passive:false});
if (treeScroll) treeScroll.addEventListener("wheel", event => wheelZoom(event, "row"), {passive:false});
headerScroller.addEventListener("wheel", event => wheelZoom(event, "column"), {passive:false});
matrixPane.addEventListener("wheel", wheelZoom, {passive:false});

/* The panes share one scroll origin: the matrix scroller is the source of truth,
   the header follows it horizontally and the tree vertically. Before this the
   tree and the matrix scrolled independently and drifted out of alignment. */
let syncingScroll = false;
function syncFromMatrix() {
  if (syncingScroll) return;
  syncingScroll = true;
  if (headerScroller) headerScroller.scrollLeft = matrixPane.scrollLeft;
  if (treeScroll) treeScroll.scrollTop = matrixPane.scrollTop;
  syncingScroll = false;
  scheduleRender(false);
}
matrixPane.addEventListener("scroll", syncFromMatrix, {passive: true});
if (headerScroller) {
  headerScroller.addEventListener("scroll", () => {
    if (syncingScroll) return;
    syncingScroll = true;
    matrixPane.scrollLeft = headerScroller.scrollLeft;
    syncingScroll = false;
    scheduleRender(false);
  }, {passive: true});
}
if (treeScroll) {
  treeScroll.addEventListener("scroll", () => {
    if (syncingScroll) return;
    syncingScroll = true;
    matrixPane.scrollTop = treeScroll.scrollTop;
    syncingScroll = false;
    scheduleRender(false);
  }, {passive: true});
}

/** Row under the pointer, in matrix coordinates. */
function rowAtPointer(event, element) {
  const rect = element.getBoundingClientRect();
  const r = Math.floor((event.clientY - rect.top + matrixScrollTop()) / cellH);
  return r >= 0 && r < leafOrder.length ? r : -1;
}
/** Visible column under the pointer, in matrix coordinates. */
function columnAtPointer(event, element) {
  const rect = element.getBoundingClientRect();
  const c = Math.floor((event.clientX - rect.left + matrixScrollLeft()) / cellW);
  return c >= 0 && c < visibleClusters.length ? c : -1;
}
function selectVisibleColumn(c) {
  if (c < 0) return;
  selected = originalClusterAtVisibleColumn(c);
  selectedSyntenyReference = false;
  populateOptions();
  treeDirty = true;
  render(true);
}

if (treeScroll) {
  treeScroll.addEventListener("mousemove", event => {
    const r = rowAtPointer(event, treeScroll);
    if (r < 0) return;
    const genomeIdx = leafOrder[r];
    const state = stateAt(selected, genomeIdx);
    const cluster = data.clusters[selected];
    const tip = tipStyle(genomeIdx);
    const tipNode = findLeafByGenome(treeRoot, genomeIdx);
    const branchInfo = tipNode
      ? `; root distance ${formatBranchLength(tipNode.distance)}, terminal branch ${formatBranchLength(tipNode.length)}`
      : "";
    const detail = colorMode.value === "cluster" && cluster
      ? `${cluster.name} / ${data.genomes[genomeIdx]}: ${labels[state]}${branchInfo}`
      : `${data.genomes[genomeIdx]}: ${tip.title}${branchInfo}`;
    tooltip.style.display = "none";
    setHover(r, -1, detail);
  });
  treeScroll.addEventListener("mouseleave", () => setHover(-1, -1, defaultStatus));
}

headerScroller.addEventListener("mousemove", event => {
  const c = columnAtPointer(event, headerScroller);
  if (c < 0) return;
  const cluster = data.clusters[originalClusterAtVisibleColumn(c)];
  if (!cluster) return;
  setHover(-1, c, `${cluster.name}: ${cluster.annotation || "No annotation provided"}`);
});
headerScroller.addEventListener("mouseleave", () => setHover(-1, -1, defaultStatus));
headerScroller.addEventListener("click", event => selectVisibleColumn(columnAtPointer(event, headerScroller)));

matrixPane.addEventListener("mousemove", event => {
  const c = columnAtPointer(event, matrixPane);
  const r = rowAtPointer(event, matrixPane);
  if (c < 0 || r < 0) {
    tooltip.style.display = "none";
    return;
  }
  const clusterIdx = originalClusterAtVisibleColumn(c);
  const genomeIdx = leafOrder[r];
  const state = stateAt(clusterIdx, genomeIdx);
  const cluster = data.clusters[clusterIdx];
  const mechanism = processKind(cluster);
  setHover(r, c, `${cluster.name} / ${data.genomes[genomeIdx]}: ${labels[state]}; ${processLabel(mechanism)}`);
  const paneRect = matrixPane.getBoundingClientRect();
  tooltip.style.display = "block";
  tooltip.style.left = `${Math.min(event.clientX - paneRect.left + 14, paneRect.width - 300)}px`;
  tooltip.style.top = `${event.clientY - paneRect.top + 14}px`;
  tooltip.innerHTML = `<strong>${escapeHtml(cluster.name)}</strong><br>${escapeHtml(data.genomes[genomeIdx])}: ${escapeHtml(labels[state])}<br><span>${escapeHtml(processLabel(mechanism))}</span><br><span>${escapeHtml(cluster.annotation || "")}</span>`;
});
matrixPane.addEventListener("mouseleave", () => {
  tooltip.style.display = "none";
  setHover(-1, -1, defaultStatus);
});
matrixPane.addEventListener("click", event => selectVisibleColumn(columnAtPointer(event, matrixPane)));

/* ---- Print and export ----------------------------------------------------
   Canvases are bitmaps sized to the on-screen viewport, so printing without
   redrawing them clips the panes to whatever the screen happened to be. Redraw
   once the print layout is in place, and again when it is torn down. */
window.addEventListener("beforeprint", () => {
  applyLayout();
  resetMatrixCache();
  syntenyKey = "";
  render(false);
});
window.addEventListener("afterprint", () => {
  applyLayout();
  resetMatrixCache();
  syntenyKey = "";
  render(false);
});

function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
/* Compose the panes as drawn into one bitmap at 2x. The report's job is to end
   up in a manuscript, and screenshotting a browser window is not a figure
   pipeline. */
function exportViewPng() {
  const scale = 2;
  const gap = 8;
  const treeBodyW = treeCanvas ? treeCanvas.clientWidth : 0;
  const treeW = treeBodyW + LINEAGE_STRIP_W;
  const headTop = headerCanvas ? headerCanvas.clientHeight : 0;
  const matrixW = canvas ? canvas.clientWidth : 0;
  const bodyH = canvas ? canvas.clientHeight : 0;
  const width = treeW + gap + matrixW;
  const height = headTop + bodyH;
  if (width < 2 || height < 2) return;
  const out = document.createElement("canvas");
  out.width = Math.round(width * scale);
  out.height = Math.round(height * scale);
  const octx = out.getContext("2d");
  octx.setTransform(scale, 0, 0, scale, 0, 0);
  octx.fillStyle = panelBg;
  octx.fillRect(0, 0, width, height);
  // A canvas with a zero dimension throws from drawImage, and any pane can be
  // collapsed by the user or by a narrow window.
  const paste = (source, x, y, w, h) => {
    if (!source || !source.width || !source.height || w < 1 || h < 1) return;
    octx.drawImage(source, x, y, w, h);
  };
  paste(headerCanvas, treeW + gap, 0, matrixW, headTop);
  paste(treeCanvas, 0, headTop, treeBodyW, bodyH);
  paste(lineageStripCanvas, treeBodyW, headTop, LINEAGE_STRIP_W, bodyH);
  paste(canvas, treeW + gap, headTop, matrixW, bodyH);
  out.toBlob(blob => {
    if (blob) downloadBlob(blob, "recore_matrix.png");
  }, "image/png");
}
/* Synteny is already vector, so export it as vector: gene arrows have to stay
   crisp at print resolution. */
function exportSyntenySvg() {
  if (!syntenySvg || syntenyHidden) return false;
  const clone = syntenySvg.cloneNode(true);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  // Inline the rules the drawing relies on, or it opens unstyled outside the page.
  const rules = [];
  for (const sheet of document.styleSheets) {
    let cssRules = null;
    try {
      cssRules = sheet.cssRules;
    } catch (_e) {
      continue;
    }
    for (const rule of cssRules || []) {
      const selector = rule.selectorText || "";
      if (selector.includes("synteny") || selector.startsWith(":root")) {
        rules.push(rule.cssText);
      }
    }
  }
  const style = document.createElementNS("http://www.w3.org/2000/svg", "style");
  style.textContent = rules.join("\n");
  clone.insertBefore(style, clone.firstChild);
  const blob = new Blob([new XMLSerializer().serializeToString(clone)], {
    type: "image/svg+xml",
  });
  downloadBlob(blob, "recore_synteny.svg");
  return true;
}
if (exportButton) {
  exportButton.addEventListener("click", () => {
    exportViewPng();
    const svg = exportSyntenySvg();
    status.textContent = svg
      ? "Exported recore_matrix.png and recore_synteny.svg."
      : "Exported recore_matrix.png.";
  });
}


if (rootMode) {
  rootMode.addEventListener("change", () => {
    applyRoot(rootMode.value);
    render(false);
    syncMenuIndicators();
    status.textContent = `Tree rerooted: ${rootMode.options[rootMode.selectedIndex].textContent}.`;
  });
}

/** Root the displayed tree on one genome, and remember it in the root menu. */
function rootOnGenome(genomeIdx) {
  const name = data.genomes[genomeIdx] || "";
  if (rootMode) {
    const option = [...rootMode.options].find(o => o.value.startsWith("genome:"))
      || rootMode.appendChild(document.createElement("option"));
    option.value = `genome:${genomeIdx}`;
    option.textContent = name;
    rootMode.value = option.value;
  }
  applyRoot(`genome:${genomeIdx}`);
  render(false);
  syncMenuIndicators();
  status.textContent = `Tree rooted on ${name}. The tree menu takes you back.`;
}

/* ---- Genome inspection ----------------------------------------------------
   A plain tip click used to reroot the whole tree, which is a big move to
   trigger by accident and left no way to simply ask about a genome. Click now
   inspects; Alt-click reroots. */
let selectedGenome = -1;

function setInspectorTab(name) {
  const wantGenome = name === "genome";
  if (clusterTab) {
    clusterTab.classList.toggle("active", !wantGenome);
    clusterTab.setAttribute("aria-selected", String(!wantGenome));
  }
  if (genomeTab) {
    genomeTab.classList.toggle("active", wantGenome);
    genomeTab.setAttribute("aria-selected", String(wantGenome));
  }
  if (details) details.hidden = wantGenome;
  if (genomeDetails) genomeDetails.hidden = !wantGenome;
}

if (clusterTab) clusterTab.addEventListener("click", () => setInspectorTab("cluster"));
if (genomeTab) genomeTab.addEventListener("click", () => setInspectorTab("genome"));

function inspectGenome(genomeIdx, {focus = true} = {}) {
  if (genomeIdx < 0) return;
  selectedGenome = genomeIdx;
  renderGenomeDetails();
  if (focus) setInspectorTab("genome");
}

function renderGenomeDetails() {
  if (!genomeDetails) return;
  if (selectedGenome < 0) {
    genomeDetails.innerHTML = "<p class='muted'>Click a tip in the phylogeny to inspect that genome.</p>";
    return;
  }
  const name = data.genomes[selectedGenome] || "";
  const totals = genomeStats[selectedGenome] || {P: 0, R: 0, N: 0, A: 0};
  const assessed = totals.P + totals.R + totals.N + totals.A;
  const lineage = lineageOfGenomeName(name);
  const node = findLeafByGenome(treeRoot, selectedGenome);
  const state = stateAt(selected, selectedGenome);
  const cluster = data.clusters[selected];
  const row = rowOfGenome[selectedGenome];
  const isReference = data.referenceGenome && data.referenceGenome.index === selectedGenome;
  genomeDetails.innerHTML = `
    <div class="cluster-head">
      <h3>${escapeHtml(name)}</h3>
      ${isReference ? "<span class='pill'>reference</span>" : ""}
    </div>
    <p class="muted">Tip ${row >= 0 ? row + 1 : "?"} of ${leafOrder.length}${lineage && lineage !== "unknown" ? ` &middot; lineage ${escapeHtml(lineage)}` : ""}</p>
    <dl class="cluster-dl">
      <dt>In the selected cluster</dt>
      <dd>${escapeHtml(cluster ? cluster.name : "none")}: <strong>${escapeHtml(labels[state])}</strong></dd>
    </dl>
    <div class="card-title">Cluster content across ${assessed.toLocaleString()} clusters</div>
    <div class="count-rows">
      ${countRow("present in input", totals.P, assessed, "p")}
      ${countRow("rescued by Recore", totals.R, assessed, "r")}
      ${countRow("supported by NSC", totals.N, assessed, "n")}
      ${countRow("absent", totals.A, assessed, "a")}
    </div>
    <dl class="cluster-dl">
      <dt>Position on the tree</dt>
      <dd>${node
        ? `root distance ${escapeHtml(formatBranchLength(node.distance))}, terminal branch ${escapeHtml(formatBranchLength(node.length))}`
        : "not on the displayed tree"}</dd>
    </dl>
    <div class="genome-actions">
      <button type="button" class="ghost-btn" data-root-here="${selectedGenome}">Root the tree here</button>
      <button type="button" class="ghost-btn" data-find-tip="${selectedGenome}">Scroll to this tip</button>
    </div>`;
}

if (genomeDetails) {
  genomeDetails.addEventListener("click", event => {
    const rootHere = event.target.closest("[data-root-here]");
    if (rootHere) {
      rootOnGenome(Number(rootHere.getAttribute("data-root-here")));
      return;
    }
    const findTip = event.target.closest("[data-find-tip]");
    if (findTip) scrollToGenome(Number(findTip.getAttribute("data-find-tip")));
  });
}

/** Centre the panes on one genome's row and ring it on the tree. */
function scrollToGenome(genomeIdx) {
  const row = rowOfGenome[genomeIdx];
  if (row < 0 || !matrixPane) return;
  const view = viewportSize(matrixViewport);
  matrixPane.scrollTop = Math.max(0, row * cellH - view.h / 2);
  if (treeScroll) treeScroll.scrollTop = matrixPane.scrollTop;
  highlightedRows = new Set([row]);
  treeDirty = true;
  render(false);
  status.textContent = `${data.genomes[genomeIdx]}: tip ${row + 1} of ${leafOrder.length}`;
}

if (treeScroll) {
  treeScroll.addEventListener("click", event => {
    const row = rowAtPointer(event, treeScroll);
    if (row < 0) return;
    const genomeIdx = leafOrder[row];
    if (event.altKey) {
      rootOnGenome(genomeIdx);
      return;
    }
    inspectGenome(genomeIdx);
  });
}

/* ---- Find palette ---------------------------------------------------------
   Clusters and genomes answer to the same shortcut, because "where is X" is one
   question and the reader should not have to know which kind of thing X is. */
const PALETTE_LIMIT = 200;
let paletteScope = "all";
let paletteRows = [];
let paletteActive = 0;

function paletteIsOpen() {
  return palette && !palette.hidden;
}

function openPalette(scope) {
  if (!palette) return;
  if (scope) setPaletteScope(scope);
  palette.hidden = false;
  if (omnibox) omnibox.setAttribute("aria-expanded", "true");
  search.select();
  search.focus();
  refreshPalette();
}

function closePalette() {
  if (!palette) return;
  palette.hidden = true;
  if (omnibox) {
    omnibox.setAttribute("aria-expanded", "false");
    omnibox.focus();
  }
}

function setPaletteScope(scope) {
  paletteScope = scope;
  palette.querySelectorAll(".palette-tabs button").forEach(button => {
    const on = button.getAttribute("data-scope") === scope;
    button.classList.toggle("active", on);
    button.setAttribute("aria-selected", String(on));
  });
}

/* Name hits beat annotation hits, and a prefix beats a hit in the middle, so the
   thing you half-remembered the start of comes first. */
function scoreHit(name, note, query) {
  const lowName = name.toLowerCase();
  if (!query) return 3;
  if (lowName.startsWith(query)) return 0;
  if (lowName.includes(query)) return 1;
  if (note && note.toLowerCase().includes(query)) return 2;
  return -1;
}

function buildPaletteRows() {
  const query = search.value.trim().toLowerCase();
  const hits = [];
  if (paletteScope !== "genome") {
    for (let i = 0; i < visibleClusters.length; i += 1) {
      const idx = visibleClusters[i];
      const cluster = data.clusters[idx];
      const note = cluster.annotation || "";
      const score = scoreHit(cluster.name, note, query);
      if (score < 0) continue;
      hits.push({
        kind: "cluster",
        score,
        idx,
        name: cluster.name,
        note: note || "No annotation provided",
        badge: `${cluster.present || 0} / ${data.genomes.length}`,
      });
    }
  }
  if (paletteScope !== "cluster") {
    for (let row = 0; row < leafOrder.length; row += 1) {
      const genomeIdx = leafOrder[row];
      const name = data.genomes[genomeIdx] || "";
      const lineage = lineageOfGenomeName(name);
      const note = lineage && lineage !== "unknown" ? `lineage ${lineage}` : "";
      const score = scoreHit(name, note, query);
      if (score < 0) continue;
      hits.push({
        kind: "genome",
        score,
        idx: genomeIdx,
        name,
        note: note || "no lineage assigned",
        badge: `tip ${row + 1}`,
      });
    }
  }
  hits.sort((a, b) => a.score - b.score || a.name.localeCompare(b.name));
  return {total: hits.length, rows: hits.slice(0, PALETTE_LIMIT)};
}

function markMatch(text, query) {
  if (!query) return escapeHtml(text);
  const at = text.toLowerCase().indexOf(query);
  if (at < 0) return escapeHtml(text);
  return `${escapeHtml(text.slice(0, at))}<mark>${escapeHtml(text.slice(at, at + query.length))}</mark>${escapeHtml(text.slice(at + query.length))}`;
}

function refreshPalette() {
  if (!paletteList) return;
  const query = search.value.trim().toLowerCase();
  const {total, rows} = buildPaletteRows();
  paletteRows = rows;
  paletteActive = 0;
  if (!rows.length) {
    paletteList.innerHTML = `<li class="palette-empty">Nothing matches ${query ? `"${escapeHtml(query)}"` : "the current filters"}.</li>`;
  } else {
    paletteList.innerHTML = rows.map((row, i) =>
      `<li class="palette-item${i === 0 ? " active" : ""}" role="option" aria-selected="${i === 0}" data-index="${i}">` +
      `<span class="palette-kind">${row.kind}</span>` +
      `<span class="palette-name">${markMatch(row.name, query)}</span>` +
      `<span class="palette-note">${escapeHtml(row.note)}</span>` +
      `<span class="palette-badge">${escapeHtml(row.badge)}</span></li>`
    ).join("");
  }
  if (paletteCount) {
    paletteCount.textContent = total > rows.length
      ? `${rows.length} of ${total.toLocaleString()}`
      : `${total.toLocaleString()} result${total === 1 ? "" : "s"}`;
  }
}

function movePaletteActive(delta) {
  if (!paletteRows.length) return;
  paletteActive = (paletteActive + delta + paletteRows.length) % paletteRows.length;
  paletteList.querySelectorAll(".palette-item").forEach((item, i) => {
    const on = i === paletteActive;
    item.classList.toggle("active", on);
    item.setAttribute("aria-selected", String(on));
    if (on) item.scrollIntoView({block: "nearest"});
  });
}

function commitPaletteRow(row) {
  if (!row) return;
  if (row.kind === "cluster") {
    selectMatrixClusterByName(row.name);
    setInspectorTab("cluster");
  } else {
    inspectGenome(row.idx);
    scrollToGenome(row.idx);
  }
  closePalette();
}

if (palette) {
  if (omnibox) omnibox.addEventListener("click", () => openPalette());
  if (paletteScrim) paletteScrim.addEventListener("click", closePalette);
  search.addEventListener("input", refreshPalette);
  palette.querySelectorAll(".palette-tabs button").forEach(button => {
    button.addEventListener("click", () => {
      setPaletteScope(button.getAttribute("data-scope") || "all");
      refreshPalette();
      search.focus();
    });
  });
  paletteList.addEventListener("click", event => {
    const item = event.target.closest(".palette-item");
    if (item) commitPaletteRow(paletteRows[Number(item.getAttribute("data-index"))]);
  });
  paletteList.addEventListener("mousemove", event => {
    const item = event.target.closest(".palette-item");
    if (!item) return;
    const index = Number(item.getAttribute("data-index"));
    if (index === paletteActive) return;
    movePaletteActive(index - paletteActive);
  });
  search.addEventListener("keydown", event => {
    if (event.key === "ArrowDown") { event.preventDefault(); movePaletteActive(1); }
    else if (event.key === "ArrowUp") { event.preventDefault(); movePaletteActive(-1); }
    else if (event.key === "Enter") { event.preventDefault(); commitPaletteRow(paletteRows[paletteActive]); }
    else if (event.key === "Escape") { event.preventDefault(); closePalette(); }
    else if (event.key === "Tab") {
      event.preventDefault();
      const order = ["all", "cluster", "genome"];
      setPaletteScope(order[(order.indexOf(paletteScope) + 1) % order.length]);
      refreshPalette();
    }
  });
}

/* ---- Shortcuts ------------------------------------------------------------ */
const isApple = /Mac|iPhone|iPad/.test(navigator.platform || navigator.userAgent || "");
if (omniboxKey) omniboxKey.textContent = isApple ? "⌘K" : "Ctrl K";
document.querySelectorAll(".status-hint kbd").forEach(node => {
  if (node.textContent === "⌘K") node.textContent = isApple ? "⌘K" : "Ctrl K";
});

function setSheetOpen(open) {
  if (!shortcutSheet) return;
  shortcutSheet.hidden = !open;
}
if (shortcutsButton) shortcutsButton.addEventListener("click", () => setSheetOpen(shortcutSheet.hidden));
const sheetScrim = document.getElementById("sheetScrim");
if (sheetScrim) sheetScrim.addEventListener("click", () => setSheetOpen(false));
document.getElementById("closeShortcuts")?.addEventListener("click", () => setSheetOpen(false));

function typingInField(target) {
  if (!target) return false;
  const tag = target.tagName;
  return tag === "INPUT" || tag === "SELECT" || tag === "TEXTAREA" || target.isContentEditable;
}

document.addEventListener("keydown", event => {
  if ((event.metaKey || event.ctrlKey) && event.key.toLowerCase() === "k") {
    event.preventDefault();
    if (paletteIsOpen()) closePalette();
    else openPalette();
    return;
  }
  if (event.key === "Escape") {
    if (paletteIsOpen()) { closePalette(); return; }
    if (shortcutSheet && !shortcutSheet.hidden) { setSheetOpen(false); return; }
    if (anyPaneMenuOpen()) { setPaneMenu(null, false); return; }
    return;
  }
  if (paletteIsOpen() || typingInField(event.target)) return;
  if (event.metaKey || event.ctrlKey || event.altKey) return;

  switch (event.key) {
    case "/":
      event.preventDefault();
      openPalette();
      return;
    case "?":
      event.preventDefault();
      setSheetOpen(shortcutSheet ? shortcutSheet.hidden : false);
      return;
    case "ArrowLeft":
      event.preventDefault();
      stepCluster(event.shiftKey ? -10 : -1);
      return;
    case "ArrowRight":
      event.preventDefault();
      stepCluster(event.shiftKey ? 10 : 1);
      return;
    default:
      break;
  }
  switch (event.key.toLowerCase()) {
    case "s":
      setSynteny(syntenyHidden);
      break;
    case "d":
      document.getElementById("openDetails")?.click();
      break;
    case "e":
      exportButton?.click();
      break;
    case "f":
      toggleFullscreen();
      break;
    case "t":
      document.getElementById("themeToggle")?.click();
      break;
    case "c":
      setTreeMode(treeMode === "phylogram" ? "cladogram" : "phylogram");
      break;
    case "r":
      if (rootMode) {
        rootMode.value = "midpoint";
        applyRoot("midpoint");
        render(false);
        syncMenuIndicators();
        status.textContent = "Tree rerooted: Midpoint.";
      }
      break;
    case "z":
      fitZoomToData();
      resetMatrixCache();
      render(true);
      status.textContent = "Zoom fitted to the data.";
      break;
    default:
      break;
  }
});

/* ---- Tree scale bar -------------------------------------------------------
   A phylogram without a scale is a picture rather than a measurement, and the
   branch lengths here are the whole argument for drawing one. */
function niceStep(raw) {
  if (!(raw > 0)) return 0;
  const magnitude = Math.pow(10, Math.floor(Math.log10(raw)));
  const scaled = raw / magnitude;
  const step = scaled >= 5 ? 5 : scaled >= 2 ? 2 : 1;
  return step * magnitude;
}

function formatDistance(value) {
  if (value >= 1) return value.toFixed(2).replace(/\.?0+$/, "");
  if (value >= 0.001) return value.toFixed(4).replace(/0+$/, "");
  return value.toExponential(1);
}

function renderTreeScaleBar() {
  if (!treeScaleBar) return;
  if (treeMode !== "phylogram" || !(maxDistance > 0)) {
    treeScaleBar.hidden = true;
    return;
  }
  const pxPerUnit = treeUsableWidth() / maxDistance;
  const step = niceStep(maxDistance / 4);
  const width = step * pxPerUnit;
  if (!(width >= 18) || !Number.isFinite(width)) {
    treeScaleBar.hidden = true;
    return;
  }
  treeScaleBar.hidden = false;
  treeScaleBar.innerHTML = `<i style="width:${Math.round(width)}px"></i><span>${escapeHtml(formatDistance(step))} substitutions/site</span>`;
}

/* ---- Whole-matrix overview ------------------------------------------------
   At 30,000 columns a scrollbar tells you how far along you are and nothing
   about what is there. Prevalence per column does, and clicking jumps. */
function renderOverview() {
  if (!overviewCanvas || !overviewStrip) return;
  // Below a few dozen columns the whole matrix is already on screen, and a
  // navigation aid for something you can see is just another band of colour.
  const useful = visibleClusters.length >= 40;
  overviewStrip.hidden = !useful;
  if (!useful) return;
  const width = overviewStrip.clientWidth;
  const height = overviewStrip.clientHeight;
  if (width < 4 || height < 4) return;
  const key = `${width}x${height}|${visibleClusters.length}|${selected}|${panelBg}`;
  if (key !== overviewKey) {
    overviewKey = key;
    const scale = canvasScale(width, height);
    const octx = overviewCanvas.getContext("2d");
    sizeCanvas(overviewCanvas, octx, width, height, scale);
    octx.clearRect(0, 0, width, height);
    const n = visibleClusters.length;
    const genomes = Math.max(1, data.genomes.length);
    if (n) {
      // One bar per device column: with more clusters than pixels, take the
      // strongest cluster in each bucket so a single high-prevalence column in
      // a crowd of empty ones is still visible.
      octx.fillStyle = colors.P;
      octx.globalAlpha = 0.45;
      for (let x = 0; x < width; x += 1) {
        const from = Math.floor(x / width * n);
        const to = Math.max(from + 1, Math.floor((x + 1) / width * n));
        let best = 0;
        for (let i = from; i < to && i < n; i += 1) {
          const present = data.clusters[visibleClusters[i]]?.present || 0;
          if (present > best) best = present;
        }
        const h = Math.max(1, Math.round(best / genomes * (height - 4)));
        octx.fillRect(x, height - 2 - h, 1, h);
      }
      octx.globalAlpha = 1;
      const at = visibleClusters.indexOf(selected);
      if (at >= 0) {
        octx.fillStyle = selectStroke;
        octx.fillRect(Math.floor(at / n * width), 0, 2, height);
      }
    }
  }
  if (overviewWindow) {
    const n = Math.max(1, visibleClusters.length);
    const view = viewportSize(matrixViewport);
    const first = matrixScrollLeft() / cellW;
    const span = Math.min(n, view.w / cellW);
    overviewWindow.style.left = `${(first / n * width).toFixed(1)}px`;
    overviewWindow.style.width = `${Math.max(2, span / n * width).toFixed(1)}px`;
  }
}

if (overviewStrip) {
  const jumpTo = clientX => {
    const rect = overviewStrip.getBoundingClientRect();
    const fraction = Math.max(0, Math.min(1, (clientX - rect.left) / Math.max(1, rect.width)));
    const view = viewportSize(matrixViewport);
    matrixPane.scrollLeft = Math.max(0, fraction * virtualW - view.w / 2);
  };
  overviewStrip.addEventListener("pointerdown", event => {
    overviewStrip.setPointerCapture(event.pointerId);
    jumpTo(event.clientX);
  });
  overviewStrip.addEventListener("pointermove", event => {
    if (event.buttons & 1) jumpTo(event.clientX);
  });
}

renderGenomeDetails();
populateOptions();
render(false);

})();
