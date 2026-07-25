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
const select = document.getElementById("clusterSelect");
const search = document.getElementById("clusterSearch");
const sizeControl = document.getElementById("cellSize");
const rowHeightControl = document.getElementById("rowHeight");
const colorMode = document.getElementById("tipColorMode");
const matrixColorMode = document.getElementById("matrixColorMode");
const labelToggle = document.getElementById("tipLabels");
const referenceButton = document.getElementById("referenceOnly");
const syntenyToggle = document.getElementById("toggleSynteny");
const syntenyZoomControl = document.getElementById("syntenyZoom");
const fullscreenButton = document.getElementById("explorerFullscreen");
const phylogramButton = document.getElementById("phylogramMode");
const cladogramButton = document.getElementById("cladogramMode");
const explorerSection = document.getElementById("explorer");
const info = document.getElementById("clusterInfo");
const details = document.getElementById("clusterDetails");
const processLegend = document.getElementById("processLegend");
const lineageLegend = document.getElementById("lineageLegend");
const status = document.getElementById("matrixStatus");
const svg = document.getElementById("treeSvg");
const headerCanvas = document.getElementById("columnHeaderCanvas");
const headerCtx = headerCanvas.getContext("2d");
const canvas = document.getElementById("matrixCanvas");
const ctx = canvas.getContext("2d");
const baseCanvas = document.createElement("canvas");
const baseCtx = baseCanvas.getContext("2d");
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
const LAYOUT_KEY = "recore.explorer.layout.v1";
const DEFAULT_TREE_W = 320;
const DEFAULT_SYNTENY_W_PCT = 38;
const MIN_TREE_W = 180;
const MIN_MATRIX_W = 220;
const MIN_SYNTENY_W = 260;
const LINEAGE_STRIP_W = 6;
// Single source for the zoom bounds. The sliders, currentCellSize/currentRowHeight
// and the wheel handler must agree, or the wheel snaps the value into its own range.
const MIN_COL_ZOOM = 0.25;
const MAX_COL_ZOOM = 44;
const MIN_ROW_ZOOM = 1;
const MAX_ROW_ZOOM = 18;
// Column names are only worth drawing once a cell is wide enough to carry a
// rotated label without the text of neighbouring columns touching.
const LABEL_MIN_CELL_W = 9;
const LABEL_BAND_H = 82;
const LABEL_FONT = "10px Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif";
const defaultStatus = "Wheel over tree zooms Y; wheel over the column header zooms columns. In the matrix body, Ctrl/Option + wheel zooms Y and Shift + wheel zooms columns.";
let selected = 0;
let hoverRow = -1;
let hoverCol = -1;
let leafOrder = [];
let visibleClusters = [];
let genomeStats = [];
let maxGenomeStats = {R:1, N:1, A:1, correction:1};
let maxDepth = 0;
let maxDistance = 0;
let treeHeight = 0;
let treeWidth = DEFAULT_TREE_W;
const headerH = 24;
const maxCanvasDim = 32000;
const maxCanvasPixels = 90000000;
let cellH = 12;
let cellW = 8;
let baseCellH = 0;
let baseCellW = 0;
let baseValid = false;
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
let treeRoot = data.tree;

/* ---- Cell decoding -------------------------------------------------------
   Cells arrive as one base64 blob per cluster holding the genome indices that
   are not absent, split by state (see cellcodec.rs). Decoding every cluster up
   front would rebuild the dense matrix we are trying to avoid, so rows are
   decoded on demand and kept in a bounded cache: rendering only ever touches
   the columns currently on screen. */
const STATE_CODES = ["P", "R", "N", "A"];
const CELL_CACHE_LIMIT = 4096;
const cellCache = new Map();

function decodeBase64(text) {
  const binary = atob(text);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i += 1) bytes[i] = binary.charCodeAt(i);
  return bytes;
}
function decodeClusterRow(encoded, nGenomes) {
  const row = new Uint8Array(nGenomes).fill(3);
  if (!encoded) return row;
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
  for (let state = 0; state < 3; state += 1) {
    const mode = bytes[pos];
    pos += 1;
    if (mode === 0 || mode === 1) {
      const count = readVarint();
      const indices = new Int32Array(count);
      let previous = 0;
      for (let i = 0; i < count; i += 1) {
        previous += readVarint();
        indices[i] = previous;
      }
      if (mode === 0) {
        for (let i = 0; i < count; i += 1) row[indices[i]] = state;
      } else {
        const absent = new Uint8Array(nGenomes);
        for (let i = 0; i < count; i += 1) absent[indices[i]] = 1;
        for (let g = 0; g < nGenomes; g += 1) if (!absent[g]) row[g] = state;
      }
    } else if (mode === 2) {
      const len = Math.ceil(nGenomes / 8);
      for (let g = 0; g < nGenomes; g += 1) {
        if (bytes[pos + (g >> 3)] & (1 << (g & 7))) row[g] = state;
      }
      pos += len;
    } else {
      return row;
    }
  }
  return row;
}
function clusterRow(clusterIdx) {
  const cached = cellCache.get(clusterIdx);
  if (cached) {
    // Refresh recency so the working set survives a pan across the matrix.
    cellCache.delete(clusterIdx);
    cellCache.set(clusterIdx, cached);
    return cached;
  }
  const row = decodeClusterRow(data.cells[clusterIdx], data.genomes.length);
  cellCache.set(clusterIdx, row);
  if (cellCache.size > CELL_CACHE_LIMIT) {
    cellCache.delete(cellCache.keys().next().value);
  }
  return row;
}
function stateAt(clusterIdx, genomeIdx) {
  const row = clusterRow(clusterIdx);
  return STATE_CODES[row ? row[genomeIdx] : 3] || "A";
}
function countsForCluster(clusterIdx) {
  const counts = {P:0, R:0, N:0, A:0};
  const row = clusterRow(clusterIdx);
  for (let i = 0; i < data.genomes.length; i += 1) {
    counts[STATE_CODES[row[i]]] += 1;
  }
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
function refreshVisibleClusters() {
  const previousSelected = selected;
  visibleClusters = data.clusters
    .map((_, idx) => idx)
    .filter(idx => passesActiveFilters(idx, true));
  if (!visibleClusters.length) {
    if (referenceOnly) {
      referenceOnly = false;
      visibleClusters = data.clusters
        .map((_, idx) => idx)
        .filter(idx => passesActiveFilters(idx, true));
    }
    if (!visibleClusters.length && mechanismFilters.size) {
      mechanismFilters.clear();
      visibleClusters = data.clusters.map((_, idx) => idx);
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
  if (referenceButton) {
    referenceButton.classList.toggle("active", referenceOnly);
    referenceButton.textContent = data.referenceGenome
      ? `Reference-absent only (${data.referenceGenome.name})`
      : "Reference-absent only";
    referenceButton.disabled = !data.referenceGenome;
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
function populateOptions() {
  const previousSelected = selected;
  const q = search.value.trim().toLowerCase();
  select.innerHTML = "";
  visibleClusters.forEach(idx => {
    const cluster = data.clusters[idx];
    const haystack = (cluster.name + " " + (cluster.annotation || "")).toLowerCase();
    if (q && !haystack.includes(q)) return;
    const opt = document.createElement("option");
    opt.value = String(idx);
    opt.textContent = cluster.annotation ? `${cluster.name} - ${cluster.annotation}` : cluster.name;
    select.appendChild(opt);
  });
  if (![...select.options].some(opt => Number(opt.value) === selected)) {
    selected = select.options.length ? Number(select.options[0].value) : 0;
  }
  if (selected !== previousSelected) {
    selectedSyntenyBlock = 0;
    selectedSyntenyReference = false;
    selectedSyntenyCarriers = new Set();
    syntenyKey = "";
  }
  select.value = String(selected);
}
function currentCellSize() {
  const value = Number(sizeControl.value || 8);
  return Math.max(MIN_COL_ZOOM, Math.min(MAX_COL_ZOOM, value));
}
function currentRowHeight() {
  const value = Number(rowHeightControl.value || 8);
  return Math.max(MIN_ROW_ZOOM, Math.min(MAX_ROW_ZOOM, value));
}
function updateCellGeometry() {
  const rows = leafOrder.length;
  const cols = visibleClusters.length;
  const colSize = currentCellSize();
  const maxWByDim = Math.max(1, Math.floor(maxCanvasDim / Math.max(1, cols)));
  const maxHByDim = Math.max(1, Math.floor(maxCanvasDim / Math.max(1, rows)));
  cellH = Math.min(currentRowHeight(), maxHByDim);
  cellW = Math.min(Math.max(1, Math.round(colSize * 0.82)), maxWByDim);
  const pixelBudget = rows * cols * cellH * cellW;
  if (pixelBudget > maxCanvasPixels) {
    const scale = Math.sqrt(maxCanvasPixels / pixelBudget);
    cellH = Math.max(1, Math.floor(cellH * scale));
    cellW = Math.max(1, Math.floor(cellW * scale));
  }
  treeHeight = Math.max(70, Math.ceil(rows * cellH));
  return {rows, cols};
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
function xForNode(node) {
  const labelPad = labelToggle.checked && leafOrder.length <= 1500 ? 150 : 24;
  const usable = Math.max(80, treeWidth - 22 - labelPad);
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
  search.value = "";
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
function syntenyLengthScale(rows) {
  const lengths = rows
    .flatMap(row => row.genes || [])
    .map(syntenyLength)
    .filter(value => value > 0)
    .sort((a, b) => a - b);
  if (!lengths.length) return 1200;
  const idx = Math.min(lengths.length - 1, Math.floor(lengths.length * 0.90));
  return Math.max(300, lengths[idx]);
}
function syntenyGeneWidth(gene, scaleLen) {
  const relative = Math.min(1, syntenyLength(gene) / Math.max(1, scaleLen));
  return (34 + relative * 132) * syntenyZoom();
}
function syntenyRowLayout(genes, left, gap, scaleLen) {
  let x = left;
  return genes.map(gene => {
    const width = syntenyGeneWidth(gene, scaleLen);
    const item = {gene, x, width};
    x += width + gap;
    return item;
  });
}
function arrowPath(x, y, w, h, strand = "+") {
  const head = Math.min(12, Math.max(7, w * 0.24));
  if (strand === "-") {
    return `M ${x + w} ${y} H ${x + head} L ${x} ${y + h / 2} L ${x + head} ${y + h} H ${x + w} Z`;
  }
  return `M ${x} ${y} H ${x + w - head} L ${x + w} ${y + h / 2} L ${x + w - head} ${y + h} H ${x} Z`;
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

  // Scroll the SVG so the target column of the first block is visible.
  const firstLayout = rowLayouts[0] || [];
  const firstTarget = firstLayout.find(item => item.gene.target);
  if (firstTarget && scrollContainer) {
    const targetX = firstTarget.x + firstTarget.width / 2;
    scrollContainer.scrollLeft = Math.max(0, targetX - scrollContainer.clientWidth / 2);
  }
  if (scrollContainer) scrollContainer.scrollTop = 0;
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
      syntenyPct: Math.round(layoutState.syntenyPct * 10) / 10,
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
    syntenyPct: Number.isFinite(stored.syntenyPct)
      ? Math.min(70, Math.max(20, stored.syntenyPct))
      : DEFAULT_SYNTENY_W_PCT,
    syntenyHiddenStored: stored.syntenyHidden,
  };
})();
function applyLayout() {
  if (!viewerSplit || !viewerBody) return;
  const bodyW = viewerBody.clientWidth || 0;
  let treeW = layoutState.treeW;
  let syntenyPct = layoutState.syntenyPct;
  if (bodyW > 0) {
    const availableForMatrix = bodyW - treeW - (syntenyHidden ? 0 : (bodyW * syntenyPct) / 100) - 12;
    if (availableForMatrix < MIN_MATRIX_W) {
      const overshoot = MIN_MATRIX_W - availableForMatrix;
      const shrinkTree = Math.min(treeW - MIN_TREE_W, overshoot * 0.5);
      treeW = Math.max(MIN_TREE_W, treeW - shrinkTree);
      const remaining = overshoot - shrinkTree;
      if (!syntenyHidden && remaining > 0) {
        const syntenyW = (bodyW * syntenyPct) / 100 - remaining;
        const minPct = (MIN_SYNTENY_W / bodyW) * 100;
        syntenyPct = Math.max(minPct, (syntenyW / bodyW) * 100);
      }
    }
    layoutState.treeW = treeW;
    layoutState.syntenyPct = syntenyPct;
  }
  viewerBody.style.setProperty("--tree-w", `${treeW}px`);
  viewerBody.style.setProperty("--synteny-w", `${syntenyPct.toFixed(2)}%`);
  if (viewerSplit) {
    viewerSplit.style.setProperty("--tree-w", `${treeW}px`);
    viewerSplit.style.setProperty("--synteny-w", `${syntenyPct.toFixed(2)}%`);
  }
  document.querySelectorAll(".synteny-control").forEach(control => {
    control.style.display = syntenyHidden ? "none" : "";
  });
  treeWidth = Math.max(180, Math.round(treeW) - LINEAGE_STRIP_W - 4);
}
function applyViewerSplit() {
  applyLayout();
}
function applySyntenyVisibility() {
  if (!viewerBody || !syntenyToggle) return;
  viewerBody.classList.toggle("synteny-hidden", syntenyHidden);
  if (viewerSplit) viewerSplit.classList.toggle("synteny-hidden", syntenyHidden);
  syntenyToggle.classList.toggle("active", !syntenyHidden);
  syntenyToggle.textContent = syntenyHidden ? "Show syntenies" : "Hide syntenies";
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
function setupDragHandle(handle, side) {
  if (!handle) return;
  let startX = 0;
  let startTreeW = 0;
  let startSyntenyW = 0;
  let bodyW = 0;
  let active = false;
  const onMove = event => {
    if (!active) return;
    const dx = event.clientX - startX;
    if (side === "tree") {
      const maxTree = Math.max(MIN_TREE_W, bodyW - MIN_MATRIX_W - (syntenyHidden ? 0 : startSyntenyW) - 12);
      layoutState.treeW = Math.max(MIN_TREE_W, Math.min(maxTree, startTreeW + dx));
    } else {
      const newSynW = Math.max(MIN_SYNTENY_W, startSyntenyW - dx);
      const minPct = (MIN_SYNTENY_W / bodyW) * 100;
      const maxPct = Math.max(minPct, ((bodyW - layoutState.treeW - MIN_MATRIX_W - 12) / bodyW) * 100);
      layoutState.syntenyPct = Math.min(maxPct, Math.max(minPct, (newSynW / bodyW) * 100));
    }
    applyLayout();
    scheduleRender(false);
  };
  const onUp = event => {
    if (!active) return;
    active = false;
    handle.classList.remove("dragging");
    document.body.classList.remove("recore-dragging");
    persistLayout();
    if (handle.hasPointerCapture && handle.hasPointerCapture(event.pointerId)) {
      handle.releasePointerCapture(event.pointerId);
    }
    treeDirty = true;
    resetMatrixCache();
    render(false);
  };
  handle.addEventListener("pointerdown", event => {
    if (event.button !== 0) return;
    bodyW = viewerBody.clientWidth || 0;
    if (bodyW <= 0) return;
    active = true;
    startX = event.clientX;
    startTreeW = layoutState.treeW;
    startSyntenyW = (bodyW * layoutState.syntenyPct) / 100;
    handle.classList.add("dragging");
    document.body.classList.add("recore-dragging");
    try { handle.setPointerCapture(event.pointerId); } catch (_e) {}
    event.preventDefault();
  });
  handle.addEventListener("pointermove", onMove);
  handle.addEventListener("pointerup", onUp);
  handle.addEventListener("pointercancel", onUp);
  handle.addEventListener("keydown", event => {
    const step = event.shiftKey ? 32 : 12;
    if (event.key === "ArrowLeft" || event.key === "ArrowRight") {
      bodyW = viewerBody.clientWidth || 0;
      const dx = (event.key === "ArrowRight" ? 1 : -1) * step;
      if (side === "tree") {
        const maxTree = Math.max(MIN_TREE_W, bodyW - MIN_MATRIX_W - (syntenyHidden ? 0 : (bodyW * layoutState.syntenyPct) / 100) - 12);
        layoutState.treeW = Math.max(MIN_TREE_W, Math.min(maxTree, layoutState.treeW + dx));
      } else {
        const newSynW = ((bodyW * layoutState.syntenyPct) / 100) - dx;
        const minPct = (MIN_SYNTENY_W / bodyW) * 100;
        const maxPct = Math.max(minPct, ((bodyW - layoutState.treeW - MIN_MATRIX_W - 12) / bodyW) * 100);
        layoutState.syntenyPct = Math.min(maxPct, Math.max(minPct, (newSynW / bodyW) * 100));
      }
      applyLayout();
      persistLayout();
      treeDirty = true;
      resetMatrixCache();
      render(false);
      event.preventDefault();
    }
  });
}
function renderLineageStrip() {
  if (!lineageStripCanvas || !lineageStripCtx) return;
  const rows = leafOrder.length;
  const height = Math.max(1, Math.ceil(rows * cellH));
  sizeCanvas(
    lineageStripCanvas,
    lineageStripCtx,
    LINEAGE_STRIP_W,
    height,
    canvasScale(LINEAGE_STRIP_W, height),
  );
  lineageStripCtx.clearRect(0, 0, LINEAGE_STRIP_W, height);
  for (let r = 0; r < rows; r += 1) {
    const genomeIdx = leafOrder[r];
    const lineage = lineageOfGenomeName(data.genomes[genomeIdx] || "");
    lineageStripCtx.fillStyle = lineageColor(lineage);
    lineageStripCtx.fillRect(0, r * cellH, LINEAGE_STRIP_W, Math.max(1, cellH));
  }
}
function treeCaption() {
  const rootLabel = hasL8Tip(treeRoot) ? "L8-rooted" : "input-rooted";
  return treeMode === "phylogram" && maxDistance > 0
    ? `phylogram, ${rootLabel}, longest root-to-tip path ${formatBranchLength(maxDistance)}`
    : `cladogram, ${rootLabel}, branch lengths ignored`;
}
function updateTopBarMeta() {
  if (topTreeMeta) {
    topTreeMeta.textContent = `${leafOrder.length} matched tips · ${treeCaption()}`;
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
function renderTree() {
  if (treeDirty) {
    const parts = [];
    drawNode(treeRoot, parts);
    svg.setAttribute("width", String(treeWidth));
    svg.setAttribute("height", String(treeHeight));
    svg.setAttribute("viewBox", `0 0 ${treeWidth} ${treeHeight}`);
    svg.innerHTML = `<line id="treeHoverGuide" x1="0" y1="0" x2="${treeWidth}" y2="0" class="tree-row-guide" style="display:none"/>${parts.join("")}`;
    treeDirty = false;
  }
  const guide = document.getElementById("treeHoverGuide");
  if (!guide) return;
  if (hoverRow >= 0) {
    const hy = hoverRow * cellH + cellH / 2;
    guide.setAttribute("y1", String(hy));
    guide.setAttribute("y2", String(hy));
    guide.style.display = "block";
  } else {
    guide.style.display = "none";
  }
}
function columnLabelBand() {
  return cellW >= LABEL_MIN_CELL_W && visibleClusters.length <= 400 ? LABEL_BAND_H : 0;
}
function renderColumnHeader() {
  const cols = visibleClusters.length;
  const selectedCol = selectedVisibleColumn();
  const band = columnLabelBand();
  const cssW = Math.max(1, cols * cellW);
  const cssH = headerH + band;
  sizeCanvas(headerCanvas, headerCtx, cssW, cssH, canvasScale(cssW, cssH));
  headerCtx.clearRect(0, 0, cssW, cssH);
  headerCtx.fillStyle = headerBg;
  headerCtx.fillRect(0, 0, cssW, cssH);
  if (matrixColorMode.value === "process") {
    visibleClusters.forEach((clusterIdx, c) => {
      headerCtx.fillStyle = processColor(processKind(data.clusters[clusterIdx]));
      headerCtx.fillRect(c * cellW, band, Math.max(1, cellW), 7);
    });
  }
  if (hoverCol >= 0) {
    headerCtx.fillStyle = hoverColTint;
    headerCtx.fillRect(hoverCol * cellW, 0, cellW, cssH);
  }
  if (selectedCol >= 0) {
    headerCtx.fillStyle = selectColTint;
    headerCtx.fillRect(selectedCol * cellW, 0, cellW, cssH);
    headerCtx.strokeStyle = selectStroke;
    headerCtx.strokeRect(selectedCol * cellW + 0.5, 0.5, Math.max(1, cellW - 1), cssH - 1);
  }
  for (let c = 0; c < cols; c += 1) {
    const clusterIdx = originalClusterAtVisibleColumn(c);
    headerCtx.fillStyle = matrixColorMode.value === "process"
      ? processColor(processKind(data.clusters[clusterIdx]))
      : gridStroke;
    headerCtx.fillRect(c * cellW, cssH - 9, Math.max(1, cellW), 6);
  }
  if (!band) return;
  headerCtx.save();
  headerCtx.font = LABEL_FONT;
  headerCtx.textBaseline = "middle";
  // After rotate(-90deg) the local +x axis points up the screen, so left-aligned
  // text grows upward from the tick row and reads bottom to top.
  headerCtx.textAlign = "left";
  for (let c = 0; c < cols; c += 1) {
    const cluster = data.clusters[originalClusterAtVisibleColumn(c)];
    if (!cluster) continue;
    headerCtx.fillStyle = c === selectedCol ? selectStroke : labelInk;
    headerCtx.save();
    headerCtx.translate(c * cellW + cellW / 2, band - 6);
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
/* Canvases are laid out in CSS pixels but were also allocated in CSS pixels, so on
   any HiDPI screen the browser upscaled a 1x bitmap and the matrix came out soft.
   Allocate at the device ratio instead, capped by the same dimension/area budget
   that keeps very large matrices inside browser canvas limits. */
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
let matrixScale = 1;
function rebuildBaseMatrix() {
  const {rows, cols} = updateCellGeometry();
  const cssW = Math.max(1, cols * cellW);
  const cssH = Math.max(1, Math.ceil(rows * cellH));
  matrixScale = canvasScale(cssW, cssH);
  baseCanvas.width = Math.max(1, Math.round(cssW * matrixScale));
  baseCanvas.height = Math.max(1, Math.round(cssH * matrixScale));
  baseCtx.setTransform(matrixScale, 0, 0, matrixScale, 0, 0);
  baseCtx.clearRect(0, 0, cssW, cssH);
  baseCtx.fillStyle = panelBg;
  baseCtx.fillRect(0, 0, cssW, cssH);
  for (let r = 0; r < rows; r += 1) {
    const genomeIdx = leafOrder[r];
    for (let c = 0; c < cols; c += 1) {
      const clusterIdx = originalClusterAtVisibleColumn(c);
      const state = stateAt(clusterIdx, genomeIdx);
      if (matrixColorMode.value === "process" && state !== "A") {
        baseCtx.fillStyle = processColor(processKind(data.clusters[clusterIdx]));
      } else {
        baseCtx.fillStyle = colors[state] || colors.A;
      }
      const cellDrawH = cellH >= 2 ? Math.max(1, cellH - 1) : Math.max(0.2, cellH);
      baseCtx.fillRect(c * cellW, r * cellH, Math.max(1, cellW - 1), cellDrawH);
    }
  }
  baseCellH = cellH;
  baseCellW = cellW;
  baseValid = true;
}
function renderMatrix(scrollToSelected) {
  const {rows, cols} = updateCellGeometry();
  const selectedCol = selectedVisibleColumn();
  if (!baseValid || baseCellH !== cellH || baseCellW !== cellW) {
    rebuildBaseMatrix();
  }
  const cssW = Math.max(1, cols * cellW);
  const cssH = Math.max(1, Math.ceil(rows * cellH));
  sizeCanvas(canvas, ctx, cssW, cssH, matrixScale);
  ctx.clearRect(0, 0, cssW, cssH);
  ctx.drawImage(baseCanvas, 0, 0, cssW, cssH);
  if (selectedCol >= 0) {
    ctx.fillStyle = "rgba(37,99,235,.13)";
    ctx.fillRect(selectedCol * cellW, 0, cellW, cssH);
  }
  if (hoverRow >= 0 && hoverRow < rows) {
    ctx.fillStyle = "rgba(245,158,11,.10)";
    ctx.fillRect(0, hoverRow * cellH, cssW, Math.max(0.25, cellH));
  }
  if (hoverCol >= 0 && hoverCol < cols) {
    ctx.fillStyle = "rgba(245,158,11,.08)";
    ctx.fillRect(hoverCol * cellW, 0, cellW, cssH);
  }
  if (hoverCol >= 0) {
    ctx.strokeStyle = "#b45309";
    ctx.lineWidth = 1;
    ctx.strokeRect(hoverCol * cellW + 0.5, 0.5, Math.max(1, cellW - 1), Math.max(1, cssH - 1));
  }
  if (hoverRow >= 0) {
    ctx.strokeStyle = "#b45309";
    ctx.lineWidth = 1;
    ctx.strokeRect(0.5, hoverRow * cellH + 0.5, Math.max(1, cssW - 1), Math.max(0.25, cellH - 0.1));
  }
  ctx.strokeStyle = selectStroke;
  ctx.lineWidth = Math.max(1, Math.min(2, cellW));
  if (selectedCol >= 0) {
    ctx.strokeRect(selectedCol * cellW + 0.5, 0.5, Math.max(1, cellW - 1), Math.max(1, cssH - 1));
  }
  if (scrollToSelected && selectedCol >= 0) {
    matrixPane.scrollLeft = Math.max(0, selectedCol * cellW - 80);
    headerScroller.scrollLeft = matrixPane.scrollLeft;
  }
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
    info.textContent = `${cluster.present}/${data.genomes.length} genomes (${Math.round(cluster.prevalence * 1000) / 10}%)`;
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
    details.innerHTML = `
      <div class="cluster-head">
        <h3>${escapeHtml(cluster.name)}</h3>
        <p class="annotation">${escapeHtml(cluster.annotation || "No annotation provided")}</p>
        <span class="prevalence num">${cluster.present} / ${data.genomes.length} genomes &middot; ${prevalencePct}</span>
      </div>
      <dl class="cluster-dl">
        <dt>Evidence states</dt>
        <dd><div class="state-counts">
          <span class="state-chip input"><strong>${counts.P}</strong> input</span>
          <span class="state-chip rescued"><strong>${counts.R}</strong> rescued</span>
          <span class="state-chip nsc"><strong>${counts.N}</strong> NSC</span>
          <span class="state-chip absent"><strong>${counts.A}</strong> absent</span>
        </div></dd>
        <dt>Mechanism</dt>
        <dd>
          <span class="process-chip" style="--process-color:${mechanismColor}"><i></i>${escapeHtml(mechanismLabel)}</span>
          <div class="process-wording" style="margin-top:4px">${escapeHtml(mechanismDescription)}</div>
        </dd>
        ${refRow}
        ${supportRow}
        ${phyloRow}
        ${detailRow}
        ${wordingRow}
      </dl>
    `;
    renderProcessLegend();
    detailsKey = nextDetailsKey;
  } else if (nextDetailsKey !== detailsKey) {
    info.textContent = "";
    details.innerHTML = "";
    processLegend.innerHTML = "";
    detailsKey = nextDetailsKey;
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
const MIN_VIEWER_H = 240;
const MAX_VIEWER_H = 920;
const FIT_TARGET_H = 560;
const FIT_TARGET_W = 720;
// Small trees left most of the viewport empty because the panes always claimed
// their maximum height. Size the viewer to whichever pane needs the most room,
// capped so a 5000-tip tree still scrolls inside a sane window.
function applyViewerHeight() {
  if (!viewerBody) return;
  const contentH = Math.max(treeHeight, leafOrder.length * cellH) + 8;
  const syntenyH = syntenyHidden ? 0 : (syntenySvg ? syntenySvg.getBoundingClientRect().height + 46 : 0);
  const viewportCap = Math.min(MAX_VIEWER_H, Math.round(window.innerHeight * 0.76));
  const needed = Math.max(contentH, syntenyH);
  const height = `${Math.round(Math.max(MIN_VIEWER_H, Math.min(viewportCap, needed)))}px`;
  if (viewerBody.style.getPropertyValue("--viewer-h") !== height) {
    viewerBody.style.setProperty("--viewer-h", height);
  }
}
function render(scrollToSelected = false) {
  renderMatrix(scrollToSelected);
  renderColumnHeader();
  renderSyntenyPanel();
  renderTree();
  renderLineageStrip();
  renderLineageLegend();
  updateTopBarMeta();
  applyViewerHeight();
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
  baseValid = false;
  hoverCol = -1;
  treeDirty = true;
}
refreshThemeColors();
treeRoot = rootByL8(data.tree);
leafOrder = [];
maxDepth = 0;
maxDistance = 0;
layout(treeRoot, 0, 0);
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
    ? Math.max(220, bodyW - layoutState.treeW - (syntenyHidden ? 0 : (bodyW * layoutState.syntenyPct) / 100) - 24)
    : FIT_TARGET_W;
  const rowFit = Math.max(MIN_ROW_ZOOM, Math.min(MAX_ROW_ZOOM, Math.floor(FIT_TARGET_H / tips)));
  // cellW is round(slider * 0.82), so undo that factor when fitting to a width.
  const colFit = Math.max(MIN_COL_ZOOM, Math.min(MAX_COL_ZOOM, Math.floor((matrixW / cols) / 0.82)));
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
setupDragHandle(syntenyHandle, "synteny");
if (typeof ResizeObserver !== "undefined") {
  // Width only: applyViewerHeight writes back to this element's height, so
  // reacting to height changes would make the observer chase its own tail.
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
    resetMatrixCache();
    render(false);
    return;
  }
  applyViewerHeight();
});
render(true);
window.addEventListener("recore:theme-change", () => {
  refreshThemeColors();
  resetMatrixCache();
  render(false);
});
search.addEventListener("input", () => { populateOptions(); treeDirty = true; render(true); });
select.addEventListener("change", () => { selected = Number(select.value || 0); selectedSyntenyBlock = 0; selectedSyntenyReference = false; syntenyKey = ""; treeDirty = true; render(true); });
colorMode.addEventListener("change", () => { treeDirty = true; render(false); });
labelToggle.addEventListener("change", () => { treeDirty = true; render(false); });
phylogramButton.addEventListener("click", () => setTreeMode("phylogram"));
cladogramButton.addEventListener("click", () => setTreeMode("cladogram"));
referenceButton.addEventListener("click", () => {
  if (!data.referenceGenome) return;
  referenceOnly = !referenceOnly;
  refreshVisibleClusters();
  populateOptions();
  resetMatrixCache();
  render(true);
});
syntenyToggle.addEventListener("click", () => {
  syntenyHidden = !syntenyHidden;
  applySyntenyVisibility();
  render(false);
});
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
  render(true);
});
fullscreenButton.addEventListener("click", async () => {
  if (!document.fullscreenElement && explorerSection.requestFullscreen) {
    await explorerSection.requestFullscreen();
  } else if (document.exitFullscreen) {
    await document.exitFullscreen();
  }
});
document.addEventListener("fullscreenchange", () => {
  fullscreenButton.textContent = document.fullscreenElement === explorerSection
    ? "Exit fullscreen"
    : "Fullscreen";
});
sizeControl.addEventListener("input", () => { resetMatrixCache(); render(false); });
rowHeightControl.addEventListener("input", () => { resetMatrixCache(); render(false); });
syntenyZoomControl.addEventListener("input", () => { syntenyKey = ""; render(false); });
matrixColorMode.addEventListener("change", () => { resetMatrixCache(); render(false); });
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
svg.addEventListener("wheel", event => wheelZoom(event, "row"), {passive:false});
headerScroller.addEventListener("wheel", event => wheelZoom(event, "column"), {passive:false});
headerCanvas.addEventListener("wheel", event => wheelZoom(event, "column"), {passive:false});
matrixPane.addEventListener("wheel", wheelZoom, {passive:false});
svg.addEventListener("mousemove", event => {
  const rect = svg.getBoundingClientRect();
  const r = Math.floor((event.clientY - rect.top) / cellH);
  if (r < 0 || r >= leafOrder.length) return;
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
svg.addEventListener("mouseleave", () => {
  setHover(-1, -1, defaultStatus);
});
headerCanvas.addEventListener("mousemove", event => {
  const rect = headerCanvas.getBoundingClientRect();
  const c = Math.floor((event.clientX - rect.left) / cellW);
  if (c < 0 || c >= visibleClusters.length) return;
  const clusterIdx = originalClusterAtVisibleColumn(c);
  const cluster = data.clusters[clusterIdx];
  setHover(-1, c, `${cluster.name}: ${cluster.annotation || "No annotation provided"}`);
});
headerCanvas.addEventListener("mouseleave", () => {
  setHover(-1, -1, defaultStatus);
});
headerCanvas.addEventListener("click", event => {
  const rect = headerCanvas.getBoundingClientRect();
  const c = Math.floor((event.clientX - rect.left) / cellW);
  if (c < 0 || c >= visibleClusters.length) return;
  selected = originalClusterAtVisibleColumn(c);
  selectedSyntenyReference = false;
  search.value = "";
  populateOptions();
  treeDirty = true;
  render(true);
});
canvas.addEventListener("mousemove", event => {
  const rect = canvas.getBoundingClientRect();
  const c = Math.floor((event.clientX - rect.left) / cellW);
  const r = Math.floor((event.clientY - rect.top) / cellH);
  if (c < 0 || r < 0 || c >= visibleClusters.length || r >= leafOrder.length) return;
  const clusterIdx = originalClusterAtVisibleColumn(c);
  const genomeIdx = leafOrder[r];
  const state = stateAt(clusterIdx, genomeIdx);
  const cluster = data.clusters[clusterIdx];
  const mechanism = processKind(cluster);
  setHover(r, c, `${cluster.name} / ${data.genomes[genomeIdx]}: ${labels[state]}; ${processLabel(mechanism)}`);
  const paneRect = matrixPane.getBoundingClientRect();
  tooltip.style.display = "block";
  tooltip.style.left = `${event.clientX - paneRect.left + matrixPane.scrollLeft + 12}px`;
  tooltip.style.top = `${event.clientY - paneRect.top + 12}px`;
  tooltip.innerHTML = `<strong>${escapeHtml(cluster.name)}</strong><br>${escapeHtml(data.genomes[genomeIdx])}: ${escapeHtml(labels[state])}<br><span>${escapeHtml(processLabel(mechanism))}</span><br><span>${escapeHtml(cluster.annotation || "")}</span>`;
});
canvas.addEventListener("mouseleave", () => {
  tooltip.style.display = "none";
  setHover(-1, -1, defaultStatus);
});
canvas.addEventListener("click", event => {
  const rect = canvas.getBoundingClientRect();
  const c = Math.floor((event.clientX - rect.left) / cellW);
  if (c < 0 || c >= visibleClusters.length) return;
  selected = originalClusterAtVisibleColumn(c);
  selectedSyntenyReference = false;
  search.value = "";
  populateOptions();
  treeDirty = true;
  render(true);
});
})();
