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
const splitControl = document.getElementById("viewerSplitSize");
const syntenyZoomControl = document.getElementById("syntenyZoom");
const fullscreenButton = document.getElementById("explorerFullscreen");
const phylogramButton = document.getElementById("phylogramMode");
const cladogramButton = document.getElementById("cladogramMode");
const explorerSection = document.getElementById("matrixExplorer");
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
const headerScroller = document.getElementById("headerScroller");
const viewerBody = document.getElementById("viewerBody");
const viewerSplit = document.getElementById("viewerSplit");
const tooltip = document.getElementById("matrixTooltip");
const syntenyPanel = document.getElementById("syntenyPanel");
const syntenySummary = document.getElementById("syntenySummary");
const syntenySvg = document.getElementById("syntenySvg");
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
const treeWidth = 410;
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

function stateAt(clusterIdx, genomeIdx) {
  return (data.states[clusterIdx] && data.states[clusterIdx][genomeIdx]) || "A";
}
function countsForCluster(clusterIdx) {
  const counts = {P:0, R:0, N:0, A:0};
  const row = data.states[clusterIdx] || "";
  for (let i = 0; i < data.genomes.length; i += 1) {
    const state = row[i] || "A";
    counts[state] = (counts[state] || 0) + 1;
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
  genomeStats = data.genomes.map(() => ({P:0, R:0, N:0, A:0}));
  data.states.forEach(row => {
    for (let gi = 0; gi < data.genomes.length; gi += 1) {
      const state = row[gi] || "A";
      genomeStats[gi][state] = (genomeStats[gi][state] || 0) + 1;
    }
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
    return `<button type="button" class="process-item${active}${filtered}" data-mechanism="${escapeHtml(mechanism)}"><span class="process-swatch" style="background:${processColor(mechanism)}"></span><div><strong>${escapeHtml(processLabel(mechanism))}</strong><small>${row.clusters} clusters; ${row.present} present cells. ${escapeHtml(processDescription(mechanism))}</small></div></button>`;
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
  return Math.max(0.25, Math.min(18, value));
}
function currentRowHeight() {
  const value = Number(rowHeightControl.value || 8);
  return Math.max(1, Math.min(18, value));
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
function drawNode(node, parts, xScale) {
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
      drawNode(child, parts, xScale);
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
  const nextKey = cluster ? `${cluster.name}|${blocks.length}|${selectedSyntenyBlock}|${selectedSyntenyReference}|${syntenyZoom()}` : "none";
  if (nextKey === syntenyKey) return;
  syntenyKey = nextKey;
  treeDirty = true;
  if (!cluster) {
    selectedSyntenyCarriers = new Set();
    syntenySummary.textContent = "No cluster selected.";
    syntenySvg.setAttribute("width", "620");
    syntenySvg.setAttribute("height", "80");
    syntenySvg.innerHTML = `<text x="14" y="36" class="synteny-label">No cluster selected.</text>`;
    return;
  }
  if (!blocks.length) {
    selectedSyntenyCarriers = new Set();
    syntenySummary.textContent = `${cluster.name}: no 5-neighbour synteny blocks in this run.`;
    syntenySvg.setAttribute("width", "620");
    syntenySvg.setAttribute("height", "88");
    syntenySvg.innerHTML = `<text x="14" y="38" class="synteny-label">No unique synteny block data for ${escapeHtml(cluster.name)}.</text>`;
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
  const rowH = 94;
  const geneH = 20;
  const gap = 12 * syntenyZoom();
  const left = 176;
  const top = 58;
  const lengthScale = syntenyLengthScale(rows);
  const rowLayouts = rows.map(row => syntenyRowLayout(row.genes, left, gap, lengthScale));
  const width = Math.max(
    740,
    ...rowLayouts.map(layout => {
      const last = layout[layout.length - 1];
      return last ? last.x + last.width + 30 : left + 120;
    }),
  );
  const height = top + rows.length * rowH + 24;
  const parts = [];
  let previousPositions = null;
  const selectedBlock = selectedSyntenyBlockData();
  const selectedCount = selectedSyntenyReference
    ? selectedSyntenyCarriers.size
    : selectedBlock ? syntenyCarrierIndices(selectedBlock).length : 0;
  const selectedLabel = selectedSyntenyReference
    ? `${data.referenceGenome ? data.referenceGenome.name : "reference"} context`
    : `#${selectedSyntenyBlock + 1}`;
  const referenceSummary = hasReferenceRow
    ? " Row 2 is the H37Rv/reference context inferred from shared neighbours."
    : " H37Rv/reference has the same displayed gene architecture as the representative block and is merged visually.";
  syntenySummary.textContent = `${cluster.name}: ${blocks.length} unique 5-neighbour block${blocks.length === 1 ? "" : "s"}; selected ${selectedLabel} marks ${selectedCount} tree tips.${referenceSummary}`;
  rows.forEach((row, rowIdx) => {
    const genes = row.genes;
    const layout = rowLayouts[rowIdx];
    const y = top + rowIdx * rowH;
    const positions = new Map();
    let targetCenter = null;
    layout.forEach(item => {
      const gene = item.gene;
      const center = {x:item.x + item.width / 2, y:y + geneH / 2};
      const key = syntenyLinkKey(gene);
      if (key && !positions.has(key)) positions.set(key, center);
      if (gene.target) targetCenter = center;
    });
    positions.forEach((pos, name) => {
      const previous = previousPositions && previousPositions.get(name);
      if (!previous) return;
      const midY = (previous.y + pos.y) / 2;
      parts.push(`<path d="M ${previous.x} ${previous.y} C ${previous.x} ${midY} ${pos.x} ${midY} ${pos.x} ${pos.y}" class="synteny-link"><title>${escapeHtml(syntenyGeneLabel(name))}</title></path>`);
    });
    const block = row.block;
    const rowClass = row.kind === "reference"
      ? selectedSyntenyReference
        ? "synteny-row-bg synteny-reference-bg synteny-row-selected"
        : "synteny-row-bg synteny-reference-bg"
      : row.blockIndex === selectedSyntenyBlock
        ? "synteny-row-bg synteny-row-selected"
        : "synteny-row-bg";
    const rowAttrs = row.kind === "block"
      ? ` class="synteny-row" data-block-index="${row.blockIndex}"`
      : ` class="synteny-reference-row" data-reference-context="1"`;
    parts.push(`<rect x="6" y="${y - 18}" width="${width - 18}" height="${rowH - 8}" rx="8" class="${rowClass}"/>`);
    parts.push(`<g${rowAttrs}>`);
    parts.push(`<rect x="6" y="${y - 18}" width="${width - 18}" height="${rowH - 8}" rx="8" class="synteny-row-hit"/>`);
    if (row.kind === "reference") {
      parts.push(`<text x="14" y="${y + 3}" class="synteny-row-label">H37Rv context</text>`);
      const referenceText = selectedSyntenyReference ? "selected reference tip" : "click to mark reference tip";
      parts.push(`<text x="14" y="${y + 20}" class="synteny-label">${escapeHtml(referenceText)}</text>`);
      parts.push(`<text x="${left}" y="${y - 25}" class="synteny-locus">reference context inferred from flanking clusters</text>`);
    } else {
      const meta = `#${block.rank} | n=${block.nGenomes} | f=${block.fraction || "NA"}`;
      const selectedText = row.blockIndex === selectedSyntenyBlock && !selectedSyntenyReference ? "selected" : "click to mark tips";
      parts.push(`<text x="14" y="${y + 3}" class="synteny-row-label">${escapeHtml(meta)}</text>`);
      parts.push(`<text x="14" y="${y + 20}" class="synteny-label">${escapeHtml(selectedText)}</text>`);
      const repBits = [
        block.representativeGenome || "",
        block.representativeGene || firstExampleLocus(block),
        block.representativeContig && block.representativeStart && block.representativeEnd ? `${block.representativeContig}:${block.representativeStart}-${block.representativeEnd}` : "",
        block.representativeStrand || "",
        block.representativeLength ? `${block.representativeLength} bp` : "",
      ].filter(Boolean).join(" ");
      parts.push(`<text x="${left}" y="${y - 25}" class="synteny-locus">${escapeHtml(clipText(repBits || firstExampleLocus(block), 82))}</text>`);
    }
    const locusY = y - 25;
    if (targetCenter) {
      parts.push(`<line x1="${left}" y1="${locusY + 4}" x2="${targetCenter.x}" y2="${targetCenter.y - geneH / 2}" class="synteny-locus-line"/>`);
    }
    layout.forEach(item => {
      const gene = item.gene;
      const x = item.x;
      const geneW = item.width;
      const fill = gene.gap ? "transparent" : gene.target ? colors.P : row.kind === "reference" ? "transparent" : "transparent";
      const cls = gene.gap ? "synteny-gene synteny-gap" : gene.target ? "synteny-gene synteny-target" : "synteny-gene";
      const labelClass = gene.target ? "synteny-target-label" : "synteny-label";
      const label = syntenyDisplayLabel(gene);
      const coord = syntenyCoordLabel(gene);
      const title = [syntenyGeneLabel(gene.cluster || gene.name), gene.geneId, coord].filter(Boolean).join(" | ");
      const geneAttrs = gene.gap || !gene.cluster || gene.cluster === "NA"
        ? ""
        : ` data-gene-cluster="${escapeHtml(gene.cluster)}"`;
      parts.push(`<path d="${arrowPath(x, y, geneW, geneH, gene.strand)}" fill="${fill}" class="${cls}"${geneAttrs}><title>${escapeHtml(title)}</title></path>`);
      parts.push(`<text x="${x + geneW / 2}" y="${y + 13}" text-anchor="middle" class="${labelClass}">${escapeHtml(clipText(label, 10))}</text>`);
      parts.push(`<text x="${x + geneW / 2}" y="${y + geneH + 13}" text-anchor="middle" class="synteny-label">${escapeHtml(clipText(syntenySecondaryLabel(gene), 12))}</text>`);
      if (gene.length) {
        parts.push(`<text x="${x + geneW / 2}" y="${y + geneH + 28}" text-anchor="middle" class="synteny-label">${escapeHtml(gene.length)} bp</text>`);
      }
    });
    const note = row.kind === "reference" ? "reference row can share neighbours without containing the focal cluster" : (block.interpretation || "").replace(/_/g, " ");
    if (note) {
      parts.push(`<text x="14" y="${y + geneH + 47}" class="synteny-label">${escapeHtml(clipText(note, 45))}</text>`);
    }
    parts.push("</g>");
    previousPositions = positions;
  });
  syntenySvg.setAttribute("width", String(width));
  syntenySvg.setAttribute("height", String(height));
  syntenySvg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  syntenySvg.innerHTML = parts.join("");
  const firstTargetIndex = blockGenes(blocks[0], cluster.name).findIndex(gene => gene.target);
  if (firstTargetIndex >= 0) {
    const targetItem = rowLayouts[0] && rowLayouts[0][firstTargetIndex];
    const targetX = targetItem ? targetItem.x + targetItem.width / 2 : left;
    syntenyPanel.scrollLeft = Math.max(0, targetX - syntenyPanel.clientWidth / 2);
  }
  syntenyPanel.scrollTop = 0;
}
function applyViewerSplit() {
  if (!viewerSplit) return;
  if (syntenyHidden) {
    viewerSplit.style.gridTemplateColumns = "minmax(0,1fr)";
  } else {
    const matrixShare = Math.max(35, Math.min(80, Number(splitControl?.value || 56)));
    const syntenyShare = 100 - matrixShare;
    viewerSplit.style.gridTemplateColumns = `minmax(360px,${matrixShare}fr) 18px minmax(320px,${syntenyShare}fr)`;
  }
  document.querySelectorAll(".synteny-control").forEach(control => {
    control.style.display = syntenyHidden ? "none" : "flex";
  });
}
function applySyntenyVisibility() {
  if (!viewerSplit || !syntenyToggle) return;
  viewerSplit.classList.toggle("synteny-hidden", syntenyHidden);
  syntenyToggle.classList.toggle("active", !syntenyHidden);
  syntenyToggle.textContent = syntenyHidden ? "Show syntenies" : "Hide syntenies";
  syntenyToggle.setAttribute(
    "aria-pressed",
    syntenyHidden ? "false" : "true",
  );
  applyViewerSplit();
  if (!syntenyHidden) {
    syntenyKey = "";
    renderSyntenyPanel();
  } else {
    selectedSyntenyCarriers = new Set();
    treeDirty = true;
  }
}
function renderTree() {
  const xScale = Math.max(24, (treeWidth - 155) / Math.max(1, maxDepth));
  if (treeDirty) {
    const parts = [];
    drawNode(treeRoot, parts, xScale);
    const rootLabel = hasL8Tip(treeRoot) ? "L8-rooted" : "input-rooted";
    if (treeMode === "phylogram" && maxDistance > 0) {
      parts.push(`<text x="14" y="13" class="tree-scale">phylogram, ${rootLabel}, branch length max ${escapeHtml(formatBranchLength(maxDistance))}</text>`);
    } else {
      parts.push(`<text x="14" y="13" class="tree-scale">cladogram, ${rootLabel}, branch lengths ignored</text>`);
    }
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
function renderColumnHeader() {
  const cols = visibleClusters.length;
  const selectedCol = selectedVisibleColumn();
  headerCanvas.width = Math.max(1, cols * cellW);
  headerCanvas.height = headerH;
  headerCanvas.style.width = `${headerCanvas.width}px`;
  headerCanvas.style.height = `${headerH}px`;
  headerCtx.clearRect(0, 0, headerCanvas.width, headerCanvas.height);
  headerCtx.fillStyle = headerBg;
  headerCtx.fillRect(0, 0, headerCanvas.width, headerH);
  if (matrixColorMode.value === "process") {
    visibleClusters.forEach((clusterIdx, c) => {
      headerCtx.fillStyle = processColor(processKind(data.clusters[clusterIdx]));
      headerCtx.fillRect(c * cellW, 0, Math.max(1, cellW), 7);
    });
  }
  if (hoverCol >= 0) {
    headerCtx.fillStyle = hoverColTint;
    headerCtx.fillRect(hoverCol * cellW, 0, cellW, headerH);
  }
  if (selectedCol >= 0) {
    headerCtx.fillStyle = selectColTint;
    headerCtx.fillRect(selectedCol * cellW, 0, cellW, headerH);
    headerCtx.strokeStyle = selectStroke;
    headerCtx.strokeRect(selectedCol * cellW + 0.5, 0.5, Math.max(1, cellW - 1), headerH - 1);
  }
  for (let c = 0; c < cols; c += 1) {
    const clusterIdx = originalClusterAtVisibleColumn(c);
    headerCtx.fillStyle = matrixColorMode.value === "process"
      ? processColor(processKind(data.clusters[clusterIdx]))
      : gridStroke;
    headerCtx.fillRect(c * cellW, headerH - 9, Math.max(1, cellW), 6);
  }
}
function rebuildBaseMatrix() {
  const {rows, cols} = updateCellGeometry();
  baseCanvas.width = Math.max(1, cols * cellW);
  baseCanvas.height = Math.max(1, Math.ceil(rows * cellH));
  baseCtx.clearRect(0, 0, baseCanvas.width, baseCanvas.height);
  baseCtx.fillStyle = panelBg;
  baseCtx.fillRect(0, 0, baseCanvas.width, baseCanvas.height);
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
  canvas.width = baseCanvas.width;
  canvas.height = baseCanvas.height;
  canvas.style.width = `${canvas.width}px`;
  canvas.style.height = `${canvas.height}px`;
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.drawImage(baseCanvas, 0, 0);
  if (selectedCol >= 0) {
    ctx.fillStyle = "rgba(37,99,235,.13)";
    ctx.fillRect(selectedCol * cellW, 0, cellW, canvas.height);
  }
  if (hoverRow >= 0 && hoverRow < rows) {
    ctx.fillStyle = "rgba(245,158,11,.10)";
    ctx.fillRect(0, hoverRow * cellH, canvas.width, Math.max(0.25, cellH));
  }
  if (hoverCol >= 0 && hoverCol < cols) {
    ctx.fillStyle = "rgba(245,158,11,.08)";
    ctx.fillRect(hoverCol * cellW, 0, cellW, canvas.height);
  }
  if (hoverCol >= 0) {
    ctx.strokeStyle = "#b45309";
    ctx.lineWidth = 1;
    ctx.strokeRect(hoverCol * cellW + 0.5, 0.5, Math.max(1, cellW - 1), Math.max(1, canvas.height - 1));
  }
  if (hoverRow >= 0) {
    ctx.strokeStyle = "#b45309";
    ctx.lineWidth = 1;
    ctx.strokeRect(0.5, hoverRow * cellH + 0.5, Math.max(1, canvas.width - 1), Math.max(0.25, cellH - 0.1));
  }
  ctx.strokeStyle = selectStroke;
  ctx.lineWidth = Math.max(1, Math.min(2, cellW));
  if (selectedCol >= 0) {
    ctx.strokeRect(selectedCol * cellW + 0.5, 0.5, Math.max(1, cellW - 1), Math.max(1, canvas.height - 1));
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
    const support = cluster.support ? `<span>support ${escapeHtml(cluster.support)}</span>` : "";
    const phylo = cluster.phylo ? `<span>phylogeny ${escapeHtml(cluster.phylo)}</span>` : "";
    const detail = cluster.mechanismDetail ? `<span>${escapeHtml(cluster.mechanismDetail)}</span>` : "";
    const wording = cluster.wording ? `<div class="process-wording">${escapeHtml(cluster.wording)}</div>` : "";
    info.textContent = `${cluster.present}/${data.genomes.length} genomes (${Math.round(cluster.prevalence * 1000) / 10}%)`;
    details.innerHTML = `<strong>${escapeHtml(cluster.name)}</strong><small>${escapeHtml(cluster.annotation || "No annotation provided")}</small><div class="state-counts"><span class="state-chip">input ${counts.P}</span><span class="state-chip">rescued ${counts.R}</span><span class="state-chip">NSC ${counts.N}</span><span class="state-chip">absent ${counts.A}</span></div><div class="process-note"><span class="process-chip" style="--process-color:${processColor(mechanism)}"><i></i>${escapeHtml(processLabel(mechanism))}</span><span>${escapeHtml(processDescription(mechanism))}</span></div>${referenceInterpretation(cluster)}<div class="process-meta">${support}${phylo}${detail}</div>${wording}`;
    renderProcessLegend();
    detailsKey = nextDetailsKey;
  } else if (nextDetailsKey !== detailsKey) {
    info.textContent = "";
    details.textContent = "No cluster selected.";
    processLegend.innerHTML = "";
    detailsKey = nextDetailsKey;
  }
}
function render(scrollToSelected = false) {
  renderMatrix(scrollToSelected);
  renderColumnHeader();
  renderSyntenyPanel();
  renderTree();
  renderLineageLegend();
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
if (leafOrder.length > 1500) {
  labelToggle.checked = false;
  sizeControl.value = "2";
  rowHeightControl.value = "2";
} else if (leafOrder.length > 800) {
  labelToggle.checked = false;
  sizeControl.value = "3";
  rowHeightControl.value = "3";
} else if (leafOrder.length > 300) {
  labelToggle.checked = false;
  sizeControl.value = "4";
  rowHeightControl.value = "4";
}
populateOptions();
applySyntenyVisibility();
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
    : "Fullscreen explorer";
});
sizeControl.addEventListener("input", () => { resetMatrixCache(); render(false); });
rowHeightControl.addEventListener("input", () => { resetMatrixCache(); render(false); });
splitControl.addEventListener("input", () => { applyViewerSplit(); });
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
    ? adjustRange(sizeControl, delta, 0.25, 18)
    : adjustRange(rowHeightControl, delta, 1, 18);
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
