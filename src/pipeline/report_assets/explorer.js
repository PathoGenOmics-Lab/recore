(() => {
const data = window.recoreExplorerData;
// No data: the boot script has already said why.
if (!data) return;
if (!data.clusters.length) {
  window.recoreExplorerFailed(new Error("the matrix has no clusters to draw"));
  return;
}
{
  const notices = document.getElementById("workspaceNotices");
  if (notices && data.notices && data.notices.length) {
    for (const text of data.notices) {
      const line = document.createElement("p");
      line.textContent = text;
      notices.appendChild(line);
    }
    notices.hidden = false;
  }
}
/* `source` is an optional declaration block to read instead of the live root:
   figures always take the light palette, whatever the screen is in. */
function readThemeVar(name, fallback, source) {
  const value = (source || getComputedStyle(document.documentElement)).getPropertyValue(name).trim();
  return value || fallback;
}
/** The light theme's tokens, read from the stylesheet rather than by flipping
    data-theme, which would restyle and transition the whole page. */
function lightThemeDeclaration() {
  for (const sheet of document.styleSheets) {
    let rules = null;
    try {
      rules = sheet.cssRules;
    } catch (_e) {
      continue;
    }
    for (const rule of rules || []) {
      if (rule.selectorText && rule.selectorText.includes(':root[data-theme="light"]')) return rule.style;
    }
  }
  return null;
}
function readVar(name, fallback) {
  return readThemeVar(name, fallback);
}
function refreshThemeColors(source) {
  const readVar = (name, fallback) => readThemeVar(name, fallback, source);
  colors.P = readVar("--present", "#b3bbc7");
  colors.R = readVar("--rescue", "#0f6fd6");
  colors.N = readVar("--nsc", "#e07b00");
  colors.K = readVar("--kept", "#c01048");
  colors.A = readVar("--absent", "#f4f5f7");
  strokes.P = colors.P;
  strokes.R = colors.R;
  strokes.N = colors.N;
  strokes.K = colors.K;
  strokes.A = readVar("--absent-stroke", "#d0d5dd");
  for (const family of Object.keys(mechanismFamilyColors)) {
    mechanismFamilyColors[family] = readVar(`--mech-${family}`, mechanismFamilyColors[family]);
  }
  burdenLow = readVar("--line", "#e4e4e7");
  qcColors.pass = readVar("--faint", "#a1a1aa");
  qcColors.review = readVar("--warning", "#d97706");
  qcColors.fail = readVar("--danger", "#dc2626");
  panelBg = readVar("--panel", "#ffffff");
  headerBg = readVar("--bg-subtle", "#f8fafc");
  hoverColTint = readVar("--warning-soft", "#fff7ed");
  hoverStroke = readVar("--muted", "#667085");
  markInk = readVar("--ink", "#101828");
  selectColTint = readVar("--accent-soft", "#eef4ff");
  selectStroke = readVar("--accent-ink", "#1d4ed8");
  gridStroke = readVar("--faint", "#98a2b3");
  labelInk = readVar("--subtle", "#475467");
}
const clusterByName = new Map(data.clusters.map(cluster => [cluster.name, cluster]));
const colors = {P:"#b3bbc7", R:"#0f6fd6", N:"#e07b00", K:"#c01048", A:"#f4f5f7"};
const strokes = {P:"#b3bbc7", R:"#0f6fd6", N:"#e07b00", K:"#c01048", A:"#d0d5dd"};
let burdenLow = "#e4e4e7";
const qcColors = {pass: "#a1a1aa", review: "#d97706", fail: "#dc2626"};
let panelBg = "#ffffff";
let headerBg = "#f8fafc";
let hoverColTint = "#fff7ed";
let hoverStroke = "#667085";
let markInk = "#101828";
let selectColTint = "#eef4ff";
let selectStroke = "#1d4ed8";
let gridStroke = "#98a2b3";
let labelInk = "#475467";
const labels = {P:"present in input", R:"corrected by Recore", N:"supported by NSC", K:"found, kept out of matrix", A:"absent"};
const shortLabels = {P:"input", R:"corrected", N:"NSC", K:"kept out", A:"absent"};
const stateDefinitions = {
  P:"Present in the input matrix.",
  R:"Absent in the input. Recore found the gene in this genome and corrected the cell.",
  N:"Absent in the input. Near-identical clusters were merged (NSC), and the merged cluster is present here.",
  K:"Recore found sequence evidence, but the matrix policy, a quality floor or a synteny or paralog guard kept it out of the corrected matrix. Click the cell for the reason.",
  A:"Absent, and Recore found no evidence for it.",
};
/* Sixteen mechanisms cannot get sixteen colours anyone can tell apart, so they
   are drawn in six families. The family hues stay off the cell-state hues and
   apart from each other under deutan and protan vision; the inspector and the
   legend still name the exact mechanism. */
const mechanismFamilies = {
  ancestral_retained_no_disruption:"retained",
  clade_associated_gain:"gain",
  clean_insertion_or_acquisition_candidate:"gain",
  deletion_loss_candidate:"loss",
  reference_or_lineage_loss_candidate:"loss",
  IS_interruption_supported:"disruption",
  pseudogenization:"disruption",
  gene_fragmentation_or_structural_break:"disruption",
  reference_gene_disruption_or_truncation:"disruption",
  duplication_or_paralog_expansion:"paralog",
  translocation_or_context_shift_candidate:"paralog",
  partial_reference_homology:"paralog",
  reference_gene_allelic_derivative:"paralog",
  annotation_or_clustering_artifact:"unresolved",
  unstable_or_homoplastic_accessory_signal:"unresolved",
  unresolved_accessory_variation:"unresolved",
};
const mechanismFamilyLabels = {
  retained:"retained",
  gain:"gain or insertion",
  loss:"loss or deletion",
  disruption:"disruption",
  paralog:"paralog or context change",
  unresolved:"artifact or unresolved",
};
const mechanismFamilyColors = {
  retained:"#667085",
  gain:"#4ca30d",
  loss:"#7a2e0e",
  disruption:"#fac515",
  paralog:"#15b79e",
  unresolved:"#1d2939",
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
  A1:"Animal-associated",
  A2:"Animal-associated",
  A3:"Animal-associated",
  A4:"Animal-associated",
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
const cellLegend = document.getElementById("cellLegend");
const status = document.getElementById("matrixStatus");
const headerCanvas = document.getElementById("columnHeaderCanvas");
const headerCtx = headerCanvas.getContext("2d");
const canvas = document.getElementById("matrixCanvas");
const matrixViewport = canvas ? canvas.parentElement : null;
const matrixSpacer = document.getElementById("matrixSpacer");
const headerSpacer = document.getElementById("headerSpacer");
const treeCanvas = document.getElementById("treeCanvas");
// Swapped for a figure's context while a figure is drawn.
let treeCtx = treeCanvas ? treeCanvas.getContext("2d") : null;
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
// Hover and selection marks live on transparent canvases above the data, so
// moving the pointer repaints a few lines instead of every pane, and the data
// canvases (and the exported figure) never carry a crosshair.
const matrixOverlay = document.getElementById("matrixOverlay");
const matrixOverlayCtx = matrixOverlay ? matrixOverlay.getContext("2d") : null;
const headerOverlay = document.getElementById("headerOverlay");
const headerOverlayCtx = headerOverlay ? headerOverlay.getContext("2d") : null;
const treeOverlay = document.getElementById("treeOverlay");
const treeOverlayCtx = treeOverlay ? treeOverlay.getContext("2d") : null;
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
const MIN_COL_ZOOM = 0.01;
const MAX_COL_ZOOM = 90;
const MIN_ROW_ZOOM = 0.02;
const MAX_ROW_ZOOM = 40;
// Target pixel span used to pick the opening zoom from the run size.
const FIT_TARGET_H = 560;
const FIT_TARGET_W = 720;
// Column names are only worth drawing once a cell is wide enough to carry a
// rotated label without the text of neighbouring columns touching.
const LABEL_MIN_CELL_W = 9;
// Tip names are 10px text: below 11px rows they overlap their neighbours.
const TIP_LABEL_MIN_ROW = 11;
const LABEL_BAND_H = 82;
// Above this a name fits across the column, and turning your head to read a
// dozen columns is a cost the layout no longer has to impose.
const LABEL_HORIZONTAL_CELL_W = 54;
const LABEL_BAND_H_FLAT = 26;
const LABEL_FONT = "10px Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif";
const defaultStatus = "Hover a cell for its evidence, click it for the full row. Ctrl or Cmd and the wheel zoom columns, Alt and the wheel zoom rows, Z fits everything. Drag the seams between panes to resize.";
// Below this width the inspector is a bottom sheet (see report.css).
const narrowQuery = window.matchMedia ? window.matchMedia("(max-width: 1100px)") : null;
let selected = 0;
let selectedGenome = -1;
// True while an exported figure is being drawn (see drawFigure).
let drawingFigure = false;
let hoverRow = -1;
let hoverCol = -1;
let leafOrder = [];
/* Matrix rows follow tree-tip order while cell sets are keyed by matrix genome
   index, so drawing needs the inverse mapping. -1 marks a genome with no tip. */
let rowOfGenome = new Int32Array(0);
let visibleClusters = [];
let genomeStats = [];
let maxGenomeStats = {R:1, N:1, A:1, correction:1};
let anyKeptOut = false;
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
const MIN_CELL_PX = 0.005;
// Below these a cell is less than a pixel, and the matrix is drawn as density:
// each device pixel shows the share of present cells it covers, with any
// corrected or kept-out cell in it marked. Drawing cells one rect at a time
// there cost millions of rects and painted neighbours over each other.
const DENSITY_BELOW_COL = 0.35;
const DENSITY_BELOW_ROW = 1;
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
/* Set once the genome card exists (it is declared further down); render() runs
   before that during start-up. */
let afterRender = null;
let treeRoot = data.tree;
let leafByGenome = new Map();
let offTreeGenomes = [];

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
// The cell states, in the order the writer stores their sets. K is evidence
// the run found that did not reach the binary matrix.
const STATE_CODES = ["P", "R", "N", "K", "A"];
const STORED_STATES = 4;
const EMPTY_SET = {kind: "list", idx: new Int32Array(0)};
const decodedRows = new Map();

function decodeBase64(text) {
  const binary = atob(text);
  const bytes = new Uint8Array(binary.length);
  for (let i = 0; i < binary.length; i += 1) bytes[i] = binary.charCodeAt(i);
  return bytes;
}
/* ---- Per-cell evidence ----------------------------------------------------
   The rows of orthology_evidence.tsv, one base64 string per cluster, decoded
   only for the clusters someone looks at. */
const evidenceData = data.evidence && data.evidence.byCluster ? data.evidence : null;
const evidenceCache = new Map();
function clusterEvidence(clusterIdx) {
  if (!evidenceData) return [];
  const cached = evidenceCache.get(clusterIdx);
  if (cached) return cached;
  const rows = decodeEvidence(clusterIdx);
  if (evidenceCache.size >= 64) evidenceCache.delete(evidenceCache.keys().next().value);
  evidenceCache.set(clusterIdx, rows);
  return rows;
}
/** Number of evidence rows of a cluster, read from the first bytes only. */
function evidenceCount(clusterIdx) {
  const encoded = evidenceData ? evidenceData.byCluster[clusterIdx] || "" : "";
  if (!encoded) return 0;
  const bytes = decodeBase64(encoded.slice(0, 8));
  let value = 0;
  let scale = 1;
  for (const byte of bytes) {
    value += (byte & 0x7f) * scale;
    if (!(byte & 0x80)) break;
    scale *= 128;
  }
  return value;
}
function decodeEvidence(clusterIdx) {
  const encoded = evidenceData.byCluster[clusterIdx] || "";
  const rows = [];
  if (encoded) {
    const bytes = decodeBase64(encoded);
    let pos = 0;
    const next = () => {
      let value = 0;
      let scale = 1;
      for (;;) {
        const byte = bytes[pos++];
        value += (byte & 0x7f) * scale;
        if (!(byte & 0x80)) return value;
        scale *= 128;
      }
    };
    const strings = evidenceData.strings;
    const count = next();
    let genome = 0;
    for (let i = 0; i < count; i += 1) {
      genome += next();
      const profile = evidenceData.profiles[next()] || [];
      const row = {genome};
      evidenceData.profileColumns.forEach((column, j) => { row[column] = strings[profile[j]] || ""; });
      evidenceData.textColumns.forEach(column => { row[column] = strings[next()] || ""; });
      const start = next();
      const end = next();
      const identity = next();
      const coverage = next();
      row.start = start ? start - 1 : null;
      row.end = end ? end - 1 : null;
      row.identity = identity ? (identity - 1) / 10 : null;
      row.coverage = coverage ? (coverage - 1) / 10 : null;
      rows.push(row);
    }
  }
  return rows;
}
/** Evidence rows for one cell; rows are sorted by genome. */
function cellEvidence(clusterIdx, genomeIdx) {
  const rows = clusterEvidence(clusterIdx);
  let lo = 0;
  let hi = rows.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (rows[mid].genome < genomeIdx) lo = mid + 1;
    else hi = mid;
  }
  const out = [];
  while (lo < rows.length && rows[lo].genome === genomeIdx) out.push(rows[lo++]);
  return out;
}
const MATRIX_ACTION_TEXT = {
  correct_presence:"used in the corrected matrix",
  audit_only_policy:"kept out by the matrix policy",
  audit_only_quality:"kept out by a quality floor",
  audit_only:"not used: rejected",
};
const PUBLICATION_TEXT = {
  publication_ready_correction:"ready for publication",
  corrected_but_review_for_publication:"corrected; review before publishing",
  corrected_paralog_review:"corrected; paralog review",
  corrected_review:"corrected; review",
  held_out_policy:"held out by the matrix policy",
  held_out_quality:"held out by a quality floor",
  audit_only_not_corrected:"audit only, not corrected",
};
function evidenceTypeLabel(row) {
  return (row.rescue_type || row.evidence_state || "evidence").split(/[:(]/)[0].trim();
}
function pctText(value) {
  return value === null || value === undefined ? "n/a" : `${value.toFixed(1)}%`;
}
function evidenceLocus(row) {
  if (!row.contig) return "";
  const span = row.start !== null && row.end !== null ? `:${formatCoord(row.start)}-${formatCoord(row.end)}` : "";
  return `${row.contig}${span}${row.strand ? ` (${row.strand})` : ""}`;
}
/** One line for the tooltip and the changed-cells list. */
function evidenceSummary(row) {
  return `${evidenceTypeLabel(row)}, ${pctText(row.identity)} identity, ${pctText(row.coverage)} coverage`;
}
function evidenceRowHtml(row) {
  const accepted = row.decision === "accepted";
  const why = row.rejection_reason ? `; ${row.rejection_reason.replace(/_/g, " ")}` : "";
  const facts = [];
  if (row.synteny) facts.push(`synteny ${escapeHtml(row.synteny)}${row.synteny_score ? ` (${escapeHtml(row.synteny_score)})` : ""}`);
  if (row.orf_integrity_signal) facts.push(`ORF ${escapeHtml(row.orf_integrity_signal.replace(/_?orf_?/gi, " ").replace(/_/g, " ").trim())}`);
  if (row.orthology_confidence_tier) facts.push(`${escapeHtml(row.orthology_confidence_tier.replace(/_confidence$/, "").replace(/_/g, " "))} confidence`);
  const notes = [];
  if (row.risk_flags) notes.push(`<dt>Flags</dt><dd>${escapeHtml(row.risk_flags.replace(/;/g, "; ").replace(/_/g, " "))}</dd>`);
  if (row.fragment_note) notes.push(`<dt>Note</dt><dd>${escapeHtml(row.fragment_note.replace(/_/g, " "))}</dd>`);
  if (row.annotation_match) notes.push(`<dt>Annotated as</dt><dd>${escapeHtml(row.annotation_match)}${row.annotation_cluster ? ` (${escapeHtml(row.annotation_cluster)})` : ""}</dd>`);
  if (row.recommended_wording) notes.push(`<dt>Wording</dt><dd>${escapeHtml(row.recommended_wording)}</dd>`);
  return `<div class="evidence-row${accepted ? "" : " rejected"}">` +
    `<div class="ev-head"><span class="ev-type">${escapeHtml(row.rescue_type || row.evidence_state || "evidence")}</span>` +
    `<span class="ev-action">${escapeHtml(MATRIX_ACTION_TEXT[row.matrix_action] || row.matrix_action || row.decision)}${escapeHtml(why)}</span></div>` +
    (evidenceLocus(row) ? `<div class="ev-locus mono">${escapeHtml(evidenceLocus(row))}</div>` : "") +
    `<div class="ev-metrics num">${pctText(row.identity)} identity &middot; ${pctText(row.coverage)} coverage</div>` +
    (facts.length ? `<div class="ev-facts">${facts.join(" &middot; ")}</div>` : "") +
    (row.publication_status ? `<div class="ev-pub">Publication: ${escapeHtml(PUBLICATION_TEXT[row.publication_status] || row.publication_status.replace(/_/g, " "))}</div>` : "") +
    (notes.length ? `<dl class="cluster-dl ev-notes">${notes.join("")}</dl>` : "") +
    `</div>`;
}
function emptyCellText(state) {
  if (state === "P") return "Present in the input matrix. Recore had nothing to correct here.";
  if (state === "N") return "Merged into a near-identical cluster by NSC; no search evidence row.";
  return "No evidence row. Either Recore searched this gap and found nothing, or the gap was outside the search (see the prevalence filter in Run details).";
}
/** The selected cell: what it is, and every evidence row behind it. */
function cellCardHtml() {
  if (selectedGenome < 0 || !data.clusters[selected]) return "";
  const state = stateAt(selected, selectedGenome);
  const rows = cellEvidence(selected, selectedGenome);
  const body = rows.length
    ? rows.map(evidenceRowHtml).join("")
    : `<p class="muted cell-empty">${escapeHtml(emptyCellText(state))}</p>`;
  return `<section class="cell-card" aria-label="Selected cell">` +
    `<div class="cell-card-head"><span class="card-title">Cell</span>` +
    `<button type="button" class="icon-btn ghost-icon cell-close" data-clear-cell="1" aria-label="Clear the selected cell" data-tip="Clear the selected cell">&times;</button></div>` +
    `<p class="cell-name"><button type="button" class="link-btn" data-find-tip="${selectedGenome}">${escapeHtml(data.genomes[selectedGenome] || "")}</button></p>` +
    `<p class="cell-state"><i class="sw ${state.toLowerCase()}"></i>${escapeHtml(labels[state])}</p>` +
    body + `</section>`;
}
/* Which lineages carry this cluster's corrections. Counting walks only the
   corrected and kept-out sets, so it costs the size of the finding, not of the
   run; it was 14 screens of scrolling with tips coloured by state. */
let genomeLineages = null;
let lineageSizes = null;
function lineageBreakdownHtml(clusterIdx) {
  if (!genomeLineages) {
    genomeLineages = data.genomes.map(name => lineageOfGenomeName(name));
    lineageSizes = new Map();
    for (const lineage of genomeLineages) lineageSizes.set(lineage, (lineageSizes.get(lineage) || 0) + 1);
  }
  if (lineageSizes.size < 2 && lineageSizes.has("unknown")) return "";
  const sets = clusterSets(clusterIdx);
  const n = data.genomes.length;
  const tally = new Map();
  const add = (set, key) => {
    if (set.kind === "list") {
      for (const g of set.idx) {
        const row = tally.get(genomeLineages[g]) || {R: 0, K: 0};
        row[key] += 1;
        tally.set(genomeLineages[g], row);
      }
    } else if (set.kind === "all" || set.kind === "complement") {
      const skip = set.kind === "complement" ? new Set(set.idx) : null;
      for (let g = 0; g < n; g += 1) {
        if (skip && skip.has(g)) continue;
        const row = tally.get(genomeLineages[g]) || {R: 0, K: 0};
        row[key] += 1;
        tally.set(genomeLineages[g], row);
      }
    }
  };
  add(sets[1], "R");
  add(sets[3], "K");
  if (!tally.size) return "";
  const rows = [...tally.entries()].sort((a, b) => (b[1].R + b[1].K) - (a[1].R + a[1].K)).slice(0, 8);
  const body = rows.map(([lineage, c]) => {
    const size = lineageSizes.get(lineage) || 1;
    const r = c.R / size * 100;
    const k = c.K / size * 100;
    const words = [c.R ? `${formatInt(c.R)} corrected` : "", c.K ? `${formatInt(c.K)} kept out` : ""].filter(Boolean).join(", ");
    return `<div class="lineage-row" data-tip="${escapeHtml(`${lineageLabel(lineage)}: ${words} of ${plural(size, "genome")}`)}">` +
      `<span class="lineage-name"><i style="background:${lineageColor(lineage)}"></i>${escapeHtml(lineage)}</span>` +
      `<span class="lineage-bar"><b class="r" style="width:${r.toFixed(1)}%"></b><b class="k" style="width:${k.toFixed(1)}%"></b></span>` +
      `<small class="num">${escapeHtml(words)} <span class="muted">/ ${formatInt(size)}</span></small></div>`;
  }).join("");
  return `<div class="card-title">Corrections by lineage</div><div class="lineage-rows">${body}</div>`;
}
/** The cells Recore changed in this cluster, findings first. */
const CHANGED_CELLS_SHOWN = 20;
function changedCellsHtml(clusterIdx) {
  const rows = clusterEvidence(clusterIdx);
  if (!rows.length) return "";
  const rank = row => (row.matrix_action === "correct_presence" ? 0 : row.decision === "accepted" ? 1 : 2);
  const sorted = rows.slice().sort((a, b) => rank(a) - rank(b) || (rowOfGenome[a.genome] - rowOfGenome[b.genome]));
  const items = sorted.slice(0, CHANGED_CELLS_SHOWN).map(row => {
    const state = stateAt(clusterIdx, row.genome);
    return `<li><button type="button" class="changed-cell" data-cell-genome="${row.genome}">` +
      `<i class="sw ${state.toLowerCase()}"></i><span class="changed-genome">${escapeHtml(data.genomes[row.genome] || "")}</span>` +
      `<small>${escapeHtml(evidenceSummary(row))}</small></button></li>`;
  }).join("");
  const more = sorted.length > CHANGED_CELLS_SHOWN
    ? `<p class="muted changed-more">and ${plural(sorted.length - CHANGED_CELLS_SHOWN, "more row")} in orthology_evidence.tsv</p>`
    : "";
  return `<div class="card-title">Evidence rows (${formatInt(rows.length)})</div><ul class="changed-cells">${items}</ul>${more}`;
}
function decodeClusterSets(encoded, nGenomes, count = STORED_STATES) {
  const sets = Array.from({length: count}, () => EMPTY_SET);
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
  for (let state = 0; state < count; state += 1) {
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
/* Drawing wants a set as matrix rows, sorted, so the visible band is a binary
   search plus a walk over what is on screen. Walking genome indices instead
   costs the whole set per column: a 50%-prevalence cluster at 10,000 genomes is
   5,000 checks to paint the 700 rows in view. Rebuilt when the tree order
   changes (rerooting moves every row). */
let rowEpoch = 0;
const rowSets = new Map();
function clusterRowSets(clusterIdx) {
  const cached = rowSets.get(clusterIdx);
  if (cached && cached.epoch === rowEpoch) return cached.rows;
  const rows = clusterSets(clusterIdx).map(set => {
    if (set.kind === "all") return null;
    const out = new Int32Array(set.idx.length);
    let n = 0;
    for (let i = 0; i < set.idx.length; i += 1) {
      const row = rowOfGenome[set.idx[i]];
      if (row >= 0) out[n++] = row;
    }
    return out.subarray(0, n).sort();
  });
  rowSets.set(clusterIdx, {epoch: rowEpoch, rows});
  return rows;
}
function buildRowOfGenome() {
  rowEpoch += 1;
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
  for (let state = STORED_STATES - 1; state >= 0; state -= 1) {
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
    K: setSize(sets[3], n),
    A: 0,
  };
  counts.A = Math.max(0, n - counts.P - counts.R - counts.N - counts.K);
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
    const totals = (data.genomeStats && data.genomeStats[gi]) || [0, 0, 0, 0, 0];
    return {P: totals[0], R: totals[1], N: totals[2], K: totals[3], A: totals[4]};
  });
  anyKeptOut = genomeStats.some(s => s.K > 0);
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
function mechanismFamily(mechanism) {
  return mechanismFamilies[mechanism] || "unresolved";
}
/** Canvas colour of a mechanism, read from the theme tokens. */
function processColor(mechanism) {
  return mechanismFamilyColors[mechanismFamily(mechanism)];
}
/** The same colour for markup, so it follows a theme change without a rebuild. */
function processCss(mechanism) {
  return `var(--mech-${mechanismFamily(mechanism)})`;
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
/* The cell legend says what the matrix is drawing now: the states this run
   has in evidence mode, the mechanism families in view in mechanism mode. It
   was static markup, so in mechanism mode green read as "rescued". */
let cellLegendKey = "";
function renderCellLegend() {
  if (!cellLegend) return;
  const byProcess = matrixColorMode.value === "process";
  const key = byProcess ? `m|${visibleClusters.length}|${mechanismFilterKey()}` : `s|${anyKeptOut}|${data.nscEnabled}`;
  if (key === cellLegendKey) return;
  cellLegendKey = key;
  const item = (cls, style, label, tip) =>
    `<span class="key" data-tip="${escapeHtml(tip)}"><i class="sw ${cls}"${style ? ` style="${style}"` : ""}></i>${escapeHtml(label)}</span>`;
  if (!byProcess) {
    cellLegend.innerHTML = shownStates()
      .map(state => item(state.toLowerCase(), "", shortLabels[state], stateDefinitions[state]))
      .join("");
    return;
  }
  const families = new Map();
  for (const idx of visibleClusters) {
    const family = mechanismFamily(processKind(data.clusters[idx]));
    families.set(family, (families.get(family) || 0) + 1);
  }
  // Six families do not fit a pane header, so they share one key that opens
  // the full, filterable list in the matrix menu.
  const present = Object.keys(mechanismFamilyLabels).filter(family => families.has(family));
  const tip = present
    .map(family => `${mechanismFamilyLabels[family]}: ${plural(families.get(family), "cluster")}`)
    .join("; ");
  const stack = present.map(family => `<i class="sw" style="background:var(--mech-${family})"></i>`).join("");
  const keys = [`<button type="button" class="key key-btn" data-open-mechanisms="1" data-tip="${escapeHtml(`Input cells take the colour of their cluster's mechanism family. ${tip}. Click for the list and filters.`)}"><span class="sw-stack">${stack}</span>mechanism</button>`];
  keys.push(item("r", "", shortLabels.R, stateDefinitions.R));
  if (anyKeptOut) keys.push(item("k", "", shortLabels.K, stateDefinitions.K));
  cellLegend.innerHTML = keys.join("");
}
function renderLineageLegend() {
  if (!lineageLegend) return;
  const mode = colorMode.value;
  if (mode === "rescue" || mode === "absence" || mode === "correction") {
    // Burden modes are a ramp: say where it starts and what its top means.
    const [, max, high, noun] = burdenOf(mode, {R: 0, N: 0, A: 0});
    lineageLegend.style.display = "flex";
    lineageLegend.innerHTML = `<span class="lineage-title">Tips coloured by ${escapeHtml(noun)} per genome</span>` +
      `<span class="lineage-chip"><i class="ramp-zero"></i><small>0</small></span>` +
      `<span class="burden-ramp" style="background:linear-gradient(90deg, ${burdenLow}, ${high})"></span>` +
      `<span class="lineage-chip"><small>${formatInt(max)}</small></span>`;
    return;
  }
  if (mode !== "lineage") {
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
    .map(lineage => `<span class="lineage-chip"><i style="background:${lineageColor(lineage)}"></i><strong>${escapeHtml(lineage)}</strong><small>${plural(counts.get(lineage), "tip")}</small></span>`)
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
  // A compact checklist: swatch, name and count, the description in the tip.
  // The clear control sits first, so it is never scrolled out of view.
  const head = mechanismFilters.size
    ? `<div class="process-head"><span>${escapeHtml(mechanismFilterSummary())}: ${plural(visibleClusters.length, "cluster")}</span><button type="button" class="process-clear" data-clear-mechanism="1">Clear</button><button type="button" class="process-clear" data-download-visible="1">Download TSV</button></div>`
    : `<div class="process-head"><span>Tick mechanisms to show only their clusters.</span><button type="button" class="process-clear" data-download-visible="1">Download TSV</button></div>`;
  const items = visibleProcessCounts().map(([mechanism, row]) => {
    const active = mechanism === selectedMechanism ? " active" : "";
    const on = mechanismFilters.has(mechanism);
    return `<button type="button" class="process-item${active}${on ? " filtered" : ""}" role="menuitemcheckbox" aria-checked="${on}" data-mechanism="${escapeHtml(mechanism)}" data-tip="${escapeHtml(processDescription(mechanism))}"><span class="process-check" aria-hidden="true">${on ? "&#10003;" : ""}</span><span class="process-swatch" style="background:${processCss(mechanism)}"></span><span class="process-name">${escapeHtml(processLabel(mechanism))}</span><small class="num">${formatInt(row.clusters)}</small></button>`;
  }).join("");
  processLegend.innerHTML = `${head}<div class="process-list">${items}</div>`;
}
function tipStyle(genomeIdx) {
  const mode = colorMode.value;
  // A figure has no selected cluster, so tips coloured by one would mean
  // something the figure cannot show.
  if (drawingFigure && (mode === "cluster" || mode === "process")) {
    return {fill: labelInk, stroke: labelInk, title: ""};
  }
  if (mode === "cluster") {
    const state = stateAt(selected, genomeIdx);
    return {fill: colors[state], stroke: strokes[state], title: labels[state]};
  }
  if (mode === "process") {
    const cluster = data.clusters[selected] || {};
    const state = stateAt(selected, genomeIdx);
    const mechanism = processKind(cluster);
    const color = processColor(mechanism);
    if (state === "A") return {fill: colors.A, stroke: strokes.A, title: `${processLabel(mechanism)}; absent`};
    if (state === "R" || state === "K") return {fill: colors[state], stroke: strokes[state], title: `${processLabel(mechanism)}; ${labels[state]}`};
    return {fill: color, stroke: color, title: `${processLabel(mechanism)}; ${labels[state]}`};
  }
  if (mode === "qc") {
    const qc = genomeQc(genomeIdx);
    const status = qc ? qc.status : "";
    const fill = status === "fail" ? qcColors.fail : status === "review" ? qcColors.review : status === "pass" ? qcColors.pass : colors.A;
    return {fill, stroke: status ? fill : strokes.A, title: status ? `QC ${status}${qc.flags ? `: ${qc.flags.replace(/;/g, "; ")}` : ""}` : "no QC row"};
  }
  if (mode === "lineage") {
    const genome = data.genomes[genomeIdx] || "";
    const lineage = lineageOfGenomeName(genome);
    const color = lineageColor(lineage);
    return {fill: color, stroke: color, title: lineageLabel(lineage)};
  }
  const stats = genomeStats[genomeIdx] || {R:0, N:0, A:0};
  const [value, max, high, noun] = burdenOf(mode, stats);
  // Zero is a hollow ring, so "none" never passes for a pale "some" in either
  // theme, and the ramp runs from the line colour to the state's own colour.
  if (!value) return {fill: panelBg, stroke: gridStroke, title: `${noun}: 0`};
  return {fill: gradientColor(value, max, burdenLow, high), stroke: high, title: `${noun}: ${formatInt(value)}`};
}
/** The value a burden tip mode shows, its maximum, its colour and its name. */
function burdenOf(mode, stats) {
  if (mode === "rescue") return [stats.R, maxGenomeStats.R, colors.R, "corrected cells"];
  if (mode === "absence") return [stats.A, maxGenomeStats.A, labelInk, "absent cells"];
  return [stats.R + stats.N, maxGenomeStats.correction, colors.R, "corrected or NSC cells"];
}
function subtreeLineage(node) {
  return node._lineage || "unknown";
}
/* Synteny sits in its own script element and is parsed the first time the panel
   needs it; until then only the per-cluster block counts are in memory. */
let syntenyStore = null;
function syntenyBlocksFor(cluster) {
  if (!cluster) return [];
  if (!syntenyStore) {
    const element = document.getElementById("syntenyData");
    try {
      syntenyStore = element ? JSON.parse(element.textContent) : {};
    } catch (error) {
      syntenyStore = {};
      status.textContent = `Synteny data could not be read: ${error.message}`;
    }
    if (element) element.remove();
  }
  return syntenyStore[cluster.name] || [];
}
function syntenyBlockCount(clusterIdx) {
  return data.syntenyCounts ? (data.syntenyCounts[clusterIdx] || 0) : 0;
}
function selectedSyntenyBlockData() {
  const cluster = data.clusters[selected];
  const blocks = syntenyBlocksFor(cluster);
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
  if (block._carriers) return block._carriers;
  const n = data.genomes.length;
  const set = decodeClusterSets(block.carriers || "", n, 1)[0];
  let out;
  if (set.kind === "all") out = Array.from({length: n}, (_, g) => g);
  else if (set.kind === "complement") {
    const skip = new Set(set.idx);
    out = [];
    for (let g = 0; g < n; g += 1) if (!skip.has(g)) out.push(g);
  } else out = Array.from(set.idx || []);
  block._carriers = out;
  return out;
}
/* ---- Tree walks ------------------------------------------------------------
   Every walk uses an explicit stack. A ladder-shaped tree of 10,000 tips is
   10,000 levels deep, and the recursive versions threw RangeError past about
   2,500 levels, which left the whole workspace blank. */
function copyNode(node) {
  return {label: node.label || "", genome: node.genome, length: Number(node.length || 0), children: []};
}
function cloneNode(node) {
  const root = copyNode(node);
  const stack = [[node, root]];
  while (stack.length) {
    const [source, target] = stack.pop();
    for (const child of source.children || []) {
      const copy = copyNode(child);
      target.children.push(copy);
      stack.push([child, copy]);
    }
  }
  return root;
}
/** Nodes with every parent before its children, leaves left to right. */
function preorder(root) {
  const out = [];
  const stack = [root];
  while (stack.length) {
    const node = stack.pop();
    out.push(node);
    const children = node.children || [];
    for (let i = children.length - 1; i >= 0; i -= 1) stack.push(children[i]);
  }
  return out;
}
function isLeaf(node) {
  return !node.children || !node.children.length;
}
function hasGenome(node) {
  return node.genome !== null && node.genome !== undefined;
}
function tipName(node) {
  return !hasGenome(node)
    ? (node.label || "")
    : (data.genomes[node.genome] || node.label || "");
}
function isL8Tip(node) {
  return isLeaf(node) && /^L8[_-]/.test(tipName(node));
}
function assignTreeIds(root, nodes, edges) {
  root._parent = null;
  for (const node of preorder(root)) {
    node._id = nodes.length;
    nodes.push(node);
    for (const child of node.children || []) child._parent = node;
  }
  for (const node of nodes) {
    for (const child of node.children || []) edges.push([node._id, child._id, Number(child.length || 0)]);
  }
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
function treeHasTip(root, test) {
  return preorder(root).some(node => isLeaf(node) && test(node));
}
/* Drops leaves that carry no genome, which rerooting turned an unmatched or
   unary branch into: such a tip drew a nameless row painted present in every
   core column. Collapses one-child nodes into their child, summing the branch
   lengths. Counts tips and marks L8 subtrees once per node, so sorting the
   children costs nothing more. */
function tidyTree(root, sortBy) {
  const order = preorder(root);
  for (let i = order.length - 1; i >= 0; i -= 1) {
    const node = order[i];
    const kept = [];
    for (let child of node.children || []) {
      if (isLeaf(child) && !hasGenome(child)) continue;
      while (child.children.length === 1) {
        const only = child.children[0];
        only.length = Number(only.length || 0) + Number(child.length || 0);
        child = only;
      }
      kept.push(child);
    }
    node.children = kept;
    node._tips = kept.length ? kept.reduce((sum, child) => sum + child._tips, 0) : (hasGenome(node) ? 1 : 0);
    node._l8 = kept.length ? kept.some(child => child._l8) : isL8Tip(node);
    if (sortBy === "l8") kept.sort((x, y) => Number(y._l8) - Number(x._l8) || y._tips - x._tips);
    else if (sortBy === "tips") kept.sort((x, y) => y._tips - x._tips);
  }
  let top = root;
  while (top.children.length === 1) top = top.children[0];
  top.length = 0;
  return top;
}
function rootByL8(node) {
  const root = cloneNode(node);
  const nodes = [];
  assignTreeIds(root, nodes, []);
  const leaves = nodes.filter(isL8Tip);
  if (!leaves.length) return root;
  const rootNode = leaves.length === 1 && leaves[0]._parent ? leaves[0]._parent : commonAncestor(leaves);
  const {adj} = graphOf(nodes);
  return buildFrom(nodes, adj, rootNode ? rootNode._id : root._id, null, 0);
}
function findLeafByGenome(_node, genomeIdx) {
  return leafByGenome.get(genomeIdx) || null;
}
/* Depth and distance on the way down; rows, row spans and the lineage of each
   subtree on the way up. One pass each, whatever the shape of the tree. */
function layoutTree(root) {
  leafOrder = [];
  leafByGenome = new Map();
  maxDepth = 0;
  maxDistance = 0;
  const order = preorder(root);
  root.depth = 0;
  root.distance = 0;
  for (const node of order) {
    if (node.depth > maxDepth) maxDepth = node.depth;
    if (node.distance > maxDistance) maxDistance = node.distance;
    if (isLeaf(node)) {
      node.row = leafOrder.length;
      leafOrder.push(node.genome);
      leafByGenome.set(node.genome, node);
      continue;
    }
    for (const child of node.children) {
      child.depth = node.depth + 1;
      child.distance = node.distance + Number(child.length || 0);
    }
  }
  for (let i = order.length - 1; i >= 0; i -= 1) {
    const node = order[i];
    if (isLeaf(node)) {
      node._span = {first: node.row, last: node.row};
      node._lineage = lineageOfGenomeName(tipName(node));
      continue;
    }
    let sum = 0;
    let first = Infinity;
    let last = -Infinity;
    const lineage = node.children[0]._lineage;
    let same = true;
    for (const child of node.children) {
      sum += child.row;
      if (child._span.first < first) first = child._span.first;
      if (child._span.last > last) last = child._span.last;
      if (child._lineage !== lineage) same = false;
    }
    node.row = sum / node.children.length;
    node._span = {first, last};
    node._lineage = same ? lineage : "mixed";
  }
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
      ? `0 / ${formatInt(visibleClusters.length)}`
      : `${formatInt(at + 1)} / ${formatInt(visibleClusters.length)}`;
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
  render(true);
}
/* The zoom sliders hold the logarithm of the zoom, so every step of the slider
   or notch of the wheel is the same ratio: 0.25 to 0.31 feels like 40 to 50. */
function readZoom(control, min, max) {
  const value = Math.exp(Number(control.value));
  return Number.isFinite(value) ? Math.max(min, Math.min(max, value)) : 8;
}
function writeZoom(control, value, min, max) {
  control.value = String(Math.log(Math.max(min, Math.min(max, value))));
}
function currentCellSize() {
  return readZoom(sizeControl, MIN_COL_ZOOM, MAX_COL_ZOOM);
}
function currentRowHeight() {
  return readZoom(rowHeightControl, MIN_ROW_ZOOM, MAX_ROW_ZOOM);
}
function clusterLabel(cluster) {
  return cluster.annotation ? `${cluster.name} - ${cluster.annotation}` : cluster.name;
}
/** Horizontal room the branches actually get, once labels are accounted for. */
function treeUsableWidth() {
  // Reserve room exactly when names are drawn (see renderTree), whatever the
  // tip count, or a zoomed-in large tree truncates every name to "L4…".
  const labelPad = labelToggle.checked && cellH >= TIP_LABEL_MIN_ROW ? 150 : 24;
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
function formatBranchLength(value) {
  const n = Number(value || 0);
  if (!Number.isFinite(n)) return "0";
  if (n === 0) return "0";
  if (Math.abs(n) < 0.001) return n.toExponential(2);
  return n.toPrecision(4);
}
/* One number format for the whole page, matching the Rust side: integers
   grouped with commas whatever the reader's locale, coordinates never grouped.
   A Spanish locale printed 12.715 for a coordinate next to 0.005 subs/site. */
function formatInt(value) {
  return Number(value).toLocaleString("en-US");
}
/** A share that never rounds a real count to 0% or 100%. */
function formatPct(value, total) {
  if (!total) return "0%";
  const pct = value / total * 100;
  if (value > 0 && pct < 0.1) return "<0.1%";
  if (value < total && pct > 99.9) return ">99.9%";
  return `${pct.toFixed(1)}%`;
}
function plural(count, noun) {
  return `${formatInt(count)} ${noun}${count === 1 ? "" : "s"}`;
}
function escapeHtml(text) {
  return String(text).replace(/[&<>"']/g, ch => ({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;","'":"&#39;"}[ch]));
}
/** Cut from the middle, keeping both ends: contig names differ at the end. */
function clipMiddle(text, maxLen) {
  const value = String(text || "");
  if (value.length <= maxLen) return value;
  const keep = Math.max(2, maxLen - 1);
  const head = Math.ceil(keep / 2);
  return `${value.slice(0, head)}…${value.slice(value.length - (keep - head))}`;
}
/** A gene in no cluster: Panaroo writes NA for it. */
function isUnclustered(gene) {
  return !gene.gap && (!gene.cluster || gene.cluster === "NA");
}
function referenceContextName() {
  return data.referenceGenome && data.referenceGenome.name
    ? `${data.referenceGenome.name} context`
    : "Reference context";
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
  return [{name:"__GAP__", cluster:"__GAP__", geneId:"reference context not resolved", strand:".", contig:"", start:"", end:"", length:"", target:false, reference:true, gap:true}];
}
function syntenyDisplayLabel(gene) {
  if (gene.gap) return gene.geneId || "context gap";
  if (isUnclustered(gene)) return "unclustered";
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
  if (bp >= 1000) return `${(bp / 1000).toFixed(1).replace(/\.0$/, "")} kb`;
  return `${Math.round(bp)} bp`;
}
function renderSyntenyPanel() {
  if (syntenyHidden) return;
  if (!syntenyPanel || !syntenySvg || !syntenySummary) return;
  const cluster = data.clusters[selected];
  const blocks = syntenyBlocksFor(cluster);
  updateSelectedSyntenyCarriers(blocks);
  const scrollContainer = syntenyPanel.querySelector(".synteny-scroll") || syntenyPanel;
  const paneWidth = Math.max(360, scrollContainer.clientWidth || syntenyPanel.clientWidth || 600);
  const nextKey = cluster
    ? `${cluster.name}|${blocks.length}|${selectedSyntenyBlock}|${selectedSyntenyReference}|${syntenyZoom().toFixed(2)}|${paneWidth}`
    : "none";
  if (nextKey === syntenyKey) return;
  syntenyKey = nextKey;

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
  // Blocks come from annotated loci, so a corrected cell (a gene the input's
  // annotation missed) has no block of its own: say so where the count is.
  const corrected = countsForCluster(selected).R;
  const correctedNote = corrected ? ` ${plural(corrected, "corrected genome")} ${corrected === 1 ? "has" : "have"} no annotated locus and no block; click their cells for the locus.` : "";
  syntenySummary.textContent = `${cluster.name}: ${plural(blocks.length, "unique 5-neighbour block")}; selected ${selectedLabel} marks ${plural(selectedCount, "tip")}.${referenceSummary}${correctedNote}`;

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
      parts.push(`<text x="${PADDING_X + 12}" y="${rowY + 36}" class="synteny-meta-title">${escapeHtml(clipText(referenceContextName(), 26))}</text>`);
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
      // The share is of genomes with an annotated locus, not of every carrier:
      // say so, or 7 of 9 printed as 78% next to a card that says 100%.
      const frac = Number(block.fraction);
      const annotated = frac > 0 && Number(block.nGenomes) > 0 ? Math.round(Number(block.nGenomes) / frac) : null;
      parts.push(
        `<text x="${PADDING_X + 12 + 26}" y="${rowY + 18}" class="synteny-meta-rank-sub">n=${block.nGenomes ?? "?"}${annotated ? ` of ${annotated} annotated` : fraction}</text>`
      );
      const titleText = block.representativeGenome
        ? `${block.representativeGenome}`
        : (firstExampleLocus(block) || cluster.name);
      parts.push(`<text x="${PADDING_X + 12}" y="${rowY + 36}" class="synteny-meta-title">${escapeHtml(clipText(titleText, 22))}</text>`);
      // Coordinates are never cut; a long contig name gives way in the middle.
      let coordText = firstExampleLocus(block) || "";
      if (block.representativeContig && block.representativeStart && block.representativeEnd) {
        const span = `:${formatCoord(block.representativeStart)}-${formatCoord(block.representativeEnd)}`;
        coordText = `${clipMiddle(block.representativeContig, Math.max(6, 30 - span.length))}${span}`;
      }
      parts.push(
        `<text x="${PADDING_X + 12}" y="${rowY + 52}" class="synteny-meta-sub"><title>${escapeHtml(`${block.representativeContig || ""}:${block.representativeStart || ""}-${block.representativeEnd || ""}${block.representativeLength ? ` (${block.representativeLength} bp)` : ""}`)}</title>${escapeHtml(coordText)}</text>`
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
      const title = isUnclustered(gene)
        ? ["Unclustered: this gene is in no cluster of the input matrix", gene.geneId, coord].filter(Boolean).join(" | ")
        : [syntenyGeneLabel(gene.cluster || gene.name), gene.geneId, coord].filter(Boolean).join(" | ");
      const geneAttrs = gene.gap || !gene.cluster || gene.cluster === "NA"
        ? ""
        : ` data-gene-cluster="${escapeHtml(gene.cluster)}"`;
      parts.push(
        `<path d="${arrowPath(x, baseY, geneW, GENE_H, gene.strand)}" class="${cls.join(" ")}"${geneAttrs}><title>${escapeHtml(title)}</title></path>`
      );
      const labelClass = gene.target
        ? "synteny-target-label"
        : isUnclustered(gene) ? "synteny-arrow-label synteny-unclustered" : "synteny-arrow-label";
      parts.push(
        `<text x="${x + geneW / 2}" y="${baseY + 14}" text-anchor="middle" class="${labelClass}">${escapeHtml(clipText(label, 11))}</text>`
      );
      // Narrow arrows keep only the name; the annotation is in the title.
      const subLabel = geneW >= 90 && !isUnclustered(gene) ? clipText(syntenySecondaryLabel(gene), 12) : "";
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
  return String(Math.round(n));
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
  // A bottom sheet (narrow windows) floats over the workspace and takes no width.
  const inspectorStyle = inspectorPane ? getComputedStyle(inspectorPane) : null;
  const inspectorShown = !!inspectorStyle && inspectorStyle.display !== "none" && inspectorStyle.position !== "fixed";
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
    resetMatrixCache();
    syntenyKey = "";
    render(false);
    event.preventDefault();
  });
}
function treeCaption() {
  // Read the mode, not the tips: an L8 tip in the tree said "L8-rooted" under
  // every root.
  const mode = rootMode ? rootMode.value : defaultRootMode;
  const rootLabel = mode === "l8"
    ? "L8-rooted"
    : mode === "midpoint"
      ? "midpoint-rooted"
      : mode.startsWith("genome:")
        ? `rooted on ${data.genomes[Number(mode.slice(7))] || "a genome"}`
        : "rooted as supplied";
  return treeMode === "phylogram" && maxDistance > 0
    ? `phylogram, ${rootLabel}, longest root-to-tip path ${formatBranchLength(maxDistance)}`
    : `cladogram, ${rootLabel}, branch lengths ignored`;
}
function updateTopBarMeta() {
  if (topTreeMeta) {
    // Short in the header, full text on hover: the caption is the first thing
    // that should give way when the pane narrows.
    // Matrix genomes the tree does not have are not drawn; say so instead of
    // letting "10 genomes" in the top bar and 4 rows on screen disagree.
    const offTree = offTreeGenomes;
    topTreeMeta.textContent = offTree.length
      ? `${plural(leafOrder.length, "tip")} · ${formatInt(offTree.length)} not on the tree`
      : plural(leafOrder.length, "tip");
    topTreeMeta.classList.toggle("warn", offTree.length > 0);
    const missing = offTree.length
      ? ` · not on the tree, so not drawn: ${offTree.slice(0, 20).join(", ")}${offTree.length > 20 ? `, and ${formatInt(offTree.length - 20)} more` : ""}`
      : "";
    topTreeMeta.setAttribute("data-tip", `${plural(leafOrder.length, "matched tip")} · ${treeCaption()}${missing}`);
  }
  if (topMatrixMeta) topMatrixMeta.textContent = `${formatInt(visibleClusters.length)} / ${plural(data.clusters.length, "cluster")}`;
  if (topSyntenyMeta) {
    const cluster = data.clusters[selected];
    const count = syntenyBlockCount(selected);
    topSyntenyMeta.textContent = cluster
      ? (count ? `${cluster.name} – ${plural(count, "unique block")}` : `${cluster.name} – no synteny data`)
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
  if (!bits.length && !data.referenceGenome) {
    // Without --reference-genome nothing was looked up, so "no equivalent
    // detected" claimed a search that never ran.
    return `<dt>Reference (audit)</dt><dd class="muted">No reference genome supplied (--reference-genome).</dd>`;
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
  const w = Math.max(1, Math.round(cssW * scale));
  const h = Math.max(1, Math.round(cssH * scale));
  // Reallocating a canvas on every frame cost a layout and a fresh bitmap;
  // only do it when the size actually changed.
  if (target.width !== w || target.height !== h) {
    target.width = w;
    target.height = h;
    target.style.width = `${cssW}px`;
    target.style.height = `${cssH}px`;
  }
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
  // The followers get some slack: their viewports differ from the matrix's by
  // a scrollbar or a strip, and without it the tree stopped a few pixels short
  // of the matrix's last row. They draw from the matrix's scroll position, so
  // the slack is never seen.
  if (headerSpacer) {
    headerSpacer.style.width = `${virtualW + 64}px`;
    headerSpacer.style.height = "1px";
  }
  if (treeSpacer) {
    treeSpacer.style.width = "1px";
    treeSpacer.style.height = `${virtualH + 64}px`;
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

function inDensityMode(cw, ch) {
  // A hair of slack: the opening fit lands exactly on the threshold, and a
  // rounding error must not tip it into density.
  return cw < DENSITY_BELOW_COL - 1e-6 || ch < DENSITY_BELOW_ROW - 1e-6;
}
function hexRgb(hex) {
  const h = String(hex).trim();
  if (h[0] !== "#" || h.length < 7) return [128, 128, 128];
  return [1, 3, 5].map(i => parseInt(h.slice(i, i + 2), 16));
}
/* The matrix as density, one device pixel at a time, written as ImageData at
   (dx, dy) device pixels. Columns whose set covers every row (or every row but
   a few) are added per pixel row in bulk, so the cost follows pixels plus set
   entries, not cells: the whole 10,000 x 30,000 matrix in one pass. */
function paintDensity(context, opts) {
  const {firstCol, lastCol, firstRow, lastRow, left, top, cw, ch, dataW, dataH, scale, dx = 0, dy = 0} = opts;
  const W = Math.max(1, Math.round(dataW * scale));
  const H = Math.max(1, Math.round(dataH * scale));
  const nPix = W * H;
  const P = new Uint32Array(nPix);
  const R = new Uint32Array(nPix);
  const K = new Uint32Array(nPix);
  const N = new Uint32Array(nPix);
  const colsAt = new Uint32Array(W);
  const bulk = new Int32Array(W);
  // Each visible row covers [y0, y1) device rows.
  const nRows = Math.max(0, lastRow - firstRow);
  const y0s = new Int32Array(nRows);
  const y1s = new Int32Array(nRows);
  const rowsAt = new Uint32Array(H);
  for (let i = 0; i < nRows; i += 1) {
    const r = firstRow + i;
    const a = Math.max(0, Math.min(H - 1, Math.floor((r * ch - top) * scale)));
    const b = Math.max(a + 1, Math.min(H, Math.floor(((r + 1) * ch - top) * scale)));
    y0s[i] = a;
    y1s[i] = b;
    for (let y = a; y < b; y += 1) rowsAt[y] += 1;
  }
  const addRow = (arr, x0, x1, r, sign) => {
    const i = r - firstRow;
    if (i < 0 || i >= nRows) return;
    for (let y = y0s[i]; y < y1s[i]; y += 1) {
      const base = y * W;
      for (let x = x0; x < x1; x += 1) arr[base + x] += sign;
    }
  };
  for (let c = firstCol; c < lastCol; c += 1) {
    const clusterIdx = originalClusterAtVisibleColumn(c);
    if (clusterIdx < 0) continue;
    const x0 = Math.max(0, Math.min(W - 1, Math.floor((c * cw - left) * scale)));
    const x1 = Math.max(x0 + 1, Math.min(W, Math.floor(((c + 1) * cw - left) * scale)));
    for (let x = x0; x < x1; x += 1) colsAt[x] += 1;
    const sets = clusterSets(clusterIdx);
    const rowsets = clusterRowSets(clusterIdx);
    const targets = [P, R, N, K];
    for (let state = 0; state < 4; state += 1) {
      const set = sets[state];
      const rows = rowsets[state];
      const arr = targets[state];
      if (set.kind === "list") {
        if (!rows || !rows.length) continue;
        let i = lowerBound(rows, firstRow);
        for (; i < rows.length && rows[i] < lastRow; i += 1) addRow(arr, x0, x1, rows[i], 1);
        continue;
      }
      // "all" or "complement": count the column in bulk, then take out the
      // excluded rows. Only input presence is ever that dense.
      if (arr !== P) {
        for (let r = firstRow; r < lastRow; r += 1) addRow(arr, x0, x1, r, 1);
        if (set.kind === "complement" && rows) {
          let i = lowerBound(rows, firstRow);
          for (; i < rows.length && rows[i] < lastRow; i += 1) addRow(arr, x0, x1, rows[i], -1);
        }
        continue;
      }
      for (let x = x0; x < x1; x += 1) bulk[x] += 1;
      if (set.kind === "complement" && rows) {
        let i = lowerBound(rows, firstRow);
        for (; i < rows.length && rows[i] < lastRow; i += 1) addRow(arr, x0, x1, rows[i], -1);
      }
    }
  }
  const absent = hexRgb(colors.A);
  const present = hexRgb(colors.P);
  const nsc = hexRgb(colors.N);
  const rescue = hexRgb(colors.R);
  const kept = hexRgb(colors.K);
  const panel = hexRgb(panelBg);
  const image = context.createImageData(W, H);
  const out = image.data;
  for (let y = 0; y < H; y += 1) {
    const rowsHere = rowsAt[y];
    for (let x = 0; x < W; x += 1) {
      const idx = y * W + x;
      const cells = colsAt[x] * rowsHere;
      let rgb = panel;
      if (cells > 0) {
        const p = Math.max(0, P[idx] + bulk[x] * rowsHere) / cells;
        const n = N[idx] / cells;
        const a = Math.max(0, 1 - p - n);
        rgb = [0, 1, 2].map(k => absent[k] * a + present[k] * p + nsc[k] * n);
        // Any finding in the pixel tints it, and the tint deepens with the
        // finding's share: one corrected cell among 500 is a hint, a pixel
        // where a twentieth of the cells were corrected is full colour.
        const finding = R[idx] ? rescue : K[idx] ? kept : null;
        if (finding) {
          const share = (R[idx] || K[idx]) / cells;
          const w = 0.3 + 0.7 * Math.min(1, share / 0.05);
          rgb = [0, 1, 2].map(k => rgb[k] * (1 - w) + finding[k] * w);
        }
      }
      const o = idx * 4;
      out[o] = rgb[0];
      out[o + 1] = rgb[1];
      out[o + 2] = rgb[2];
      out[o + 3] = 255;
    }
  }
  context.putImageData(image, dx, dy);
}
/* Paints the cells of columns [firstCol, lastCol) and rows [firstRow, lastRow),
   with the matrix origin at (left, top) and cw x ch cells. The screen and the
   exported figure both draw through here, so they cannot disagree. */
function paintCells(context, opts) {
  const {firstCol, lastCol, firstRow, lastRow, left, top, cw, ch, dataW, dataH} = opts;
  // Absent cells have their own fill, so the matrix has an edge and an empty
  // run does not look like an empty pane.
  context.fillStyle = colors.A;
  context.fillRect(0, 0, dataW, dataH);
  // Cells never shrink below a hairline, so a dense matrix stays legible instead
  // of fading out when several clusters share a pixel column.
  const drawW = Math.max(1, cw > 2 ? cw - 1 : cw);
  const drawH = Math.max(1, ch > 2 ? ch - 1 : ch);
  // A gutter exists only when drawH < ch. At exactly 2px rows there is none,
  // and drawing row by row there was pure cost.
  const gutterRows = drawH < ch;
  const byProcess = matrixColorMode.value === "process";

  // Runs of adjacent rows are one rect each; row gutters, when rows are tall
  // enough to have them, are cut afterwards as one background line per visible
  // row. Drawing every cell as its own rect cost 340,000 fillRect calls for a
  // dense window at 3px rows; this is a few thousand.
  const fillRows = (x, a, b) => {
    context.fillRect(x, a * ch - top, drawW, (b - a) * ch);
  };
  const paintSet = (set, rows, x) => {
    if (set.kind === "list") {
      let i = lowerBound(rows, firstRow);
      while (i < rows.length && rows[i] < lastRow) {
        const from = rows[i];
        let to = from + 1;
        i += 1;
        while (i < rows.length && rows[i] === to && to < lastRow) { to += 1; i += 1; }
        fillRows(x, from, to);
      }
      return;
    }
    // "all" or "complement": paint the visible runs between excluded rows. No
    // background punch-out, so a later state never erases an earlier one.
    const excluded = set.kind === "all" ? null : rows;
    const n = excluded ? excluded.length : 0;
    let i = n ? lowerBound(excluded, firstRow) : 0;
    let from = firstRow;
    while (from < lastRow) {
      const stop = i < n ? Math.min(excluded[i], lastRow) : lastRow;
      if (stop > from) fillRows(x, from, stop);
      from = stop + 1;
      i += 1;
    }
  };
  // Input and NSC cells first for the whole window, then evidence kept out of
  // the matrix, then rescued cells last. Below one pixel per column neighbours
  // share pixels, and painting column by column let the next cluster's cells
  // cover the corrections, which are the cells this report exists to show.
  for (const states of [[0, 2], [3], [1]]) {
    for (let c = firstCol; c < lastCol; c += 1) {
      const clusterIdx = originalClusterAtVisibleColumn(c);
      if (clusterIdx < 0) continue;
      const sets = clusterSets(clusterIdx);
      const rowsets = clusterRowSets(clusterIdx);
      // Below a pixel per column, snap to whole pixels: a fractional x is
      // antialiased into a blend of neighbours, and a green cell blended into
      // blue reads as neither.
      const x = cw < 1 ? Math.floor(c * cw - left) : c * cw - left;
      const processFill = byProcess
        ? processColor(processKind(data.clusters[clusterIdx]))
        : null;
      for (const state of states) {
        const set = sets[state];
        if (set.kind === "list" && !set.idx.length) continue;
        // In mechanism mode the input cells take the mechanism's colour, and
        // the cells Recore changed keep theirs, so they still read as findings.
        const own = colors[STATE_CODES[state]] || colors.A;
        context.fillStyle = processFill && (state === 0 || state === 2) ? processFill : own;
        paintSet(set, rowsets[state], x);
      }
    }
  }

  if (gutterRows) {
    context.fillStyle = panelBg;
    for (let r = firstRow; r < lastRow; r += 1) {
      context.fillRect(0, r * ch - top + drawH, dataW, ch - drawH);
    }
  }
  if (drawW < cw) {
    // Column gutters too, or absent cells would run together into a band.
    context.fillStyle = panelBg;
    for (let c = firstCol; c < lastCol; c += 1) {
      context.fillRect(c * cw - left + drawW, 0, cw - drawW, dataH);
    }
  }
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
  const dataW = Math.min(view.w, virtualW - left);
  const dataH = Math.min(view.h, virtualH - top);
  const [firstCol, lastCol] = visibleColumnRange(view.w);
  const [firstRow, lastRow] = visibleRowRange(view.h);
  if (rows && cols && dataW > 0 && dataH > 0) {
    const opts = {firstCol, lastCol, firstRow, lastRow, left, top, cw: cellW, ch: cellH, dataW, dataH};
    if (inDensityMode(cellW, cellH)) paintDensity(ctx, {...opts, scale});
    else paintCells(ctx, opts);
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
  // Gate on the cell, not the run: only the columns in view are labelled, so a
  // 30,000-cluster run zoomed to 10px columns labels about 80 of them. The old
  // count cap meant a large run never showed a single cluster name.
  if (cellW < LABEL_MIN_CELL_W) return 0;
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

  // Ticks that touch and share a colour are one bar: a rect per column cost a
  // fillRect per cluster in view, thousands of them, on every repaint.
  if (!byProcess && cellW <= 2) {
    headerCtx.fillStyle = gridStroke;
    const x0 = firstCol * cellW - left;
    headerCtx.fillRect(x0, cssH - 9, Math.max(1, (lastCol - firstCol) * cellW), 6);
  }
  for (let c = firstCol; c < lastCol && (byProcess || cellW > 2); c += 1) {
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
/* Longest prefix that fits, found by halving rather than one character per
   measureText call: 120-character names cost thousands of calls a header.
   Results are cached per name and width; the font never changes. */
const fitCache = new Map();
function fitText(context, text, maxWidth) {
  const value = String(text || "");
  const key = `${Math.round(maxWidth)}|${value}`;
  const cached = fitCache.get(key);
  if (cached !== undefined) return cached;
  let out = value;
  if (context.measureText(value).width > maxWidth) {
    let lo = 1;
    let hi = value.length - 1;
    out = "…";
    while (lo <= hi) {
      const mid = (lo + hi) >> 1;
      const candidate = `${value.slice(0, mid)}…`;
      if (context.measureText(candidate).width <= maxWidth) {
        out = candidate;
        lo = mid + 1;
      } else {
        hi = mid - 1;
      }
    }
  }
  if (fitCache.size > 20000) fitCache.clear();
  fitCache.set(key, out);
  return out;
}
function fitLabel(text, maxWidth) {
  return fitText(headerCtx, text, maxWidth);
}

/* The tree is a canvas rather than SVG: 10,000 tips is roughly 40,000 SVG nodes,
   which the browser will not lay out at an interactive frame rate. Only the tips
   inside the visible band are drawn, and internal edges are clipped to it. */
/* ---- Rerooting ------------------------------------------------------------
   The displayed root drives every ancestral reading in the panel, so it has to
   be the reader's choice rather than whatever the Newick happened to carry.
   assignTreeIds already yields an adjacency list, so rooting anywhere is a
   matter of picking the node and rebuilding downwards from it. */
function graphOf(nodes) {
  const adj = nodes.map(() => []);
  for (const node of nodes) {
    for (const child of node.children || []) {
      const length = Number(child.length || 0);
      adj[node._id].push({id: child._id, length});
      adj[child._id].push({id: node._id, length});
    }
  }
  return {adj};
}
function treeGraph(node) {
  const root = cloneNode(node);
  const nodes = [];
  assignTreeIds(root, nodes, []);
  return {root, nodes, adj: graphOf(nodes).adj};
}
function buildFrom(nodes, adj, rootId, parentId = null, length = 0) {
  const root = copyNode(nodes[rootId]);
  root.length = length;
  const stack = [[rootId, parentId, root]];
  while (stack.length) {
    const [id, from, target] = stack.pop();
    for (const edge of adj[id]) {
      if (edge.id === from) continue;
      const child = copyNode(nodes[edge.id]);
      child.length = edge.length;
      target.children.push(child);
      stack.push([edge.id, id, child]);
    }
  }
  return root;
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
   the root goes exactly halfway along the path between them, splitting the
   edge that straddles the midpoint. Rooting on the nearer end of that edge put
   every root distance off by up to half an edge. */
function midpointTree(nodes, adj) {
  if (nodes.length <= 1) return buildFrom(nodes, adj, 0);
  const first = farthestFrom(adj, 0).best;
  const {best: second, dist, prev} = farthestFrom(adj, first);
  const half = dist[second] / 2;
  if (!(half > 0)) {
    // No lengths to halve: root beside the tip, never on it, so it stays a tip.
    const next = adj[first][0];
    return buildFrom(nodes, adj, next ? next.id : first);
  }
  let cursor = second;
  while (prev[cursor] >= 0 && dist[prev[cursor]] > half) cursor = prev[cursor];
  const parent = prev[cursor];
  const eps = 1e-12 * Math.max(1, half);
  if (parent < 0 || Math.abs(dist[cursor] - half) <= eps) return buildFrom(nodes, adj, cursor);
  if (Math.abs(dist[parent] - half) <= eps) return buildFrom(nodes, adj, parent);
  return {
    label: "",
    genome: null,
    length: 0,
    children: [
      buildFrom(nodes, adj, cursor, parent, dist[cursor] - half),
      buildFrom(nodes, adj, parent, cursor, half - dist[parent]),
    ],
  };
}
function rerootTree(mode) {
  if (mode === "input") return tidyTree(cloneNode(data.tree), "none");
  if (mode === "l8") return tidyTree(rootByL8(data.tree), "l8");
  const {nodes, adj} = treeGraph(data.tree);
  if (mode === "midpoint") return tidyTree(midpointTree(nodes, adj), "tips");
  if (mode.startsWith("genome:")) {
    const genomeIdx = Number(mode.slice(7));
    const id = nodes.findIndex(n => n.genome === genomeIdx);
    if (id >= 0) {
      // Root on the tip's parent, so the chosen genome hangs off the root as an
      // outgroup rather than becoming a degenerate root itself.
      const parent = adj[id][0];
      return tidyTree(buildFrom(nodes, adj, parent ? parent.id : id, null, 0), "tips");
    }
  }
  return tidyTree(hasL8 ? rootByL8(data.tree) : cloneNode(data.tree), hasL8 ? "l8" : "none");
}
/** Rebuild layout state after the displayed root changes. */
function applyRoot(mode) {
  treeRoot = rerootTree(mode);
  layoutTree(treeRoot);
  offTreeGenomes = data.genomes.filter((_, gi) => !leafByGenome.has(gi));
  buildRowOfGenome();
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
  // Only tips in the visible band are drawn, so the tip count is no reason to
  // hide names once a row is tall enough to carry one.
  const showLabels = labelToggle.checked && cellH >= TIP_LABEL_MIN_ROW;

  treeCtx.lineWidth = 1;
  treeCtx.strokeStyle = gridStroke;
  drawTreeNode(treeRoot, top, firstRow, lastRow, view.h, showLabels);
  renderLineageStrip(view.h, top, firstRow, lastRow);
}
function drawTreeNode(root, top, firstRow, lastRow, viewH, showLabels) {
  const stack = [root];
  while (stack.length) {
    const node = stack.pop();
    // Whole subtrees outside the band cost nothing.
    if (node._span.last < firstRow || node._span.first > lastRow) continue;
    const x = xForNode(node);
    const y = node.row * cellH + cellH / 2 - top;
    if (!isLeaf(node)) {
      for (const child of node.children) {
        if (child._span.last < firstRow || child._span.first > lastRow) continue;
        const cx = xForNode(child);
        const cy = child.row * cellH + cellH / 2 - top;
        treeCtx.strokeStyle = treeEdgeColor(child);
        treeCtx.beginPath();
        treeCtx.moveTo(x, cy);
        treeCtx.lineTo(cx, cy);
        treeCtx.moveTo(x, y);
        treeCtx.lineTo(x, cy);
        treeCtx.stroke();
        stack.push(child);
      }
      continue;
    }
    if (node.row < firstRow || node.row >= lastRow) continue;
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
}
function fitTipLabel(text) {
  const budget = Math.max(24, treeWidth - xForNode({depth: maxDepth, distance: maxDistance, children: []}) - 10);
  return fitText(treeCtx, text, budget);
}
/** Row span covered by a subtree, memoised: used to skip offscreen branches. */
function nodeRowSpan(node) {
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
  const share = formatPct(value, total);
  return `<div class="count-row"><span class="sw ${key}"></span>` +
    `<span class="count-name">${escapeHtml(label)}</span>` +
    `<b class="num">${formatInt(value)}</b><small class="num">${share}</small></div>`;
}

/** The states this run can have, in drawing order. NSC only when it ran, and
    kept-out evidence only when there is any: a row of zeros for a step that
    never ran reads as a result. */
function shownStates() {
  const states = ["P", "R"];
  if (data.nscEnabled) states.push("N");
  if (anyKeptOut) states.push("K");
  states.push("A");
  return states;
}
function stateCountRows(counts, total) {
  return shownStates()
    .map(state => countRow(labels[state], counts[state] || 0, total, state.toLowerCase()))
    .join("");
}
/** Prevalence as a bar, because "16 / 16" and "3 / 16" should not look alike. */
function prevalenceBar(present, total, pct) {
  const width = total ? Math.max(1, Math.round(present / total * 100)) : 0;
  return `<div class="prevalence-bar"><i style="width:${width}%"></i></div>` +
    `<p class="prevalence num">${formatInt(present)} / ${plural(total, "genome")} &middot; ${pct}</p>`;
}

/* The inspector reflects the selected cluster; keyed so it is not rebuilt on
   every scroll frame. */
function renderClusterDetails() {
  const cluster = data.clusters[selected];
  const nextDetailsKey = cluster
    ? `${selected}|${selectedGenome}|${referenceOnly}|${mechanismFilterKey()}|${visibleClusters.length}`
    : "none";
  if (nextDetailsKey !== detailsKey && cluster) {
    const counts = countsForCluster(selected);
    const mechanism = processKind(cluster);
    const mechanismLabel = processLabel(mechanism);
    const mechanismDescription = processDescription(mechanism);
    const mechanismColor = processCss(mechanism);
    const prevalencePct = formatPct(cluster.present, data.genomes.length);
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
    details.innerHTML = `${cellCardHtml()}
      <div class="cluster-head">
        <h3>${escapeHtml(cluster.name)}</h3>
        <p class="annotation">${escapeHtml(cluster.annotation || "No annotation provided")}</p>
      </div>
      ${prevalenceBar(cluster.present, total, prevalencePct)}
      <div class="card-title">Evidence in ${plural(total, "genome")}</div>
      <div class="count-rows">${stateCountRows(counts, total)}</div>
      ${lineageBreakdownHtml(selected)}
      ${changedCellsHtml(selected)}
      <div class="card-title">Mechanism</div>
      <span class="process-chip" style="--process-color:${mechanismColor}"><i></i>${escapeHtml(mechanismLabel)}</span>
      <p class="process-wording">${escapeHtml(mechanismDescription)}</p>
      <dl class="cluster-dl">
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
  renderLineageLegend();
  renderCellLegend();
  renderTreeScaleBar();
  renderOverview();
  updateTopBarMeta();
  renderFilterChips();
  renderOverlays();
  if (afterRender) afterRender();
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
/* Marks a column or row with lines that sit in the gutters when cells have
   gutters, and just outside the cells when they do not, so no cell is tinted:
   a selected absent cell still reads as absent. */
function markColumn(context, x, w, h, width) {
  if (w >= 3) {
    context.fillRect(x - width, 0, width, h);
    context.fillRect(x + w - 1, 0, width, h);
    return;
  }
  const bx = Math.floor(x);
  context.fillRect(bx - width, 0, width, h);
  context.fillRect(bx + Math.max(1, Math.ceil(w)), 0, width, h);
}
function markRow(context, y, h, w, width) {
  if (h >= 3) {
    context.fillRect(0, y - width, w, width);
    context.fillRect(0, y + h - 1, w, width);
    return;
  }
  const by = Math.floor(y);
  context.fillRect(0, by - width, w, width);
  context.fillRect(0, by + Math.max(1, Math.ceil(h)), w, width);
}
function clearOverlay(target, context, base) {
  if (!target || !context || !base) return false;
  const w = base.clientWidth;
  const h = base.clientHeight;
  if (target.width !== base.width || target.height !== base.height) {
    sizeCanvas(target, context, w, h, base.width / Math.max(1, w));
  } else {
    context.setTransform(1, 0, 0, 1, 0, 0);
    context.clearRect(0, 0, target.width, target.height);
    context.setTransform(base.width / Math.max(1, w), 0, 0, base.width / Math.max(1, w), 0, 0);
  }
  return true;
}
function renderOverlays() {
  const left = matrixScrollLeft();
  const top = matrixScrollTop();
  const selectedCol = selectedVisibleColumn();
  if (clearOverlay(matrixOverlay, matrixOverlayCtx, canvas)) {
    const o = matrixOverlayCtx;
    const w = Math.min(canvas.clientWidth, virtualW - left);
    const h = Math.min(canvas.clientHeight, virtualH - top);
    if (hoverRow >= 0 && h > 0) {
      o.fillStyle = hoverStroke;
      markRow(o, hoverRow * cellH - top, cellH, w, 1);
    }
    if (hoverCol >= 0 && h > 0) {
      o.fillStyle = hoverStroke;
      markColumn(o, hoverCol * cellW - left, cellW, h, 1);
    }
    if (selectedCol >= 0 && h > 0) {
      o.fillStyle = markInk;
      markColumn(o, selectedCol * cellW - left, cellW, h, cellW >= 3 ? 2 : 1);
    }
  }
  if (clearOverlay(headerOverlay, headerOverlayCtx, headerCanvas) && hoverCol >= 0) {
    const o = headerOverlayCtx;
    o.fillStyle = hoverStroke;
    markColumn(o, hoverCol * cellW - left, cellW, headerCanvas.clientHeight, 1);
  }
  if (clearOverlay(treeOverlay, treeOverlayCtx, treeCanvas) && hoverRow >= 0) {
    const o = treeOverlayCtx;
    o.fillStyle = hoverStroke;
    markRow(o, hoverRow * cellH - top, cellH, treeCanvas.clientWidth, 1);
  }
}
let overlayQueued = false;
function scheduleOverlays() {
  if (overlayQueued) return;
  overlayQueued = true;
  requestAnimationFrame(() => {
    overlayQueued = false;
    renderOverlays();
  });
}
function setHover(row, col, text) {
  const changed = hoverRow !== row || hoverCol !== col;
  hoverRow = row;
  hoverCol = col;
  if (text) status.textContent = text;
  if (changed) scheduleOverlays();
}
function setTreeMode(mode) {
  treeMode = mode;
  phylogramButton.classList.toggle("active", mode === "phylogram");
  cladogramButton.classList.toggle("active", mode === "cladogram");
  render(false);
}
function resetMatrixCache() {
  hoverCol = -1;
}
refreshThemeColors();
const hasL8 = treeHasTip(data.tree, isL8Tip);
const defaultRootMode = hasL8 ? "l8" : "input";
if (rootMode) {
  // No L8 tip, no L8 outgroup to offer.
  const l8Option = [...rootMode.options].find(option => option.value === "l8");
  if (l8Option && !hasL8) l8Option.remove();
}
// Lineage colouring only means something when some tip names a lineage.
if (!treeHasTip(data.tree, node => lineageOfGenomeName(tipName(node)) !== "unknown")) {
  const lineageOption = [...colorMode.options].find(option => option.value === "lineage");
  if (lineageOption) lineageOption.remove();
}
applyRoot(defaultRootMode);
if (rootMode) rootMode.value = defaultRootMode;
computeGenomeStats();
refreshVisibleClusters();
// Open at a zoom that suits the run. A 16-genome, 12-cluster matrix at the old
// fixed sliders was a thumbnail in a half-empty pane; a 3000-tip tree needed the
// opposite. Column width is fitted to the space the matrix pane actually has.
/* `whole` fits every row and column, as density when cells fall below a pixel.
   Without it (the opening view) cells stay at least the density threshold, so
   the first paint draws cells, not a raster of the whole run. */
function fitZoomToData(whole = false) {
  const tips = Math.max(1, leafOrder.length);
  const cols = Math.max(1, visibleClusters.length);
  const bodyW = viewerBody ? viewerBody.clientWidth : 0;
  // The inspector takes a column only when it is not a bottom sheet.
  const inspectorTakesWidth = !(narrowQuery && narrowQuery.matches);
  const matrixW = bodyW > 0
    ? Math.max(220, bodyW - layoutState.treeW - (inspectorTakesWidth ? layoutState.inspectorW : 0) - 24)
    : FIT_TARGET_W;
  const viewH = matrixViewport ? matrixViewport.clientHeight : 0;
  const targetH = viewH > 120 ? viewH - 4 : FIT_TARGET_H;
  // Whole pixels when the rows fit at one or more; otherwise exactly the share
  // that fits, so Z always shows every row (as density below a pixel).
  const rowFit = targetH / tips >= 1
    ? Math.max(1, Math.min(MAX_ROW_ZOOM, Math.floor(targetH / tips)))
    : whole ? Math.max(MIN_ROW_ZOOM, targetH / tips) : DENSITY_BELOW_ROW;
  // cellW is round(slider * 0.82), so undo that factor when fitting to a width.
  // Round down to the slider step (0.25), not to whole units: 1,118 columns in
  // 788px want 0.86, and flooring that to 0 sent the fit to the 0.25 minimum,
  // leaving half the pane empty after Z.
  const colStep = 0.25;
  const exactCol = (matrixW - 2) / cols / 0.82;
  const colFit = exactCol >= colStep
    ? Math.max(colStep, Math.min(MAX_COL_ZOOM, Math.floor(exactCol / colStep) * colStep))
    : whole ? Math.max(MIN_COL_ZOOM, exactCol) : Math.max(exactCol, DENSITY_BELOW_COL / 0.82);
  writeZoom(rowHeightControl, rowFit, MIN_ROW_ZOOM, MAX_ROW_ZOOM);
  writeZoom(sizeControl, colFit, MIN_COL_ZOOM, MAX_COL_ZOOM);
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
if (cellLegend) {
  cellLegend.addEventListener("click", event => {
    if (!event.target.closest("[data-open-mechanisms]")) return;
    event.stopPropagation();
    setPaneMenu(matrixMenu, true);
    if (processLegend) processLegend.scrollIntoView({block: "nearest"});
  });
}
document.addEventListener("pointerdown", event => {
  if (!anyPaneMenuOpen()) return;
  if (event.target.closest(".pane-menu, #treeMenuBtn, #matrixMenuBtn, [data-open-mechanisms]")) return;
  setPaneMenu(null, false);
});

colorMode.addEventListener("change", () => render(false));
labelToggle.addEventListener("change", () => render(false));
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
    render(false);
    return;
  }
  const row = event.target.closest(".synteny-row[data-block-index]");
  if (!row) return;
  selectedSyntenyBlock = Number(row.getAttribute("data-block-index") || 0);
  selectedSyntenyReference = false;
  syntenyKey = "";
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
/** The clusters in view, with what the matrix shows for each, as TSV. */
function visibleClustersTsv() {
  const lines = [["cluster", "annotation", "mechanism", "present_in_input", "corrected", "nsc", "kept_out", "absent"].join("\t")];
  for (const idx of visibleClusters) {
    const cluster = data.clusters[idx];
    const c = countsForCluster(idx);
    lines.push([cluster.name, cluster.annotation || "", processKind(cluster), c.P, c.R, c.N, c.K, c.A]
      .map(v => String(v).replace(/[\t\n]/g, " ")).join("\t"));
  }
  return `${lines.join("\n")}\n`;
}
/* Active mechanism filters as removable chips beside the matrix count, so a
   filter left on is visible with the menu closed. */
function renderFilterChips() {
  let box = document.getElementById("filterChips");
  if (!box && topMatrixMeta) {
    box = document.createElement("span");
    box.id = "filterChips";
    box.className = "filter-chips";
    topMatrixMeta.after(box);
  }
  if (!box) return;
  const key = mechanismFilterKey();
  if (box.dataset.key === key) return;
  box.dataset.key = key;
  box.innerHTML = [...mechanismFilters].sort().map(mechanism =>
    `<button type="button" class="filter-chip" data-remove-mechanism="${escapeHtml(mechanism)}" data-tip="Remove this filter"><i style="background:${processCss(mechanism)}"></i>${escapeHtml(processLabel(mechanism))} &times;</button>`).join("");
}
function applyMechanismFilters() {
  refreshVisibleClusters();
  populateOptions();
  resetMatrixCache();
  syncMenuIndicators();
  // A filter changes how many columns there are; fit them to the pane again.
  fitZoomToData(true);
  if (matrixPane) matrixPane.scrollLeft = 0;
  render(true);
}
document.addEventListener("click", event => {
  const chip = event.target.closest("[data-remove-mechanism]");
  if (!chip) return;
  mechanismFilters.delete(chip.getAttribute("data-remove-mechanism"));
  applyMechanismFilters();
});
processLegend.addEventListener("click", event => {
  if (event.target.closest("[data-download-visible]")) {
    downloadBlob(new Blob([visibleClustersTsv()], {type: "text/tab-separated-values"}), "recore_visible_clusters.tsv");
    return;
  }
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
  applyMechanismFilters();
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
    // Sorting is asked for to see the top of the new order; following the old
    // selection used to land the view thousands of columns into it.
    if (matrixPane) matrixPane.scrollLeft = 0;
    if (visibleClusters.length) {
      selected = originalClusterAtVisibleColumn(0);
      selectedSyntenyReference = false;
      populateOptions();
    }
    render(false);
  });
}
/* A plain wheel scrolls, everywhere: it used to zoom rows over the tree, so the
   tree could not be scrolled and one tick switched on the slow path. Ctrl or
   Cmd and the wheel (a trackpad pinch arrives that way) zoom columns, Alt and
   the wheel zoom rows, by 1.25 a notch, keeping the cell under the pointer
   where it is. Motion is read from whichever axis carries it, since macOS sends
   Shift-wheel as deltaX, and accumulated, so a trackpad swipe of 30 events is a
   few steps and not thirty. */
const ZOOM_STEP = 1.25;
const WHEEL_NOTCH_PX = 40;
let wheelAccum = 0;
let wheelAccumAt = 0;
function zoomAround(axis, factor, clientX, clientY) {
  const rect = matrixPane.getBoundingClientRect();
  if (axis === "column") {
    const px = Math.max(0, Math.min(rect.width, clientX - rect.left));
    const at = (matrixPane.scrollLeft + px) / cellW;
    const before = currentCellSize();
    writeZoom(sizeControl, before * factor, MIN_COL_ZOOM, MAX_COL_ZOOM);
    if (currentCellSize() === before) return false;
    updateCellGeometry();
    matrixPane.scrollLeft = Math.max(0, at * cellW - px);
  } else {
    const py = Math.max(0, Math.min(rect.height, clientY - rect.top));
    const at = (matrixPane.scrollTop + py) / cellH;
    const before = currentRowHeight();
    writeZoom(rowHeightControl, before * factor, MIN_ROW_ZOOM, MAX_ROW_ZOOM);
    if (currentRowHeight() === before) return false;
    updateCellGeometry();
    matrixPane.scrollTop = Math.max(0, at * cellH - py);
  }
  return true;
}
function wheelZoom(event) {
  const columnZoom = event.ctrlKey || event.metaKey;
  const rowZoom = !columnZoom && event.altKey;
  if (!columnZoom && !rowZoom) {
    // Shift and a vertical wheel scroll sideways. macOS already turns it into
    // deltaX; other systems and synthetic events do not, so do it here.
    if (event.shiftKey && Math.abs(event.deltaY) > Math.abs(event.deltaX)) {
      event.preventDefault();
      const unit = event.deltaMode === 1 ? 16 : event.deltaMode === 2 ? 400 : 1;
      matrixPane.scrollLeft += event.deltaY * unit;
    }
    return;
  }
  event.preventDefault();
  event.stopPropagation();
  const unit = event.deltaMode === 1 ? 16 : event.deltaMode === 2 ? 400 : 1;
  const motion = (Math.abs(event.deltaY) >= Math.abs(event.deltaX) ? event.deltaY : event.deltaX) * unit;
  const now = performance.now();
  if (now - wheelAccumAt > 300 || Math.sign(motion) !== Math.sign(wheelAccum)) wheelAccum = 0;
  wheelAccumAt = now;
  wheelAccum += motion;
  if (Math.abs(wheelAccum) < WHEEL_NOTCH_PX) return;
  const factor = wheelAccum < 0 ? ZOOM_STEP : 1 / ZOOM_STEP;
  wheelAccum = 0;
  if (zoomAround(columnZoom ? "column" : "row", factor, event.clientX, event.clientY)) {
    resetMatrixCache();
    render(false);
  }
}
viewerBody.addEventListener("wheel", wheelZoom, {passive:false});

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
    // The matrix moved the header; this is that move coming back, not a new one.
    if (syncingScroll || Math.abs(headerScroller.scrollLeft - matrixPane.scrollLeft) < 1) return;
    syncingScroll = true;
    matrixPane.scrollLeft = headerScroller.scrollLeft;
    syncingScroll = false;
    scheduleRender(false);
  }, {passive: true});
}
if (treeScroll) {
  treeScroll.addEventListener("scroll", () => {
    if (syncingScroll || Math.abs(treeScroll.scrollTop - matrixPane.scrollTop) < 1) return;
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
  openInspectorSheet();
  const next = originalClusterAtVisibleColumn(c);
  // A new cluster starts on its first synteny block, as stepping does; the
  // previous cluster's block number used to carry over.
  if (next !== selected) selectedSyntenyBlock = 0;
  selected = next;
  selectedSyntenyReference = false;
  populateOptions();
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
  // The cell's own evidence first, the cluster's mechanism after and labelled as
  // the cluster's, so it does not read as the reason for this cell.
  const evidence = state === "R" || state === "K" ? cellEvidence(clusterIdx, genomeIdx)[0] : null;
  const stateLine = evidence ? `${labels[state]}: ${evidenceSummary(evidence)}` : labels[state];
  setHover(r, c, `${data.genomes[genomeIdx]} / ${cluster.name}: ${stateLine}`);
  const key = `${clusterIdx}|${genomeIdx}`;
  if (tooltip.dataset.key !== key) {
    tooltip.dataset.key = key;
    tooltip.innerHTML = `<strong>${escapeHtml(data.genomes[genomeIdx])}</strong> / ${escapeHtml(cluster.name)}` +
      `<span class="tt-state"><i class="sw ${state.toLowerCase()}"></i>${escapeHtml(stateLine)}</span>` +
      `<span class="tt-mech">Cluster mechanism: ${escapeHtml(processLabel(mechanism))}</span>`;
  }
  tooltip.style.display = "block";
  // Below and right of the pointer, flipped at the pane edges.
  const paneRect = matrixPane.getBoundingClientRect();
  const x = event.clientX - paneRect.left;
  const y = event.clientY - paneRect.top;
  const w = tooltip.offsetWidth;
  const h = tooltip.offsetHeight;
  const left = x + 14 + w > paneRect.width ? Math.max(4, x - 14 - w) : x + 14;
  const top = y + 14 + h > paneRect.height ? Math.max(4, y - 14 - h) : y + 14;
  tooltip.style.transform = `translate(${Math.round(left)}px, ${Math.round(top)}px)`;
});
/* ---- Changed cells, for the drawer ----------------------------------------
   Every evidence row as a list a reader can copy or download, built from the
   same per-cluster data the cell card reads. The preview is built when the
   drawer first opens; the full text only when it is asked for. */
const CELL_STATE_WORDS = {P: "present", R: "corrected", N: "NSC", K: "kept out", A: "absent"};
const CORRECTION_COLUMNS = ["genome", "cluster", "cell_state", "decision", "rescue_type", "identity", "coverage", "contig", "start", "end", "strand", "synteny", "synteny_score", "orf_integrity_signal", "matrix_action", "publication_status", "rejection_reason", "risk_flags"];
function correctionRecord(clusterIdx, row) {
  const value = column => {
    if (column === "genome") return data.genomes[row.genome] || "";
    if (column === "cluster") return data.clusters[clusterIdx].name;
    if (column === "cell_state") return CELL_STATE_WORDS[stateAt(clusterIdx, row.genome)] || "";
    if (column === "identity" || column === "coverage") return row[column] === null ? "" : row[column].toFixed(1);
    if (column === "start" || column === "end") return row[column] === null ? "" : String(row[column]);
    return row[column] || "";
  };
  return CORRECTION_COLUMNS.map(value);
}
function correctionsTsv() {
  const lines = [CORRECTION_COLUMNS.join("\t")];
  for (let ci = 0; ci < data.clusters.length; ci += 1) {
    if (!evidenceCount(ci)) continue;
    for (const row of decodeEvidence(ci)) {
      lines.push(correctionRecord(ci, row).map(v => String(v).replace(/[\t\n]/g, " ")).join("\t"));
    }
  }
  return `${lines.join("\n")}\n`;
}
const CORRECTIONS_PREVIEW = 50;
let correctionsBuilt = false;
function buildCorrectionsPreview() {
  const box = document.getElementById("correctionsTable");
  if (!box || correctionsBuilt || !evidenceData) return;
  correctionsBuilt = true;
  let total = 0;
  for (let ci = 0; ci < data.clusters.length; ci += 1) total += evidenceCount(ci);
  if (!total) {
    box.innerHTML = "<p class='muted'>No evidence rows: Recore changed no cell in this run.</p>";
    return;
  }
  // Corrections first, then kept out, then rejected, in matrix order.
  const picked = [[], [], []];
  for (let ci = 0; ci < data.clusters.length && picked[0].length < CORRECTIONS_PREVIEW; ci += 1) {
    if (!evidenceCount(ci)) continue;
    for (const row of decodeEvidence(ci)) {
      const rank = row.matrix_action === "correct_presence" ? 0 : row.decision === "accepted" ? 1 : 2;
      if (picked[rank].length < CORRECTIONS_PREVIEW) picked[rank].push([ci, row]);
    }
  }
  const shown = [...picked[0], ...picked[1], ...picked[2]].slice(0, CORRECTIONS_PREVIEW);
  const body = shown.map(([ci, row]) => {
    const state = stateAt(ci, row.genome);
    return `<tr data-cell-cluster="${ci}" data-cell-genome="${row.genome}" tabindex="0">` +
      `<td><span class="state-chip ${state.toLowerCase()}">${escapeHtml(CELL_STATE_WORDS[state])}</span></td>` +
      `<td class="mono">${escapeHtml(data.genomes[row.genome] || "")}</td>` +
      `<td class="mono">${escapeHtml(data.clusters[ci].name)}</td>` +
      `<td>${escapeHtml(row.rescue_type || row.evidence_state || "")}</td>` +
      `<td class="num">${pctText(row.identity)}</td><td class="num">${pctText(row.coverage)}</td></tr>`;
  }).join("");
  box.innerHTML = `<p class="muted">${plural(total, "evidence row")}${total > shown.length ? `; the first ${formatInt(shown.length)} are shown, corrections first. Copy or download the rest.` : "."} Click a row to open that cell.</p>` +
    `<div class="table-wrap"><table class="corrections"><thead><tr><th>Cell</th><th>Genome</th><th>Cluster</th><th>Type</th><th class="num">Identity</th><th class="num">Coverage</th></tr></thead><tbody>${body}</tbody></table></div>`;
}
function openCellFromDrawer(row) {
  const clusterIdx = Number(row.getAttribute("data-cell-cluster"));
  const genomeIdx = Number(row.getAttribute("data-cell-genome"));
  document.getElementById("closeDetails")?.click();
  selected = clusterIdx;
  selectedGenome = genomeIdx;
  selectedSyntenyReference = false;
  populateOptions();
  setInspectorTab("cluster");
  render(true);
  scrollToGenome(genomeIdx);
}
// The drawer is written after this script, so its controls are reached by
// delegation rather than looked up now.
document.addEventListener("click", event => {
  const row = event.target.closest("#correctionsTable tr[data-cell-cluster]");
  if (row) {
    openCellFromDrawer(row);
    return;
  }
  if (event.target.closest("#openDetails")) {
    requestAnimationFrame(buildCorrectionsPreview);
    return;
  }
  const copyButton = event.target.closest("#copyCorrections");
  if (copyButton) {
    const copy = window.recoreCopyText || (() => Promise.resolve(false));
    copy(correctionsTsv()).then(ok => {
      copyButton.textContent = ok ? "Copied" : "Copy failed: use Download";
      setTimeout(() => { copyButton.textContent = "Copy as TSV"; }, 1800);
    });
    return;
  }
  if (event.target.closest("#downloadCorrections")) {
    downloadBlob(new Blob([correctionsTsv()], {type: "text/tab-separated-values"}), "recore_changed_cells.tsv");
  }
});
document.addEventListener("keydown", event => {
  const row = event.target.closest && event.target.closest("#correctionsTable tr[data-cell-cluster]");
  if (row && (event.key === "Enter" || event.key === " ")) {
    event.preventDefault();
    openCellFromDrawer(row);
  }
});

/* ---- Keyboard cell cursor -------------------------------------------------
   The matrix is a canvas, so without this a keyboard could reach a row only by
   typing its genome into find. The cursor reuses the hover marks and readout,
   and says where it is through the live region. */
const liveRegion = document.getElementById("liveRegion");
function announce(text) {
  if (!liveRegion) return;
  liveRegion.textContent = "";
  // A fresh node each time, so the same text twice is still read out.
  requestAnimationFrame(() => { liveRegion.textContent = text; });
}
const cursor = {row: -1, col: -1};
function cellReadout(row, col) {
  const clusterIdx = originalClusterAtVisibleColumn(col);
  const genomeIdx = leafOrder[row];
  const cluster = data.clusters[clusterIdx];
  if (!cluster || genomeIdx === undefined) return "";
  const state = stateAt(clusterIdx, genomeIdx);
  const evidence = state === "R" || state === "K" ? cellEvidence(clusterIdx, genomeIdx)[0] : null;
  return `${data.genomes[genomeIdx]}, ${cluster.name}: ${labels[state]}${evidence ? `, ${evidenceSummary(evidence)}` : ""}`;
}
function moveCursor(row, col) {
  const rows = leafOrder.length;
  const cols = visibleClusters.length;
  if (!rows || !cols) return;
  cursor.row = Math.max(0, Math.min(rows - 1, row));
  cursor.col = Math.max(0, Math.min(cols - 1, col));
  const view = viewportSize(matrixViewport);
  const y = cursor.row * cellH;
  const x = cursor.col * cellW;
  if (y < matrixPane.scrollTop) matrixPane.scrollTop = y;
  else if (y + cellH > matrixPane.scrollTop + view.h) matrixPane.scrollTop = y + cellH - view.h;
  if (x < matrixPane.scrollLeft) matrixPane.scrollLeft = x;
  else if (x + cellW > matrixPane.scrollLeft + view.w) matrixPane.scrollLeft = x + cellW - view.w;
  const text = cellReadout(cursor.row, cursor.col);
  setHover(cursor.row, cursor.col, text);
  announce(text);
}
matrixPane.addEventListener("focus", () => {
  if (cursor.row >= 0) return;
  const row = selectedGenome >= 0 ? rowOfGenome[selectedGenome] : Math.floor(matrixPane.scrollTop / cellH);
  moveCursor(Math.max(0, row), Math.max(0, selectedVisibleColumn()));
});
matrixPane.addEventListener("blur", () => setHover(-1, -1, defaultStatus));
matrixPane.addEventListener("keydown", event => {
  if (event.metaKey || event.altKey) return;
  const view = viewportSize(matrixViewport);
  const pageRows = Math.max(1, Math.floor(view.h / cellH) - 1);
  const pageCols = Math.max(1, Math.floor(view.w / cellW) - 1);
  const {row, col} = cursor.row >= 0 ? cursor : {row: 0, col: Math.max(0, selectedVisibleColumn())};
  let next = null;
  switch (event.key) {
    case "ArrowDown": next = [row + (event.shiftKey ? pageRows : 1), col]; break;
    case "ArrowUp": next = [row - (event.shiftKey ? pageRows : 1), col]; break;
    case "ArrowRight": next = [row, col + (event.shiftKey ? pageCols : 1)]; break;
    case "ArrowLeft": next = [row, col - (event.shiftKey ? pageCols : 1)]; break;
    case "PageDown": next = [row + pageRows, col]; break;
    case "PageUp": next = [row - pageRows, col]; break;
    case "Home": next = event.ctrlKey ? [row, 0] : [0, col]; break;
    case "End": next = event.ctrlKey ? [row, visibleClusters.length - 1] : [leafOrder.length - 1, col]; break;
    case "Enter":
    case " ":
      event.preventDefault();
      event.stopPropagation();
      selectedGenome = leafOrder[row];
      selectVisibleColumn(col);
      announce(`Opened ${cellReadout(row, col)}`);
      return;
    default:
      return;
  }
  event.preventDefault();
  event.stopPropagation();
  moveCursor(next[0], next[1]);
});
matrixPane.addEventListener("mouseleave", () => {
  tooltip.style.display = "none";
  tooltip.dataset.key = "";
  setHover(-1, -1, defaultStatus);
});
matrixPane.addEventListener("click", event => {
  // A cell is a cluster and a genome. Selecting only the cluster left an open
  // genome card describing some other genome, in some other cluster.
  const c = columnAtPointer(event, matrixPane);
  const r = rowAtPointer(event, matrixPane);
  if (r >= 0) selectedGenome = leafOrder[r];
  if (c >= 0) selectVisibleColumn(c);
  else if (r >= 0) renderGenomeDetails();
});

/* ---- Print and export ----------------------------------------------------
   Canvases are bitmaps sized to the on-screen viewport, so printing without
   redrawing them clips the panes to whatever the screen happened to be. Redraw
   once the print layout is in place, and again when it is torn down. */
/* Printing swaps the live panes for the figure: the whole matrix, in the
   light palette, with tips on their rows, a key and a caption. Redrawing the
   canvases for paper ran under screen styles (dark, the pane's screen width)
   and printed 17 of 24 columns. */
function preparePrint() {
  const main = document.getElementById("explorer");
  if (!main || document.getElementById("printFigure")) return;
  const image = document.createElement("img");
  image.id = "printFigure";
  image.className = "print-figure";
  image.alt = "Corrected matrix with the phylogeny";
  image.src = figureCanvas(true).toDataURL("image/png");
  main.before(image);
}
function finishPrint() {
  document.getElementById("printFigure")?.remove();
}
window.addEventListener("beforeprint", preparePrint);
window.addEventListener("afterprint", finishPrint);
if (window.matchMedia) {
  const printQuery = window.matchMedia("print");
  const onPrintChange = event => (event.matches ? preparePrint() : finishPrint());
  if (printQuery.addEventListener) printQuery.addEventListener("change", onPrintChange);
}

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
/* ---- Figure export ---------------------------------------------------------
   One bitmap with the tree, the lineage strip, the column names and the matrix
   drawn on one geometry, so tip N sits on row N. It used to paste the screen's
   canvases, stretching the tree to the matrix height (rows drifted up to 38 px
   off their tips), baking in the hover band, the selection tint and the dark
   theme, and upscaling a 1x bitmap. The whole matrix goes in when it fits;
   otherwise the window on screen, and the caption names the ranges. */
const FIGURE_SCALE = 3;
const FIGURE_MAX_DIM = 16000;
const FIGURE_MAX_PIXELS = 90e6;
const FIGURE_WHOLE_MAX = 4000;
function figureRange(forceWhole = false) {
  const rows = leafOrder.length;
  const cols = visibleClusters.length;
  if (cols * cellW <= FIGURE_WHOLE_MAX && rows * cellH <= FIGURE_WHOLE_MAX) {
    return {c0: 0, c1: cols, r0: 0, r1: rows, whole: true, cw: cellW, ch: cellH};
  }
  if (forceWhole) {
    // Paper takes the whole run: shrink the cells (to density if need be).
    return {
      c0: 0, c1: cols, r0: 0, r1: rows, whole: true,
      cw: Math.min(cellW, 1600 / Math.max(1, cols)),
      ch: Math.min(cellH, 1100 / Math.max(1, rows)),
    };
  }
  const view = viewportSize(matrixViewport);
  const left = matrixScrollLeft();
  const top = matrixScrollTop();
  return {
    c0: Math.max(0, Math.floor(left / cellW)),
    c1: Math.min(cols, Math.ceil((left + view.w) / cellW)),
    r0: Math.max(0, Math.floor(top / cellH)),
    r1: Math.min(rows, Math.ceil((top + view.h) / cellH)),
    whole: false,
    cw: cellW,
    ch: cellH,
  };
}
function drawFigure(forceWhole = false) {
  const {c0, c1, r0, r1, whole, cw, ch} = figureRange(forceWhole);
  const margin = 16;
  const gap = 8;
  const treeW = treeWidth;
  const stripW = colorMode.value === "lineage" ? LINEAGE_STRIP_W : 0;
  const band = columnLabelBand();
  const headH = band ? band + 10 : 10;
  const matrixW = Math.max(1, (c1 - c0) * cw);
  const matrixH = Math.max(1, (r1 - r0) * ch);
  const keyH = 22;
  const captionH = 34;
  const width = margin * 2 + treeW + stripW + gap + matrixW;
  const height = margin * 2 + headH + matrixH + keyH + captionH;
  const scale = Math.max(1, Math.min(
    FIGURE_SCALE,
    FIGURE_MAX_DIM / Math.max(width, height),
    Math.sqrt(FIGURE_MAX_PIXELS / (width * height)),
  ));
  const out = document.createElement("canvas");
  out.width = Math.round(width * scale);
  out.height = Math.round(height * scale);
  const f = out.getContext("2d");
  f.setTransform(scale, 0, 0, scale, 0, 0);
  f.fillStyle = panelBg;
  f.fillRect(0, 0, width, height);
  const bodyTop = margin + headH;
  const matrixLeft = margin + treeW + stripW + gap;

  // Tree: the screen's drawing code on the figure's context and rows.
  const keepCtx = treeCtx;
  const keepHighlights = highlightedRows;
  const keepCarriers = selectedSyntenyCarriers;
  const keepCellH = cellH;
  highlightedRows = new Set();
  selectedSyntenyCarriers = new Set();
  try {
    // The tree's rows follow the figure's row height, which paper may shrink.
    cellH = ch;
    treeCtx = f;
    f.save();
    f.beginPath();
    f.rect(margin, bodyTop, treeW, matrixH);
    f.clip();
    f.translate(margin, bodyTop);
    f.lineWidth = 1;
    drawTreeNode(treeRoot, r0 * ch, r0, r1, matrixH, labelToggle.checked && ch >= TIP_LABEL_MIN_ROW);
    f.restore();
  } finally {
    cellH = keepCellH;
    treeCtx = keepCtx;
    highlightedRows = keepHighlights;
    selectedSyntenyCarriers = keepCarriers;
  }
  if (stripW) {
    for (let r = r0; r < r1; r += 1) {
      f.fillStyle = lineageColor(lineageOfGenomeName(data.genomes[leafOrder[r]] || ""));
      f.fillRect(margin + treeW, bodyTop + (r - r0) * ch, stripW, Math.max(1, ch));
    }
  }

  // Matrix and column names.
  f.save();
  f.translate(matrixLeft, bodyTop);
  const figureCells = {firstCol: c0, lastCol: c1, firstRow: r0, lastRow: r1, left: c0 * cw, top: r0 * ch, cw, ch, dataW: matrixW, dataH: matrixH};
  if (inDensityMode(cw, ch)) {
    paintDensity(f, {...figureCells, scale, dx: Math.round(matrixLeft * scale), dy: Math.round(bodyTop * scale)});
  } else {
    paintCells(f, figureCells);
  }
  f.restore();
  if (band) {
    f.save();
    f.font = LABEL_FONT;
    f.fillStyle = labelInk;
    f.textBaseline = "middle";
    const flat = columnLabelsFlat();
    f.textAlign = flat ? "center" : "left";
    for (let c = c0; c < c1; c += 1) {
      const cluster = data.clusters[originalClusterAtVisibleColumn(c)];
      if (!cluster) continue;
      const x = matrixLeft + (c - c0) * cw + cw / 2;
      if (flat) {
        f.fillText(fitLabel(cluster.name, cw - 8), x, margin + band / 2);
      } else {
        f.save();
        f.translate(x, margin + band - 4);
        f.rotate(-Math.PI / 2);
        f.fillText(fitLabel(cluster.name, band - 12), 0, 0);
        f.restore();
      }
    }
    f.restore();
  }

  // Key: the states drawn, then the tree's scale.
  let x = margin;
  const keyY = bodyTop + matrixH + 8;
  f.font = LABEL_FONT;
  f.textBaseline = "middle";
  f.textAlign = "left";
  for (const state of shownStates()) {
    f.fillStyle = colors[state];
    f.fillRect(x, keyY, 10, 10);
    if (state === "A") {
      f.strokeStyle = strokes.A;
      f.strokeRect(x + 0.5, keyY + 0.5, 9, 9);
    }
    f.fillStyle = labelInk;
    f.fillText(labels[state], x + 14, keyY + 5);
    x += 14 + f.measureText(labels[state]).width + 14;
  }
  if (treeMode === "phylogram" && maxDistance > 0) {
    const step = niceStep(maxDistance * 0.32);
    const barW = step * (treeUsableWidth() / maxDistance);
    if (barW >= 14 && Number.isFinite(barW)) {
      x += 10;
      f.fillStyle = labelInk;
      f.fillRect(x, keyY + 5, barW, 1.2);
      f.fillText(`${formatDistance(step)} substitutions/site`, x + barW + 6, keyY + 5);
    }
  }

  // Caption: what is shown, out of what, under which settings.
  const policy = (document.querySelector(".readout b.tag") || {}).textContent || "";
  const version = (document.querySelector(".version-chip") || {}).textContent || "";
  const shown = whole
    ? `${plural(r1 - r0, "genome")} x ${plural(c1 - c0, "cluster")}`
    : `rows ${formatInt(r0 + 1)}-${formatInt(r1)} of ${formatInt(leafOrder.length)} genomes, columns ${formatInt(c0 + 1)}-${formatInt(c1)} of ${formatInt(visibleClusters.length)} clusters`;
  const caption = [
    `Recore corrected matrix: ${shown}`,
    `${treeCaption()}`,
    policy ? `${policy} policy` : "",
    version ? `recore ${version}` : "",
  ].filter(Boolean).join(" · ");
  f.fillStyle = labelInk;
  f.fillText(caption, margin, keyY + 26);
  return out;
}
/** A figure in the light palette, restored however the drawing ends. */
function figureCanvas(forceWhole = false) {
  const light = lightThemeDeclaration();
  try {
    drawingFigure = true;
    if (light) refreshThemeColors(light);
    return drawFigure(forceWhole);
  } finally {
    drawingFigure = false;
    if (light) refreshThemeColors();
  }
}
function exportViewPng() {
  const out = figureCanvas();
  return new Promise(resolve => {
    out.toBlob(blob => {
      if (blob) downloadBlob(blob, "recore_matrix.png");
      resolve(Boolean(blob));
    }, "image/png");
  });
}
/* Synteny is already vector, so export it as vector: gene arrows have to stay
   crisp at print resolution. */
/* Every drawn property is copied from the live element as a plain attribute:
   PDF tools and illustration software do not resolve var(), and an embedded
   <style> with the page's tokens drew black bars in MuPDF. rgba() becomes a hex
   colour plus an opacity, which every renderer understands. */
const SVG_STYLE_PROPS = ["fill", "fill-opacity", "stroke", "stroke-opacity", "stroke-width", "stroke-dasharray", "opacity", "font-family", "font-size", "font-weight", "font-style", "text-anchor", "dominant-baseline"];
function plainColor(value) {
  const match = /^rgba?\(([^)]+)\)$/.exec(String(value).trim());
  if (!match) return [value, null];
  const parts = match[1].split(/[ ,/]+/).filter(Boolean).map(Number);
  const hex = "#" + parts.slice(0, 3).map(v => Math.round(v).toString(16).padStart(2, "0")).join("");
  return [hex, parts.length > 3 ? parts[3] : null];
}
function exportSyntenySvg() {
  if (!syntenySvg || syntenyHidden) return false;
  const clone = syntenySvg.cloneNode(true);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  const live = [syntenySvg, ...syntenySvg.querySelectorAll("*")];
  const copies = [clone, ...clone.querySelectorAll("*")];
  const light = lightThemeDeclaration();
  // Read the styles under the light theme. The switch and the switch back
  // happen inside this one task, so nothing is painted in between, and
  // transitions are held off so none starts from the brief change.
  const root = document.documentElement;
  const theme = root.getAttribute("data-theme");
  root.classList.add("theme-frozen");
  if (theme !== "light") root.setAttribute("data-theme", "light");
  try {
  live.forEach((element, i) => {
    const target = copies[i];
    const computed = getComputedStyle(element);
    for (const prop of SVG_STYLE_PROPS) {
      const value = computed.getPropertyValue(prop);
      if (!value || value === "none" && prop !== "fill" && prop !== "stroke") continue;
      target.setAttribute(prop, value);
    }
    // A colour's own alpha folds into its opacity attribute; a transparent
    // label box written as plain #000000 drew black.
    for (const prop of ["fill", "stroke"]) {
      const [hex, alpha] = plainColor(computed.getPropertyValue(prop));
      if (!target.hasAttribute(prop) || hex === computed.getPropertyValue(prop)) continue;
      target.setAttribute(prop, alpha === 0 ? "none" : hex);
      if (alpha !== null && alpha > 0 && alpha < 1) {
        const own = Number(computed.getPropertyValue(`${prop}-opacity`)) || 1;
        target.setAttribute(`${prop}-opacity`, String(Math.round(own * alpha * 1000) / 1000));
      }
    }
    target.removeAttribute("class");
    target.removeAttribute("style");
  });
  } finally {
    if (theme === null) root.removeAttribute("data-theme");
    else if (theme !== "light") root.setAttribute("data-theme", theme);
    getComputedStyle(root).color;
    root.classList.remove("theme-frozen");
  }
  clone.querySelectorAll("style").forEach(node => node.remove());
  // A figure is read on paper: a white ground and a title naming the cluster.
  const width = Number(syntenySvg.getAttribute("width")) || syntenySvg.clientWidth;
  const height = Number(syntenySvg.getAttribute("height")) || syntenySvg.clientHeight;
  const ns = "http://www.w3.org/2000/svg";
  const group = document.createElementNS(ns, "g");
  group.setAttribute("transform", "translate(0 28)");
  while (clone.firstChild) group.appendChild(clone.firstChild);
  const ground = document.createElementNS(ns, "rect");
  ground.setAttribute("width", String(width));
  ground.setAttribute("height", String(height + 28));
  ground.setAttribute("fill", readThemeVar("--panel", "#ffffff", light));
  const title = document.createElementNS(ns, "text");
  title.setAttribute("x", "14");
  title.setAttribute("y", "19");
  title.setAttribute("font-family", "Inter, Helvetica, Arial, sans-serif");
  title.setAttribute("font-size", "13");
  title.setAttribute("font-weight", "600");
  title.setAttribute("fill", readThemeVar("--ink", "#18181b", light));
  title.textContent = `${(data.clusters[selected] || {}).name || ""}: unique synteny blocks (Recore ${(document.querySelector(".version-chip") || {}).textContent || ""})`;
  clone.appendChild(ground);
  clone.appendChild(title);
  clone.appendChild(group);
  clone.setAttribute("width", String(width));
  clone.setAttribute("height", String(height + 28));
  clone.setAttribute("viewBox", `0 0 ${width} ${height + 28}`);
  const blob = new Blob([new XMLSerializer().serializeToString(clone)], {
    type: "image/svg+xml",
  });
  downloadBlob(blob, "recore_synteny.svg");
  return true;
}
if (exportButton) {
  exportButton.addEventListener("click", async () => {
    const files = [];
    if (await exportViewPng()) files.push("recore_matrix.png");
    // Safari starts only one download per user gesture tick; the second has to
    // come a moment later, or it is dropped while the status claims it was saved.
    if (!syntenyHidden && syntenySvg) {
      await new Promise(resolve => setTimeout(resolve, 400));
      if (exportSyntenySvg()) files.push("recore_synteny.svg");
    }
    status.textContent = files.length ? `Exported ${files.join(" and ")}.` : "Nothing was exported.";
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
if (details) {
  details.addEventListener("click", event => {
    if (event.target.closest("[data-clear-cell]")) {
      selectedGenome = -1;
      highlightedRows = new Set();
      render(false);
      return;
    }
    const cell = event.target.closest("[data-cell-genome]");
    if (cell) {
      // A row of the evidence list opens that cell and brings its tip into view.
      selectedGenome = Number(cell.getAttribute("data-cell-genome"));
      scrollToGenome(selectedGenome);
      return;
    }
    const findTip = event.target.closest("[data-find-tip]");
    if (findTip) scrollToGenome(Number(findTip.getAttribute("data-find-tip")));
  });
}
if (genomeTab) genomeTab.addEventListener("click", () => setInspectorTab("genome"));

/* In a narrow window the inspector is a bottom sheet: a selection opens it. */
function openInspectorSheet() {
  if (inspectorPane && narrowQuery && narrowQuery.matches) inspectorPane.classList.add("sheet-open");
}
function closeInspectorSheet() {
  if (inspectorPane) inspectorPane.classList.remove("sheet-open");
}
document.getElementById("closeInspector")?.addEventListener("click", closeInspectorSheet);
function inspectGenome(genomeIdx, {focus = true} = {}) {
  if (genomeIdx < 0) return;
  selectedGenome = genomeIdx;
  renderGenomeDetails();
  if (focus) setInspectorTab("genome");
  openInspectorSheet();
}

function genomeQc(genomeIdx) {
  const row = data.genomeQc && data.genomeQc[genomeIdx];
  return row && row[0] ? {status: row[0], flags: row[1] || ""} : null;
}
function genomeQcHtml(genomeIdx) {
  const qc = genomeQc(genomeIdx);
  if (!qc) return "";
  const flags = qc.flags ? qc.flags.split(";").filter(Boolean).map(flag => flag.replace(/_/g, " ")) : [];
  return `<p class="qc-line"><span class="qc-chip ${escapeHtml(qc.status)}">QC ${escapeHtml(qc.status)}</span>${flags.length ? ` <small class="muted">${escapeHtml(flags.join("; "))}</small>` : ""}</p>`;
}
function renderGenomeDetails() {
  if (!genomeDetails) return;
  if (selectedGenome < 0) {
    genomeDetails.innerHTML = "<p class='muted'>Click a tip in the phylogeny to inspect that genome.</p>";
    return;
  }
  const name = data.genomes[selectedGenome] || "";
  const totals = genomeStats[selectedGenome] || {P: 0, R: 0, N: 0, K: 0, A: 0};
  const assessed = totals.P + totals.R + totals.N + totals.K + totals.A;
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
    <p class="muted">Tip ${row >= 0 ? formatInt(row + 1) : "?"} of ${formatInt(leafOrder.length)}${lineage && lineage !== "unknown" ? ` &middot; lineage ${escapeHtml(lineage)}` : ""}</p>
    ${genomeQcHtml(selectedGenome)}
    <dl class="cluster-dl">
      <dt>In the selected cluster</dt>
      <dd>${escapeHtml(cluster ? cluster.name : "none")}: <strong>${escapeHtml(labels[state])}</strong>${(() => {
        const row = cluster ? cellEvidence(selected, selectedGenome)[0] : null;
        return row ? `<br><span class="muted">${escapeHtml(evidenceSummary(row))}</span>` : "";
      })()}</dd>
    </dl>
    <div class="card-title">Cluster content across ${plural(assessed, "cluster")}</div>
    <div class="count-rows">${stateCountRows(totals, assessed)}</div>
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

/* The run-details drawer ships with every report, the viewer only when there is
   a tree. It asks by name and we answer if the cluster is here. */
window.addEventListener("recore:select-cluster", event => {
  const name = event.detail && event.detail.cluster;
  if (!name) return;
  if (selectMatrixClusterByName(name)) {
    setInspectorTab("cluster");
    event.detail.handled = true;
  }
});

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
/* Lower is better. An exact name, with its case, beats the same name in
   another case: Enter on "Group_1" opened group_1 when both existed. */
function scoreHit(name, note, query, raw = query) {
  const lowName = name.toLowerCase();
  if (!query) return 5;
  if (name === raw) return 0;
  if (lowName === query) return 1;
  if (lowName.startsWith(query)) return 2;
  if (lowName.includes(query)) return 3;
  if (note && note.toLowerCase().includes(query)) return 4;
  return -1;
}

function buildPaletteRows() {
  const raw = search.value.trim();
  const query = raw.toLowerCase();
  const hits = [];
  if (paletteScope !== "genome") {
    for (let i = 0; i < visibleClusters.length; i += 1) {
      const idx = visibleClusters[i];
      const cluster = data.clusters[idx];
      const note = cluster.annotation || "";
      const score = scoreHit(cluster.name, note, query, raw);
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
      const score = scoreHit(name, note, query, raw);
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
      ? `${formatInt(rows.length)} of ${formatInt(total)}`
      : plural(total, "result");
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
    if (inspectorPane && inspectorPane.classList.contains("sheet-open")) { closeInspectorSheet(); return; }
    if (shortcutSheet && !shortcutSheet.hidden) { setSheetOpen(false); return; }
    if (anyPaneMenuOpen()) { setPaneMenu(null, false); return; }
    return;
  }
  if (event.defaultPrevented || paletteIsOpen() || typingInField(event.target)) return;
  if (event.metaKey || event.ctrlKey || event.altKey) return;
  // Nothing acts behind an open dialog: the drawer handles its own keys, and
  // the shortcut sheet only answers to the key that closes it.
  if (document.getElementById("runDetails")?.classList.contains("open")) return;
  if (shortcutSheet && !shortcutSheet.hidden && event.key !== "?") return;
  // A focused seam uses the arrows to resize itself.
  if (event.target.closest && event.target.closest("[role='separator']")) return;

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
      fitZoomToData(true);
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
  const usable = treeUsableWidth();
  const pxPerUnit = usable / maxDistance;
  // Aim the bar at about a third of the drawing width and then round to a nice
  // number, rather than rounding a fixed fraction of the tree depth: with tip
  // labels on, a 320px pane leaves 138px to draw in, and maxDistance/4 rounded
  // down to a bar of 17px, one pixel under the old floor, so no bar at all.
  const step = niceStep(maxDistance * 0.32);
  const width = step * pxPerUnit;
  if (!(width >= 14) || !Number.isFinite(width)) {
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
      // One bar per device column: the mean prevalence of the clusters in it
      // (the maximum filled every column at 30,000 clusters and said nothing),
      // and a mark on top where any of them has evidence behind a change.
      for (let x = 0; x < width; x += 1) {
        const from = Math.floor(x / width * n);
        const to = Math.max(from + 1, Math.floor((x + 1) / width * n));
        let sum = 0;
        let count = 0;
        let changed = 0;
        for (let i = from; i < to && i < n; i += 1) {
          const idx = visibleClusters[i];
          sum += data.clusters[idx]?.present || 0;
          count += 1;
          if (evidenceCount(idx)) changed += 1;
        }
        const h = Math.max(1, Math.round(sum / Math.max(1, count) / genomes * (height - 6)));
        octx.globalAlpha = 0.55;
        octx.fillStyle = colors.P;
        octx.fillRect(x, height - 2 - h, 1, h);
        if (changed) {
          octx.globalAlpha = 0.35 + 0.65 * Math.min(1, changed / Math.max(1, count) * 4);
          octx.fillStyle = colors.R;
          octx.fillRect(x, 0, 1, 3);
        }
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

/* The genome card names the selected cluster, so it follows every change of
   selection, not only a click on a tip. */
let genomeCardKey = "";
afterRender = () => {
  const key = `${selected}|${selectedGenome}`;
  if (key === genomeCardKey) return;
  genomeCardKey = key;
  renderGenomeDetails();
};

renderGenomeDetails();
populateOptions();
render(false);
window.recoreExplorerRunning = true;

})();
