# Evidence states

The multistate matrix, `gene_presence_absence_recore_states.csv`, has the same
shape as the input presence/absence matrix but stores a biological state per cell
instead of a 0 or a 1. It is the file to read for audit and for figures; the
binary matrix is the file to feed downstream tools.

| State | Meaning | Corrected in the binary matrix under `balanced` |
|-------|---------|--------------------------------------------------|
| `PRESENT_INPUT` | present in the original matrix | already present |
| `ABSENT` | searched, and no accepted evidence found | no |
| `NOT_SEARCHED_LOW_PREVALENCE` | never searched: the cluster is in fewer genomes than `--min-prevalence` | no |
| `NOT_SEARCHED_NO_CLUSTER_SEQUENCE` | never searched: the cluster has no sequence in the pan-genome reference | no |
| `NOT_SEARCHED_NO_FASTA` | never searched: the genome has no FASTA file | no |
| `STRICT_PRESENT` | same-locus evidence from annotation or synteny | yes |
| `ALIAS_CLUSTER` | the same biological locus assigned to another cluster | yes |
| `SEQUENCE_ONLY` | full-length sequence evidence, without stronger annotation context | yes |
| `MISSING_CDS_UNVALIDATED` | full-length sequence in a locus with no CDS annotation, and no synteny support | no |
| `PARTIAL_HOMOLOGY` | partial homologous sequence | no |
| `PSEUDOGENE` | disrupted gene evidence | no |
| `HOMOLOG_PARALOG` | a homologous or paralogous hit | no |
| `REJECTED_CONTEXT` | candidate rejected by context or paralog checks | no |

The last five rows are the point of the tool. Under `balanced`, none of them
becomes a presence call, and all of them stay visible with the reason attached.

## Absent or not searched

`ABSENT` is a finding: Recore searched that gap and nothing cleared the
thresholds. A gap outside the search is not, and says why instead. Under the
default `--min-prevalence` of 0.1 that is every gap of a cluster in fewer than
one genome in ten, which in most pangenomes is most of the empty cells. Before
this distinction both read `ABSENT`.

When two reasons apply, the cluster's is written, in the order the pipeline
filters gaps: prevalence, then a missing cluster sequence, then a missing FASTA.
That is the order of the gap lines in `summary.txt`, which count the same cells
except those NSC later filled or merged away.

## Reading the states in the report

The HTML viewer collapses these thirteen states into the six a matrix can
actually be read in:

| Cell | States |
|------|--------|
| input | `PRESENT_INPUT` |
| corrected | the states corrected by the active policy |
| NSC | presence merged in from a nucleotide-redundant cluster |
| kept out | evidence the policy, a quality floor or a guard kept out of the matrix |
| absent | `ABSENT` and the other states that are not corrected |
| not searched, hatched | the three `NOT_SEARCHED_*` states |

Switching the matrix to mechanism colouring instead shows the evolutionary call
from `evolutionary_events.tsv`, which answers a different question: not "is it
there?" but "what kind of event is this?".

## Which policy corrects which state

See [Matrix policies](../guide/matrix-policies.md) for the full table. In short:
`strict` corrects only `STRICT_PRESENT` and `ALIAS_CLUSTER`; `balanced` adds
`SEQUENCE_ONLY`; `sensitive` also admits partial and pseudogene evidence that
clears the adaptive guard.
