# Evidence states

The multistate matrix, `gene_presence_absence_recore_states.csv`, has the same
shape as the input presence/absence matrix but stores a biological state per cell
instead of a 0 or a 1. It is the file to read for audit and for figures; the
binary matrix is the file to feed downstream tools.

| State | Meaning | Corrected in the binary matrix under `balanced` |
|-------|---------|--------------------------------------------------|
| `PRESENT_INPUT` | present in the original matrix | already present |
| `ABSENT` | no accepted evidence found | no |
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

## Reading the states in the report

The HTML viewer collapses these ten states into the four colours a matrix can
actually be read in:

| Colour | States |
|--------|--------|
| input | `PRESENT_INPUT` |
| rescued | the states corrected by the active policy |
| NSC | presence merged in from a nucleotide-redundant cluster |
| absent | everything not corrected |

Switching the matrix to mechanism colouring instead shows the evolutionary call
from `evolutionary_events.tsv`, which answers a different question: not "is it
there?" but "what kind of event is this?".

## Which policy corrects which state

See [Matrix policies](../guide/matrix-policies.md) for the full table. In short:
`strict` corrects only `STRICT_PRESENT` and `ALIAS_CLUSTER`; `balanced` adds
`SEQUENCE_ONLY`; `sensitive` also admits partial and pseudogene evidence that
clears the adaptive guard.
