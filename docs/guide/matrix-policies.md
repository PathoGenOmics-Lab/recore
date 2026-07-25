# Matrix policies

Recore always records every piece of detected evidence in the audit files. The
matrix policy controls one thing only: which accepted evidence is written into
the **binary** corrected matrix.

Nothing is ever thrown away by a policy. A cell that `strict` declines to correct
is still in `gene_presence_absence_recore_states.csv` with its state, and still
in `orthology_evidence.tsv` with the reason.

## Choosing one

| Goal | Policy | What it corrects | What stays audit-only |
|------|--------|------------------|------------------------|
| Paper-ready binary PA matrix | `balanced` | strict same-locus calls, alias clusters, and full-length sequence evidence | pseudogenes, partial homology, context-rejected hits |
| Maximum conservatism | `strict` | only strict same-locus and alias-cluster evidence | full sequence-only hits without annotation or synteny support |
| Exploratory rescue | `sensitive` | accepted strict, alias, sequence-only, partial and pseudogene evidence that passes the quality guard | rejected contexts, weak accepted evidence, and homologous or paralogous traps |

`balanced` is the recommendation for anything that will be published. It rescues
the false absences that motivated the tool without letting a disrupted homologue
become a clean presence call.

Reach for `strict` when the downstream analysis is sensitive to false positives
in a way that outweighs the false absences, for example a core-genome definition
that feeds a phylogeny.

Reach for `sensitive` when you are exploring, and read the multistate matrix
alongside it rather than the binary one.

## The adaptive orthology guard

Policy is not the only gate. Every accepted rescue is also scored by an adaptive
guard that combines:

- sequence identity and query coverage;
- locus and synteny evidence;
- ORF integrity;
- the evidence state;
- cluster prevalence;
- representative sequence length.

Rare clusters and very short clusters have to clear a **higher** sequence bar,
because they are the easiest to over-call from weak homology. Same-locus evidence
such as a passing synteny context can override that stricter sequence floor for
`STRICT_PRESENT` and `ALIAS_CLUSTER` calls, but weak sequence-only evidence is
held back as `audit_only_quality`.

The guard shows its work in `orthology_evidence.tsv`:

| Column | Meaning |
|--------|---------|
| `adaptive_identity_floor` | the identity this cluster had to clear |
| `adaptive_coverage_floor` | the coverage this cluster had to clear |
| `adaptive_sequence_pass` | whether the sequence evidence cleared both |
| `orthology_confidence_score` | combined score |
| `orthology_confidence_tier` | the score bucketed for reading |
| `risk_flags` | what made this call risky |

So a cell that was not corrected under `balanced` can be traced to either the
policy or the guard, and the two are distinguishable.

## Reporting it

The wording that matches what `balanced` actually does:

> We used Recore with the balanced matrix policy, which rescues full-length false
> absences while retaining partial or disrupted homologues as multistate
> audit-only evidence.

If a specific cluster needs its own sentence, `orthology_evidence.tsv` carries
per-cluster wording guidance, and the run report shows it in the inspector under
**recommended wording**.
