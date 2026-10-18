# Citation

If you use Recore, please cite:

> Ruiz-Rodriguez P, Coscolla M. (2026). Recore: post-hoc rescue and quality
> control of pangenome presence/absence matrices. *In preparation.*

```bibtex
@misc{recore,
  author       = {Ruiz-Rodriguez, Paula and Coscolla, Mireia},
  title        = {Recore: post-hoc rescue and quality control of pangenome
                  presence/absence matrices},
  year         = {2026},
  note         = {In preparation},
  howpublished = {\url{https://github.com/PathoGenOmics-Lab/recore}}
}
```

Please also cite the tools Recore depends on for the run you report:

- **minimap2**, for every sequence search:
  Li H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
  *Bioinformatics* 34(18):3094-3100.
- **MMseqs2**, if you used `--nsc-engine mmseqs2`:
  Steinegger M, Soeding J. (2017). MMseqs2 enables sensitive protein sequence
  searching for the analysis of massive data sets. *Nature Biotechnology*
  35:1026-1028.
- the pangenome tool that produced the input matrix, for example **Panaroo**:
  Tonkin-Hill G, et al. (2020). Producing polished prokaryotic pangenomes with
  the Panaroo pipeline. *Genome Biology* 21:180.

## Contact

Paula Ruiz-Rodriguez, paula.ruiz.rodriguez@csic.es

PathoGenOmics Lab, [github.com/PathoGenOmics-Lab](https://github.com/PathoGenOmics-Lab)

## License

GNU Affero General Public License v3.0 (`AGPL-3.0-only`). See
[LICENSE](https://github.com/PathoGenOmics-Lab/recore/blob/main/LICENSE).
