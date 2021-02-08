
2020-01-12

# External validation re-run with TCGA and METABRIC

## Handling of multi-transcript probes


> **What transcripts are used for annotating mutations?**
> 
> Prior to loading a study into cBioPortal, we run all mutation data through a standard pipeline (see above), which re-annotates all mutations to the canonical UniProt transcript,

[Standard cBioportal probe list](https://github.com/mskcc/vcf2maf/blob/main/data/isoform_overrides_uniprot)

Imported by

    wget https://github.com/mskcc/vcf2maf/raw/main/data/isoform_overrides_uniprot


***

## Differential gene expression analysis for METABRIC

- additional filter of FC >1.5 / <0.5 required?


## Technical concordance between Aligent and RNA-seq

- Pearson correlation (coeff. > 0.7)
  - BH correction for multiple testing
  - FDR <0.01

### METABRIC uses log2 scaled intensities

External validation re-run:
- both datasets using log2
- convert to z-score shouldn't matter, but will log2 be better than converting all to z-scores?

***
