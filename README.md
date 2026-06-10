# vanWeperenV_2026_EpiduralRTX

Analysis code accompanying:

> van Weperen V.Y.H.\*, Hoang J.D.\*, Jani N., Avasthi S., Chan C.A., Cao K.,
> Lokhandwala Z.A., Emamimeybodi M., Atmani K., Vaseghi M.
> *Spinal nociceptive denervation impedes chronic autonomic remodeling after
> myocardial infarction.* (\*equal contribution)

This repository contains the R code used to perform the **downstream** Gene
Ontology (GO) pathway enrichment analysis of plasma differentially expressed
proteins (DEPs). The enrichment tables it produces were used to generate the
associated figure panel in GraphPad Prism.

> **Scope.** This repository covers the GO pathway enrichment **only**. Plasma
> sample preparation, LC-MS/MS acquisition, and identification of differentially
> expressed proteins (DIA-NN, FragPipe-Analyst) were performed by the UCLA
> Proteome Research Center, which provided the DEP table (`Data/dep.xlsx`) used
> here as input. The final figure was produced in GraphPad Prism from the
> exported enrichment tables. The raw and processed proteomics data are deposited
> in the ProteomeXchange Consortium via the PRIDE partner repository (accession
> **[TO BE ADDED]**, https://www.ebi.ac.uk/pride/), not in this repository.

---

## Overview

Plasma proteomes from infarcted animals treated with epidural resiniferatoxin
(cRTX+MI) versus vehicle (MI) were compared by data-independent acquisition mass
spectrometry. The DEP identification (upstream of this repository) was carried
out by the core facility using DIA-NN and FragPipe-Analyst. Starting from the
resulting DEP table, the code here tests GO term over-representation using
[RITAN](https://doi.org/10.7717/peerj.6994). The script:

1. Reads the core-facility DEP table and splits proteins into up- and
   down-regulated sets.
2. Runs term enrichment against the GO database for each set.
3. Writes the significant enrichment tables, which were then imported into
   GraphPad Prism to produce the published figure.

---

## Repository structure

```
vanWeperenV_2026_EpiduralRTX/
├── README.md
├── LICENSE
├── pathway_enrichment.R           # analysis script
└── Data/
    └── dep.xlsx                   # differentially expressed proteins (input)
```

Output files (`up_pathways.xlsx`, `down_pathways.xlsx`) are written to the
repository root on run.

---

## Input format

`Data/dep.xlsx` — one row per protein, **produced by the UCLA Proteome Research
Center** and exported from FragPipe-Analyst. Provided here as the input to the
downstream analysis. Required columns:

| Column | Description |
| --- | --- |
| `Protein ID` | UniProt accession |
| `Gene Name` | Gene symbol (used for GO mapping) |
| `X1_vs_X2_log2 fold change` | log2 fold change for the contrast |
| `X1_vs_X2_p.val` | nominal p-value |
| `X1_vs_X2_p.adj` | adjusted p-value |
| `significant`, `X1_vs_X2_significant`, `imputed`, `num_NAs` | FragPipe-Analyst flags |

`X1_vs_X2` denotes the cRTX+MI vs Vehicle+MI contrast.

---

## Requirements

- R (≥ 4.2 recommended)
- CRAN: `tidyverse`, `readxl`, `rstudioapi`, `WriteXLS`
- Bioconductor: `RITAN`, `RITANdata`

All packages are installed automatically on first run if absent (Bioconductor
packages via `BiocManager`).

---

## Usage

From RStudio, open `pathway_enrichment.R` and source it; the working directory is
set to the script location automatically.

From the command line:

```bash
Rscript pathway_enrichment.R
```

When run outside RStudio, ensure the working directory is the repository root so
that `Data/dep.xlsx` resolves correctly.

### Parameters

Adjust at the top of the script:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `logfc_thr` | `0.5` | \|log2 fold change\| cutoff defining a DEP |
| `p_thr` | `0.05` | p-value cutoff for DEPs and for enrichment |
| `resources` | `"GO"` | RITAN geneset resource |

---

## Output

| File | Contents |
| --- | --- |
| `up_pathways.xlsx` | Significant GO terms enriched among up-regulated DEPs |
| `down_pathways.xlsx` | Significant GO terms enriched among down-regulated DEPs |

Each table reports, per GO term, the gene ratio, protein count, and enrichment
p-value. These tables were imported into GraphPad Prism to produce the published
figure panel.

---

## Data availability

The raw and processed plasma mass spectrometry data are deposited in the
ProteomeXchange Consortium via the PRIDE partner repository under accession
**[TO BE ADDED]** (https://www.ebi.ac.uk/pride/). This repository contains the
analysis code and the DEP table used as its input only.

---

## Citation

If you use this code, please cite the associated publication (above) and RITAN:

> Zimmermann M.T., Kabat B., Grill D.E., Kennedy R.B., Poland G.A.
> RITAN: rapid integration of term annotation and network resources.
> *PeerJ* 7:e6994 (2019). https://doi.org/10.7717/peerj.6994

---

## License

Released under the MIT License. See [LICENSE](LICENSE).
