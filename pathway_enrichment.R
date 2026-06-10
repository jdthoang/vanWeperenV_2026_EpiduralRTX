################################################################################
## Gene Ontology pathway enrichment of plasma differentially expressed proteins
##
## Companion code for:
##   van Weperen et al. "Spinal nociceptive denervation impedes chronic autonomic
##   remodeling after myocardial infarction"
##
## Description
##   Performs over-representation (term enrichment) analysis of differentially
##   expressed plasma proteins (DEPs) against the Gene Ontology (GO) database
##   using RITAN, separately for up- and down-regulated proteins. For each set it
##   reports, per GO term, the gene ratio (proteins in the DEP list mapping to the
##   term / term set size), the protein count, and the enrichment p-value. The
##   exported enrichment tables were used to generate the published figure panel
##   in GraphPad Prism.
##
## Scope
##   This script covers the DOWNSTREAM pathway enrichment only. Plasma sample
##   preparation, LC-MS/MS acquisition, and DEP identification (DIA-NN,
##   FragPipe-Analyst) were performed by the UCLA Proteome Research Center, which
##   provided the DEP table (Data/dep.xlsx) used here as input. The final figure
##   was produced in GraphPad Prism from the enrichment tables exported below.
##   Raw and processed proteomics data are in ProteomeXchange/PRIDE under
##   accession [TO BE ADDED], not in this repository.
##
## Input
##   Data/dep.xlsx  -- one row per protein, PRODUCED BY THE CORE FACILITY, with
##     columns:
##     Protein ID, Gene Name,
##     X1_vs_X2_log2 fold change, X1_vs_X2_p.val, X1_vs_X2_p.adj,
##     significant, X1_vs_X2_significant, imputed, num_NAs
##   ("X1_vs_X2" denotes the contrast as exported by FragPipe-Analyst /
##   DIA-NN label-free quantification; here cRTX+MI vs Vehicle+MI.)
##
## Output
##   up_pathways.xlsx, down_pathways.xlsx -- full significant enrichment tables
##                                           (subsequently plotted in Prism)
##
## Author: Jonathan D. Hoang
## License: MIT (see LICENSE)
################################################################################


# ---- 1. Dependencies ---------------------------------------------------------
# CRAN packages install via install.packages(); Bioconductor packages (RITAN,
# RITANdata) install via BiocManager. BiocManager itself is bootstrapped if
# absent.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_pkgs <- c("tidyverse", "readxl", "rstudioapi", "WriteXLS")
bioc_pkgs <- c("RITAN", "RITANdata")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE)
}

library(RITAN)
library(RITANdata)   # provides geneset_list (GO, etc.)
library(tidyverse)
library(readxl)


# ---- 2. Parameters -----------------------------------------------------------
# Edit these to re-run with different thresholds or input/output locations.

dep_file  <- file.path("Data", "dep.xlsx")  # differentially expressed proteins
logfc_thr <- 0.5    # |log2 fold change| cutoff for a DEP
p_thr     <- 0.05   # nominal p-value cutoff for a DEP and for enrichment
resources <- "GO"   # RITAN geneset resource(s); manuscript uses Gene Ontology

# Resolve the working directory to this script's location (RStudio) so that the
# relative paths above are stable regardless of where R was launched. Falls back
# to the current working directory when run via Rscript / a non-RStudio session.
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}


# ---- 3. Load and split differentially expressed proteins ---------------------
# DEPs are defined as proteins passing both the fold-change and p-value cutoffs.
# They are split by the sign of the log2 fold change into up- and down-regulated
# sets (relative to the cRTX+MI vs Vehicle+MI contrast). Gene symbols are upper-
# cased to match RITAN's GO geneset identifiers.

dep <- read_xlsx(dep_file)

# Standardize the two columns used downstream to syntactic names.
dep <- dep %>%
  rename(log2fc = `X1_vs_X2_log2 fold change`,
         pval   = `X1_vs_X2_p.val`)

is_dep   <- abs(dep$log2fc) > logfc_thr & dep$pval < p_thr
up_dep   <- toupper(dep$`Gene Name`[is_dep & dep$log2fc > 0])
down_dep <- toupper(dep$`Gene Name`[is_dep & dep$log2fc < 0])

message(sprintf("Up-regulated DEPs: %d | Down-regulated DEPs: %d",
                length(up_dep), length(down_dep)))


# ---- 4. Term enrichment (over-representation) against GO ----------------------
# RITAN::term_enrichment tests each set against every term in `resources`.
# GeneRatio = n (DEPs in term) / n.set (term size). Significant terms (p < p_thr)
# are ordered by p-value and written out; these tables were imported into
# GraphPad Prism to produce the published figure.

enrich <- function(genes) {
  e <- term_enrichment(genes, resources = resources)
  e$GeneRatio <- e$n / e$n.set
  e <- e[e$p < p_thr, ]
  e[order(e$p, decreasing = FALSE), ]
}

e_up   <- enrich(up_dep)
e_down <- enrich(down_dep)

WriteXLS::WriteXLS(e_up,   "up_pathways.xlsx")
WriteXLS::WriteXLS(e_down, "down_pathways.xlsx")

message("Done. Wrote: up_pathways.xlsx, down_pathways.xlsx")
