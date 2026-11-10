# CNmimicry

Code accompanying the MSc thesis *Systematically Predicting Regulatory Alterations Mimicking Copy Number Alterations* (Hebrew University of Jerusalem, 2026). The repository contains code only; data, results, figures and tables are not included.

CNmimicry identifies amplified-like diploid (ALD) and deleted-like diploid (DLD) expression states in TCGA tumors and examines their associations with tumor genotype, promoter methylation, somatic mutations, and overall survival. The workflow combines expression and copy-number data with GISTIC gene sets, generates cohort-level summaries, and provides scripts for downstream analyses.

## Requirements

- R and the packages listed in [software_versions.txt](software_versions.txt). The recorded environment uses R 4.4.0 on Windows. 
- Required packages are distributed through CRAN and Bioconductor. The version record documents the local installation; it does not install packages automatically.
- Internet access for retrieving uncached TCGA data through TCGAbiolinks and uncached Broad GDAC GISTIC resources.
- Sufficient disk space for downloaded data, prepared matrices, and analysis-specific caches. Methylation and somatic mutation downloads can be substantial.

Open `Workflow.Rproj`, or start R with the Workflow root as the working directory. Run all scripts from that directory; their `source()` calls use project-relative paths.

## Project Structure

| Path | Contents |
| --- | --- |
| `config.R` | Cohorts, core thresholds, GISTIC selection, CNV normalization, and project paths. |
| `run_all.R` | Core workflow entry point. |
| `make_excel.R` | ALD and DLD Excel summaries of core results. |
| `R/` | Core data preparation, GISTIC helpers, and ALD/DLD calculations. |
| `analysis/` | Downstream analysis scripts and helpers (caches and results are generated locally). |
| `cache/` | Prepared cohort expression/CNV matrices and parsed GISTIC gene sets. |
| `GDCdata/` | Core TCGA downloads. |
| `resources/` | GDAC GISTIC source files and biological annotation resources. |
| `results_cnvnorm_noGISTICthr/` | Normalized core outputs without the thresholded-GISTIC intersection. |
| `excel/` | Core summary workbooks. |
| `software_versions.txt` | Recorded R and directly used package versions. |

Directories such as `cache/`, `GDCdata/`, `resources/`, `results_cnvnorm_noGISTICthr/` and `excel/` are created locally when the workflow runs and are not tracked in this repository. Analysis-specific caches and downloads also reside under `analysis/`; these are distinct from the core cache and download directories.

## Primary Configuration

Review `config.R` before running. The primary configuration is:

| Setting | Value |
| --- | --- |
| Cohorts | 20 TCGA cohorts defined in `cancers`; BRCA is restricted to ER-positive tumors. |
| General expression filter | TPM >1 in >5% of cohort samples, applied during data preparation. |
| CNV normalization | Sample-wise median normalization enabled. |
| GISTIC significant-gene restriction | Enabled, using `wide_peak` with peak q-value <=0.25. |
| GDAC resource run | `analyses__2016_01_28`. |
| Thresholded-GISTIC intersection | Disabled. |
| ALD amplified-group median requirement | Strictly >1 TPM (`ald_min_amplified_median_tpm = 1`). |
| Copy-number alteration prevalence | Amplification or deletion in >3% of cohort samples (`prevalence_thresholds = c(0.03)`). |

The general expression filter and ALD amplified-group median requirement serve different purposes. The former filters expression data upstream; the latter determines which genes are eligible for ALD classification within a cohort. It does not restrict true amplification or change DLD classification.

The copy-number alteration prevalence threshold selects eligible amplified or deleted genes. It is distinct from the recurrent ALD/DLD gene-hit definition: mimicry calls in >20% of the gene's diploid samples.

Core processing iterates over all configured prevalence thresholds. Many downstream analyses use `prevalence_thresholds[1]`. Some analysis-specific settings, such as minimum group sizes and mutation categories, are defined in the corresponding scripts rather than `config.R`.

Setting `ald_min_amplified_median_tpm` to `NULL` disables the ALD median requirement. Generate matching core results before running downstream analyses with any different configuration. Changing settings does not regenerate existing results automatically.

Prepared expression/CNV caches are reused when present. Changes to upstream expression-filter settings require rebuilding the affected prepared cache; existing cache filenames do not encode those thresholds. Analysis-specific cache behavior is documented in the module READMEs and script headers.

## Running the Core Workflow

Rerunning a script with the same output settings can overwrite its existing outputs.

From the project root:

```r
source("run_all.R")
source("make_excel.R")
```

`run_all.R` prepares or reuses cohort data, selects eligible gene sets, computes ALD/DLD calls, and writes per-cohort outputs. `make_excel.R` reads the matching core summaries and writes ALD and DLD workbooks.

For the primary configuration, core files are written under:

```text
results_cnvnorm_noGISTICthr/with_SIG_wide_peak/aldTPMgt1/<cohort>/prev_0.03/
```

Each directory contains `summary.csv`, sample-level ALD/DLD summaries, long-format call files, amplified/deleted gene lists, the ALD-eligible gene list, and per-gene prevalence objects. The ALD score denominator uses genes eligible for ALD classification; true-amplification analyses retain the original amplified-gene set.

## Downstream Analyses

Run the required modules after generating their matching core inputs. `run_all.R` does not execute downstream analyses. The table below is a directory of modules, not a universal execution order; the module READMEs and script headers describe their dependencies.

| Module | Script or documentation |
| --- | --- |
| Permutation evaluation | `analysis/permutations/run_permutations.R` |
| IDH-associated burden and gene enrichment in LGG | `analysis/idh/run_IDH_LGG_ALD_DLD_analysis.R` |
| CTCF-associated burden and gene enrichment in UCEC | `analysis/ctcf/run_CTCF_UCEC_ALD_DLD_analysis.R` |
| Promoter methylation | [Module README](analysis/methylation/README.md) |
| Survival screen and follow-ups | [Module README](analysis/survival/README.md) |
| DLD mutation associations | `analysis/mutation/run_DLD_mutation_gene_association.R` |
| DLD mutation composition | `analysis/mutation/run_DLD_mutation_type_composition.R` |
| Subtype pilot | `analysis/subtype/run_subtype_ALD_DLD_pilot.R` |

Permutation evaluation uses 1,000 simulations per cohort and seed 123, as specified in its script.

The subtype pilot is restricted to predefined LGG, SARC, and UCEC subtype groups; select `pilot_cohort` in the script for each run. It requires the corresponding whole-cohort survival results.

## Sensitivity Analyses

Generate matching core results for each sensitivity configuration before running the corresponding analyses:

- Narrower GISTIC peak-region filtering: keep `use_gistic_significant = TRUE` and set `gistic_sig_gene_source = "genes_in_peak"`.
- No significant-peak restriction: set `use_gistic_significant = FALSE` to generate the `no_SIG` core results for the peak-restriction comparison, alongside the `wide_peak` and `genes_in_peak` results.
- Copy-number prevalence robustness: set `prevalence_thresholds = c(0.03, 0.05, 0.10)` to generate the three core result sets for the prevalence-threshold comparison.

Keep the other settings consistent with the intended comparison. Restore the primary configuration before rerunning the primary analyses.

## Scope of This Repository

This repository includes the core workflow and the downstream analysis modules. Figure and table builders, generated results and downloaded data are not included. Some survival follow-up scripts expect a significant-hit summary table that was produced by a figure builder in the full project; their headers and the [survival README](analysis/survival/README.md) describe the expected inputs.

## Output Naming

Output paths and filenames distinguish relevant analysis settings. The `aldTPMgt1` tag identifies core runs requiring an amplified-group median >1 TPM. DLD mutation and methylation analyses read DLD inputs from the configured core run but retain DLD output names without an ALD TPM tag, because their calculations are independent of that requirement.

