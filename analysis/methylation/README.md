# Promoter methylation analysis

This module tests whether ALD and DLD genes show promoter methylation differences in the same samples where the mimicry state is observed.

## Scope

- The analysis includes parallel ALD and DLD promoter methylation screens.
- For ALD, eligible amplified genes are tested in ALD tumors versus other diploid tumors for the same gene.
- When `ald_min_amplified_median_tpm` is enabled in `config.R`, the ALD screen reads the corresponding `aldTPMgt<value>` core subfolder and `ald_eligible_genes.rds`. ALD results, input summaries, and cross-cancer summaries receive the same TPM tag; `NULL` retains baseline paths and names.
- ALD input summaries distinguish the original amplified-gene count from the ALD-eligible count. Both screens read the core folder selected by the TPM setting in config.R. The filter does not change promoter matrix caches or DLD calls, gene lists, calculations, or output names.
- For DLD, eligible deleted genes are tested in DLD tumors versus other diploid tumors for the same gene.
- True amplified or deleted tumors are not used as the reference group.
- The statistical test is a Wilcoxon rank-sum test.
- P values are adjusted within each cohort and analysis mode using BH FDR.
- The analysis does not adjust for subtype, purity, stage, age, or other covariates.

## Current methylation settings

- Platform: Illumina Human Methylation 450
- Sample type: Primary Tumor
- Promoter regions: TSS200, TSS1500, 5'UTR
- Promoter methylation summary: median beta value across promoter-associated probes for each gene and patient
- Minimum group size: at least 10 ALD/DLD tumors and at least 10 other diploid reference tumors
- Delta methylation cutoff used for summaries: absolute delta methylation >= 0.20

Promoter methylation matrix cache filenames include the cohort, promoter regions, platform, and sample type. These matrices are reused across ALD/DLD tests and group-size or delta-cutoff settings. Result filenames additionally include the group-size thresholds, delta methylation cutoff, and core workflow setting tags to distinguish analysis outputs.

## Run order

Run scripts from the Workflow project root; script paths below are relative to that root.

1. Run the core CNmimicry workflow first so ALD/DLD calls and eligible amplified/deleted gene sets exist.
2. Run `analysis/methylation/scripts/run_methylation_ALD_screen.R`.
3. Run `analysis/methylation/scripts/run_methylation_DLD_screen.R`.
4. Run `analysis/methylation/scripts/summarize_methylation_ALD_results.R`.
5. Run `analysis/methylation/scripts/summarize_methylation_DLD_results.R`.

The screening scripts write:

- Per-cohort ALD or DLD promoter methylation Wilcoxon results.
- Per-cohort input-overlap summaries.
- Setting-aware promoter methylation matrix caches.

The summary scripts write cross-cancer summaries under `analysis/methylation/outputs/fdr_0.05`, with the ALD TPM tag in filenames when the requirement is enabled.

## Input-overlap summaries

For each cohort, the screening scripts save input-summary CSV files. These record the number of methylation patients, CNV patients, ALD/DLD-call patients, patient overlaps, eligible genes with promoter methylation data, and genes tested after the minimum sample filters.

These files are intended for reproducibility checks and for explaining why sample or gene counts differ between cohorts.
