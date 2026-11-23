# Fixed 13-gene expression scoring for lower-grade glioma

This tool calculates the fixed 13-gene expression score used in the TCGA-LGG and CGGA analyses. It uses base R only; no additional packages are required. Tested with R 4.5.3.

## Files

| File | Type and contents |
|---|---|
| `score_signature.R` | Scoring functions and command-line entry point. |
| `score_model.rda` | The `score_model` parameter list: 13 original-precision coefficients, gene mappings, expression-processing rules and original-cohort score references. This is a fixed parameter object, not a fitted `glmnet` or clinical Cox object. |
| `example_expression.rda` | The `example_expression` numeric matrix: 13 genes by one CGGA693 sample, labelled `Example_1`. |
| `README.md` | Input requirements, commands, formulas and worked example. |

Keep the four files in one directory. There is no R package to install.

## Run the worked example

From that directory:

```sh
Rscript score_signature.R example_expression.rda CGGA693 CGGA693
```

The second argument specifies expression identifiers/platform; the optional third argument explicitly requests the matching original cohort's score reference. The command prints:

```text
    sample          raw_score             score_z risk_group
 Example_1 1.2228083012082467 0.70339263643909866       High
```

For the raw score alone, omit the reference argument:

```sh
Rscript score_signature.R example_expression.rda CGGA693
```

No output files are written automatically.

The same functions can be called from R:

```r
source("score_signature.R")
model <- load_score_model("score_model.rda")
x <- read_expression_rda("example_expression.rda")
score_signature(x, "CGGA693", model)  # raw score only
score_signature(x, "CGGA693", model, reference = "CGGA693")
model$features  # all 13 coefficients and their identifier mappings
```

## Expression input

Supply a numeric matrix with **genes in rows and samples in columns**. Gene identifiers must be row names; sample column names must be present and unique. The matrix may contain just the 13 model genes or additional genes. Selected model-gene values must be finite, non-negative and nonmissing. Measured zeros are valid; missing model genes are not replaced with zero.

For the command-line interface, save exactly one such matrix in an `.rda` file using R's `save()` function. The matrix object name can be chosen freely. For the R interface, pass the matrix directly.

**Input must not already be log-transformed.** Both platforms undergo `log2(expression + 1)` once. The tool cannot infer from numeric values whether a matrix has already been logged: supplying logged values as ordinary input would transform them again. The R argument `input_scale` accepts only `"untransformed"`; declaring another scale produces an error, not a reverse transformation.

- **TCGA:** use the study's FPKM-UQ expression scale (upper-quartile-normalised FPKM), with Ensembl gene IDs. This is not a converter from read counts, TPM or another expression scale. Ensembl version suffixes are removed. For duplicated target IDs, retain the row with the highest mean **untransformed** expression across all supplied sample columns; ties retain the first input row. A warning reports the selected row. Exact reproduction of duplicate resolution requires the original preprocessing sample inventory (534 input columns before the final 511-case analysis selection), or rows already resolved using that rule.
- **CGGA693 / CGGA325:** use the original mRNAseq_693 / mRNAseq_325 **RSEM-genes** expression-matrix values, releases labelled **20200506**, with gene symbols as row names. Do not substitute a different RSEM output scale. The EPOP coefficient matches **C17orf96** in both CGGA matrices. Duplicate model symbols produce an error and must be resolved explicitly before scoring.

The identifier mappings are stored in `score_model$features`; no external annotation file is needed. There is no per-gene z-scoring, cross-platform quantile normalisation or batch correction in this tool.

## Fixed score and calculation details

For each sample:

```text
raw_score = sum_j coefficient_j * log2(expression_j + 1)
```

The coefficients are the 13 nonzero TCGA Cox-LASSO coefficients at `lambda.min`, preserved as the original R double values, rather than rounded manuscript coefficients. No model is fitted during scoring. Gene matching and duplicate resolution precede log transformation.

To inspect the complete 13-gene worked calculation in R:

```r
calculation <- score_signature(x, "CGGA693", model,
                               reference = "CGGA693", details = TRUE)
calculation$contributions
calculation$scores
```

`contributions` has one row per gene and sample, containing `sample`, `gene` (Ensembl ID), `symbol`, `input_gene_id`, `expression`, `log2_expression_plus1`, `coefficient` and `contribution`. `mapping_audit` additionally records selected input row numbers and match counts. These are returned in R, not stored as separate example-results files.

For EPOP in the supplied example:

| Model symbol | Input identifier | Expression | log2(expression+1) | Coefficient | Contribution |
|---|---|---:|---:|---:|---:|
| EPOP | C17orf96 | 14.44 | 3.9486008474933558 | -0.061143919374343082 | -0.2414329318605965 |

Summing all 13 contributions gives the raw score **1.2228083012082467**.

## Optional score standardisation and grouping

The references in `score_model$reference_parameters` describe the **full original analysis cohorts, before clinical complete-case restriction**:

| Reference | n | Mean raw score | Sample SD | Median raw score |
|---|---:|---:|---:|---:|
| TCGA | 511 | 0.1171609258818581 | 0.75608604126647294 | -0.1404758319828866 |
| CGGA693 | 271 | 0.5481727749293579 | 0.95911656069390805 | 0.41015503265843312 |
| CGGA325 | 137 | 1.424597888950194 | 0.98173376807527524 | 1.1740596634664893 |

Only when the matching reference is explicitly requested:

```text
score_z = (raw_score - mean_raw) / sd_raw
risk_group = High if raw_score >= median_raw, otherwise Low
```

SD uses the sample denominator n-1. CGGA693 uses its own reference, and CGGA325 its own; neither uses the TCGA mean/SD. A reference that does not match `platform` is rejected. The program does not estimate a reference from the supplied input or call `scale()` on a single sample.

For `Example_1`, the CGGA693 reference gives **score_z = 0.70339263643909866** and **High**. The supplied references reproduce the original cohorts; they are not automatically estimated reference values or universal clinical thresholds for a new dataset.

## Output and scope

By default, `score_signature()` returns a data frame with `sample` and `raw_score`. Explicit reference selection adds `score_z` and `risk_group`. With `details = TRUE`, it returns a list containing `scores`, `contributions` and `mapping_audit`.

The tool starts from quantified expression and ends at the fixed expression score, with optional reference standardisation/grouping. It does not perform gene selection, LASSO training, Cox fitting or performance evaluation. It does not calculate the separate age/grade/subtype-plus-score predictor or absolute survival probabilities.

## Data source and verification

The example is one actual CGGA mRNAseq_693 profile restricted to the 13 model genes. Its original sample identifier is replaced by `Example_1`; no clinical or survival fields are included. CGGA: https://www.cgga.org.cn/download.jsp. TCGA-LGG FPKM-UQ source: https://portal.gdc.cancer.gov/projects/TCGA-LGG.

The four-file tool was checked against the original TCGA511, CGGA693271 and CGGA325137 saved scores (919 samples in total). In each cohort, the maximum absolute raw-score difference was 4.885e-15, and all original cohort-median group labels matched. The maximum absolute score-z difference was 6.6613381477509392e-15. The coefficients, mappings and reference values were unchanged; no model was refitted.
