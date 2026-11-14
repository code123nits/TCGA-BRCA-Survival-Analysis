# Evaluating the Added Predictive Value of RNA-Sequencing Gene Expression for Overall Survival in Breast Cancer

Nitsuh Wubishet Boru | Exploratory TCGA-BRCA research manuscript | 13 November 2026

**Status:** Research draft; not peer reviewed; model fitting not independently reproduced in this package.

## Abstract

Background: It remains unclear whether tumor RNA expression improves survival prediction beyond routinely available clinical variables in this cohort. We compared elastic-net Cox models using age and recorded pathological stage with models that also included 500 training-selected RNA-expression features. The matched TCGA-BRCA cohort comprised 1,072 patients (149 recorded deaths); 857 were assigned to training and 215 to testing (30 deaths). Archived test-set Harrell C-indices were 0.6712 for the clinical model, 0.6175 for the combined model at alpha 0.03, and 0.5854 at alpha 0.05. Paired combined-minus-clinical differences were -0.0537 (95% bootstrap interval -0.2126 to 0.0972) and -0.0858 (-0.2379 to 0.0481). An earlier archived three-fold training-CV analysis produced similar test C-indices for clinical and combined models (0.5832 and 0.5854). The original fixed-penalty selection is not documented, the test data were inspected repeatedly, and no fresh model refit or external validation has been completed for this package. These exploratory comparisons do not establish an incremental predictive benefit from the selected expression features.

## 1. Introduction

Tumor expression measurements can reveal biological differences that are not fully described by age or stage. Whether those measurements improve out-of-sample survival prediction is a separate question. Adding hundreds of correlated features to a cohort with relatively few observed deaths can make regularized models sensitive to feature selection and penalty choice. This project examines the incremental predictive performance of 500 high-variance RNA features beyond a simple clinical comparator in TCGA-BRCA. It is an exploratory analysis and reproducibility audit, not a clinical prediction tool.

## 2. Data and methods

The archived analysis matched TCGA-BRCA clinical and tumor RNA-sequencing records for 1,072 patients, including 149 recorded deaths. An 80:20 stratified split (random state 42) assigned 857 patients (119 deaths) to training and 215 (30 deaths) to testing. Overall survival used recorded survival time and death status. The expression matrix contained 19,962 protein-coding features transformed as log2(TPM + 1). Raw-data construction and file-selection audits are described in documentation/DATA_ACCESS.md, but the raw GDC build script and source manifest are not included in this package.

Age was median-imputed and standardized using training observations. Pathological stage was encoded using recorded categories, with Stage I as reference, yielding 13 clinical columns in total. The 500 highest-variance genes were selected using training data and standardized using training statistics. The combined model had 513 columns. CoxnetSurvivalAnalysis was specified with elastic-net mixing parameter l1_ratio=0.5. Archived fixed-penalty specifications were clinical alpha 0.001, clinical plus RNA alpha 0.03, and clinical plus RNA alpha 0.05. The original selection rationale for these values is unknown.

Harrell's concordance index summarized ranking discrimination. The repository includes saved test predictions and a patient-level paired bootstrap (2,000 resamples, seed 42) to obtain percentile intervals for C-indices and paired differences. These intervals condition on the archived predictions and do not incorporate model-fitting, feature-selection, or penalty-selection uncertainty. The originally held-out test split has been inspected more than once and cannot be described as untouched confirmatory validation.

## 3. Archived results

The clinical model's test C-index was 0.6712 (95% bootstrap interval 0.5243 to 0.8040). The combined model at alpha 0.03 had a C-index of 0.6175 (0.4848 to 0.7420), and the combined model at alpha 0.05 had 0.5854 (0.4380 to 0.7363). Relative to the clinical model, paired differences were -0.0537 (-0.2126 to 0.0972) and -0.0858 (-0.2379 to 0.0481). Both intervals crossed zero. These numbers were checked against the saved predictions, not reproduced by refitting the original models.

The archived clinical fit retained 13 nonzero coefficients. The combined alpha 0.03 and alpha 0.05 fits retained 58 and 21 coefficients, respectively. These are penalty-dependent feature-retention counts, not evidence that the selected genes are mechanistically relevant.

In an earlier three-fold training-cohort cross-validation analysis, the clinical alpha was 0.01 and the combined alpha was 0.05. The corresponding archived test C-indices were 0.5832 and 0.5854. This comparison differs from the fixed-penalty comparison and should not be merged with it. The clinical alpha was at the edge of the tested grid, which omitted alpha 0.001. A newer repeated-CV script is provided but has not been run on the real inputs.

A feasibility check of SCAN-B (GSE96058) found differences in clinical-variable availability and expression representation. No SCAN-B predictions or external-validation C-indices were computed.

## 4. Discussion

Under the archived fixed-penalty specifications, the clinical model had a higher test C-index than either model with 500 RNA features. This difference does not establish that RNA expression lacks prognostic value. The uncertainty intervals are broad, only 30 test patients had a recorded death, and the relative ordering changed when the penalties were selected in an earlier exploratory cross-validation analysis.

Variance selection identifies highly variable genes, not necessarily genes associated with survival. Penalization can also change which clinical or molecular coefficients remain nonzero. The absence of stage coefficients in the combined archived fits therefore does not show that stage is biologically unimportant.

The main lesson of this analysis is methodological: conclusions about incremental predictive value depend on a transparent comparison, reliable preprocessing, prespecified tuning, and independent evaluation. The saved predictions make part of the analysis auditable, but they do not establish that the full fitting pipeline can be reproduced.

## 5. Limitations

The analysis has several limitations. The test set contains 30 deaths, so estimates are imprecise. Fixed-penalty choices lack documented original selection rationale, and the test set was inspected repeatedly. Fine-grained stage categories may be sparse. The timing and selection of some tumor specimens relative to the survival origin are incompletely documented. The bootstrap does not account for training uncertainty. The full raw-to-processed pipeline and a fresh real-data model refit are not verified in this package. Proportional-hazards assumptions, calibration, time-dependent performance, clinical utility, and external validation were not assessed for the archived primary comparison.

## 6. Conclusion

This exploratory TCGA-BRCA comparison did not demonstrate improved archived test-set discrimination after adding 500 variance-selected RNA features to age and stage. The result is sensitive to penalty selection and limited by a small number of test events and incomplete reproducibility. Further work would require verified input provenance, training-only model selection, an independently evaluated cohort, and appropriate survival-model diagnostics.

## Reproducibility and data availability

The repository provides archived predictions with sequential sample labels in place of TCGA identifiers, metrics, figures, model code, bootstrap code, and hashes for three processed inputs. The processed inputs are not redistributed, and the raw GDC construction script and manifest are absent. Archived prediction-based calculations can be audited; fresh fitting and end-to-end reproduction remain outstanding. The data originate from NCI GDC TCGA-BRCA; GEO GSE96058 was considered for external-validation feasibility only. This manuscript is a research draft, not a peer-reviewed publication or journal-ready validation study.

## References

[1] The Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. Nature. 2012;490:61–70. doi:10.1038/nature11412.
[2] Cox DR. Regression models and life-tables. Journal of the Royal Statistical Society, Series B. 1972;34:187–220.
[3] Tibshirani R. The lasso method for variable selection in the Cox model. Statistics in Medicine. 1997;16:385–395.
[4] Zou H, Hastie T. Regularization and variable selection via the elastic net. Journal of the Royal Statistical Society, Series B. 2005;67:301–320.
[5] Collins GS et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. BMJ. 2024;385:e078378.
Bibliographic details should be checked against publisher records before journal submission.

