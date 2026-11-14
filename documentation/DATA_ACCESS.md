# Data access and expected input files

The inputs are **not** redistributed in this repository (the processed expression matrix is ~168 MB, and
TCGA/GDC and GEO data are subject to their own access terms). `frozen_refit.py` and `repeated_cv.py` need
three files in one directory (`--data`):

| File | Contents |
|---|---|
| `tcga_brca_survival_dataset.csv` | One row per patient. First (unnamed) column = TCGA patient ID. Columns used: `event` (0/1 death), `survival_days`, `age_at_diagnosis`, `pathologic_stage` (stage label as recorded; missing is treated as `Unknown`). After matching to expression: 1,072 patients, 149 events. |
| `tcga_brca_matched_expression.pkl` | pandas DataFrame, rows = the same patient IDs, columns = 19,962 protein-coding genes, values = log2(TPM-unstranded + 1). Pickle format: only load files you created yourself. |
| `test_patient_ids.csv` | The 215 held-out patient IDs (header row, then one ID per line). Can be regenerated with `train_test_split(..., test_size=0.2, stratify=event, random_state=42)`. The scripts assert that it matches. |

## Where the data come from
* Clinical and survival data and RNA-seq: NCI Genomic Data Commons, project TCGA-BRCA (open-access data).
* Expression: GDC STAR-Counts files, using the TPM-unstranded column, transformed as log2(TPM + 1).
* External-feasibility cohort (not used for model scoring): GEO GSE96058 (SCAN-B).

## Processing steps recorded in the project notes (verify against your audit files)
* 1,095 source RNA files were selected; MD5 checks matched for all of them and the processed matrix was
  reconstructed exactly (21,858,390 values).
* Some patients had several candidate RNA files (an earlier draft recorded 11 patients with 27 candidate
  files, and five patients with multiple distinct 01A files chosen alphabetically, all in the training set).
  Alphabetical choice is reproducible but is not a biological quality criterion.

## Still missing from this package (TODO for the author)
The script that built the two data files from raw GDC downloads, the GDC file manifest (file UUIDs and MD5
checksums), and the audit outputs are not included. Until they are added, other people can rerun the model
code only if they obtain equivalent processed files, and the raw-to-matrix step cannot be reproduced from
this repository. Do not describe the pipeline as reproducible end to end until these are added.
