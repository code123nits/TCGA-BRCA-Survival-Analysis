# Raw GDC to processed-matrix build: status

**Not verified in this repository.** The two processed inputs (`tcga_brca_survival_dataset.csv` and
`tcga_brca_matched_expression.pkl`) are described in `DATA_ACCESS.md`, but the script that builds them from raw GDC
downloads, the GDC file manifest (file UUIDs and MD5 checksums) and the audit outputs are not included here.
An earlier project audit reported MD5 matches for 1,095 source files and exact reconstruction of the matrix;
those audit outputs have not been re-verified in this package.

To close this gap: add the build script and manifest, rerun the reconstruction, and record the result here.
