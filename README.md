---
output:
  html_document: default
  pdf_document: default
---

<!-- 
README for running multisample DPClust pipeline  
Created by: Miaomiao Gao  
Contact: miaomiao.gao@postgrad.manchester.ac.uk  
-->

# Pipeline for Multisample DPClust

This pipeline is designed to cluster somatic SNVs from multiple tumor samples of a single individual and reconstruct the phylogenetic tree.

---

## Prerequisites

Before running the pipeline, you must prepare the following:

### 1. Run Battenberg and SNV Calling

- Use **Battenberg** to call CNVs and a mutation caller (e.g., Mutect2) to generate an SNV VCF file for each tumor sample.
- Prepare the paths to your BAM and SNV files: `Bam_file_path` and `snv_path`.

The expected structure is that the VCF files are stored under `${snv_path}/${sample_name}/0001.vcf`.  
If your folder structure differs, you can modify:
- `snv_path` in `00_generate_loci_AF.R`
- `Bam_file_path` in `01_allele_count.txt`

---

### 2. Create the `run_info.txt` File

This file stores sample metadata for one case. You can define the case name arbitrarily.  
The `useBBfolder` column is used to locate subclone and purity information from the Battenberg output.

Example:

| case              | Tumour                      | useBBfolder                      | Gender | Normal                    |
|-------------------|-----------------------------|----------------------------------|--------|----------------------------|
| DW111230_original | DW5_11_C009453T1FTa_S1      | DW5_11_C009453T1FTa_S1           | Male   | DW5_35_C009453T2Wa_S25     |
| DW111230_original | DW5_12_C009453T1FTb_S2      | DW5_12_C009453T1FTb_S2           | Male   | DW5_35_C009453T2Wa_S25     |
| DW111230_original | DW5_30_C009453T3FTa_S20     | DW5_30_C009453T3FTa_S20_reft_WGD | Male   | DW5_35_C009453T2Wa_S25     |

Make sure you define the following:
- `BB_file_path`: path to Battenberg outputs
- `run_info_path`: path to the `run_info.txt` file  
Note: The `subclone.txt` file is expected at `${BB_file_path}/${sample_name}/subclone.txt`.

---

### 3. Configure `pipeline.sh`

Edit `pipeline.sh` and set the following parameters:
- `case`: e.g., `DW111230_original`
- `sample_num`: number of tumor samples (e.g., 3)
- `run_in_parallel`: number of 3-sample combinations from total samples.  
  For example:
  - 8 samples → 56 combinations
  - 3 samples → 1 combination

You can either pass these values as parameters or hardcode them in the `pipeline.sh` script.  
See the provided `pipeline.sh` for examples.

---

### 4. Run the Pipeline

Execute the pipeline using:

```bash
sh pipeline.sh 1
```
1 here refers to the patient number (you can define your own indexing).



## Pipeline Overview

The pipeline includes the following 7 main steps:

### Step 0: `00_generate_loci_AF.R`
Generates a unified list of loci from all tumor samples of a given case.

### Step 1: `01_allele_count.txt`
Counts the reference and alternative read counts at each locus for each sample.

> **Tip**: If you encounter a "Permission Denied" error, run:
> ```bash
> chmod 777 alleleCounter
> ```

### Step 2: `02_GenerateMutWT.py`
Generates `wt_count` and `mut_count` files required for downstream clustering.

### Step 3: `03_Get_DP_Input.txt`
Creates input files compatible with DPClust.

### Step 4: `04_runDP_GS.txt`
Performs Gibbs sampling using DPClust to infer cluster structures.

### Step 5: `05_parallel_assign_mutations.txt`
Assigns mutations to clusters using combinations of 3 samples at a time.

### Step 6: `06_combine.txt`
Combines all assignment results across combinations to generate a final clustering outcome.

> **Note**:  
> In the single-sample DPClust workflow, Steps 4–6 are executed together.  
> However, for multisample analysis (sample number > 3), Step 5 is run on all 3-sample combinations to avoid memory issues, and Step 6 aggregates the results.

---

## Notes

- Ensure that all file paths (e.g., `snv_path`, `Bam_file_path`, `BB_file_path`) are correctly specified in the relevant scripts.
- The final outputs include mutation cluster assignments, inferred tree topologies, and associated metrics.
- The pipeline assumes that DPClust and all required dependencies are properly installed in your environment.

---

## Contact

For questions or troubleshooting, please contact:

**Miaomiao Gao**  
📧 miaomiao.gao@postgrad.manchester.ac.uk
