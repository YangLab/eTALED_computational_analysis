# eTALED_computational_analysis
This project provides relavant code for analyzing the efficient ABE of mitochondrial DNA

Author: Bao-Qing Gao

Email: gaobaoqing2019@sibs.ac.cn

---
## 1. Features  
- **Editing Efficiency Calculation** `scripts_for_base_substitution_analysis`: Evaluates engineered ABE editing efficiency at target sites with the deep sequencing data.
- **RNA Off-Target Effect Analysis** `scripts_for_RNA_editing_analysis`: Detects potential off-target editing events with the whole-transcriptome sequencing data.
- **Whole Nuclear Genome Off-Target Effect Analysis** `scripts_for_genome_editing_analysis`: Detects potential off-target editing events with the whole nuclear genome sequencing data.
- **Whole Mitochondrial Genome Off-Target Effect Analysis** `scripts_for_mitochondrial_genome_editing_analysis`: Detects potential off-target editing events in whole mitochondrial genome with the deep sequencing data.

---

## 2. Prerequisites  
Before running this analysis, install the following tools:  
- **CFBI**: [https://github.com/YangLab/CFBI](https://github.com/YangLab/CFBI)  
- **RADAR**: [https://github.com/YangLab/RADAR](https://github.com/YangLab/RADAR)
- **BEIDOU**: [https://github.com/YangLab/BEIDOU](https://github.com/YangLab/BEIDOU)
---

## 3. Usage  

### **PART 1: Editing Efficiency Calculation**  
- **Directory**: `scripts_for_base_substitution_analysis`  

  All original deep sequencing data used in this study can be accessed at the NCBI BioProject under the accession code **PRJNA1223321**. The provided code here serves as an example for reference.

- #### **Constructing Index Files**  
  Before running the analysis, you need to construct index files. Refer to the `index/` directory, which contains:  
  
  - `20240220_deepseq_sample_name_mapping.csv`:
    
    | ID | Sample 1 | Target site 1 | Sample 2 | Target site 2 | ... |
    | -------------- | -------- | ------ | -------- | ------ | ------ |

  - `ref_20240220.txt`:
    
    | Information                                |
    |--------------------------------------------|
    | Target site 1                             |
    | Sequence of target site 1                 |
    | Editing window for target site 1: Editing window sequence 1 |
    | None                                           |
    | None                                             |
    | Target site 2                             |
    | Sequence of target site 2                 |
    | Editing window for target site 2: Editing window sequence 2 |
    | None                                            |
    | None                                            |
    | ...                                        |

- #### **Calculating Editing Efficiency**  
  Run  
  ```bash
  sh CFBI_main.sh

### **PART 2: RNA Off-Target Effect Analysis**  
- **Directory**: `scripts_for_RNA_editing_analysis`

  All original whole-transcriptome sequencing data used in this study can be accessed at the NCBI BioProject under the accession code **PRJNA1221787**. The provided code here serves as an example for reference.

- #### **Quality Control & Trim Reads**  
  Run  
  ```bash
  sh work_QC_Trim.sh

- #### **Whole-Transcriptome Off-Target Effect Analysis**  
  Run  
  ```bash
  sh work_RADAR_parallel.sh

### **PART 3: Whole Nuclear Genome Off-Target Effect Analysis**  
- **Directory**: `scripts_for_nuclear_genome_editing_analysis`

  All original whole-genome sequencing data used in this study can be accessed at the NCBI BioProject under the accession code **PRJNA1221509**. The provided code here serves as an example for reference.

- #### **Quality Control & Trim Reads**  
  Run  
  ```bash
  sh work_QC_Trim.sh

- #### **Whole Nuclear Genome Off-Target Effect Analysis**  
  Run  
  ```bash
  sh work_BEIDOU.sh

  
### **PART 4: Whole Mitochondrial Genome Off-Target Effect Analysis**  
- **Directory**: `scripts_for_mitochondrial_genome_editing_analysis`  

  All original deep sequencing data used in this study can be accessed at the NCBI BioProject under the accession code **PRJNA1223321**. The provided code here serves as an example for reference.

- #### **Preparing Input Files**  
  Before running the analysis, you need to prepare input files. Refer to the `input_files/` directory, which contains:  
  
  - `CFBI_output_files_test`:
    The CFBI outputs from PART 1 were used to calculate mitochondrial genome-wide average off-target editing frequency. CFBI output files for several samples are provided in this directory for testing.
    
  - `sample_edit_site.txt`:
    
    | Sample ID | Editing positon | Editing window |
    | -------------- | -------- | ------ |

- #### **Calculating Average Off-target Editing Frequency**  
  Run  
  ```Rscript
  Rscript mit_offtarget.r

---

## 4. Citation
Fan Y#, Xu W#, Gao BQ#, Qin H, Wu X, Wei J, Ni Q, Zhou L, Xiang J, Wu J, Yang B, Yang L and Chen J*. Leveraging base excision repair for efficient adenine base editing of mitochondrial DNA. ***Nat Biotechnol***, 2025 (Accepted in principle).
