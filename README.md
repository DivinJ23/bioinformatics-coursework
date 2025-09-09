# Bioinformatics Assignment – RNA-seq Reanalysis of Blanco-Melo et al. (2020)

## 📖 Overview
This repository contains my bioinformatics assignment based on the paper by **Blanco-Melo et al. (2020)**:  
[“Imbalanced Host Response to SARS-CoV-2 Drives Development of COVID-19”](https://doi.org/10.1016/j.cell.2020.04.026).

The aim of this assignment was to re-analyze a subset of the RNA-seq data presented in the paper, focusing on infections caused by RSV, IAV, and SARS-CoV-2. The assignment involves both coding in R and a 1500-word written report.

---

## 📂 Repository Contents
- **`bioinformatics_assignment.R`** – R script used to perform the RNA-seq data processing and visualizations.
- **`bioinformatics_assignment.docx`** – Final written report (1500 words).
- **`README.md`** – This file, providing context for the project.

---

## 🧪 Assignment Structure
The report is divided into the following sections:

### **Section A – Background on RNA-seq (600 words, 40 marks)**
- Explanation of RNA-seq methodology.
- Experimental design considerations (controls, replicates).
- Pre-sequencing steps.
- Data preparation before running DESeq2.

### **Section B – R Coding (150 words + code, 10 marks)**
- Well-formatted R code demonstrating:
  - PCA plot generation.
  - Sample distance matrix.
  - Differential expression (volcano plot).
- Short explanations accompanying each code block.

### **Section C – Analysis Results (600 words, 40 marks)**
- Presentation of data re-analysis, including:
  - PCA plot of A549 cell infections (RSV, IAV, SARS-CoV-2).
  - Sample distance matrix.
  - Volcano plot of RSV vs. SARS-CoV-2 response.
  - Heatmap of interferon gene expression.
- Narrative-driven interpretation of the figures.

### **Section D – Next Steps (150 words, 10 marks)**
- Suggested follow-up experiments.
- Discussion of the advantages and limitations of RNA-seq.

---

## ⚙️ Tools & Packages
This project was completed in **R** using the following packages:
- `DESeq2`
- `ggplot2`
- `pheatmap`
- `dplyr`
- `RColorBrewer`

---

## 📌 Notes
- Report follows university formatting requirements (12pt font, 1.5 line spacing, 2 cm margins).
- Student anonymity preserved – only student number included in the report.
- No screenshots were used; all code and plots were generated directly in R.

---

## 📜 Reference
Blanco-Melo D, Nilsson-Payant BE, Liu W-C, et al. (2020).  
*Imbalanced Host Response to SARS-CoV-2 Drives Development of COVID-19.*  
Cell, 181(5), 1036–1045.e9.  
https://doi.org/10.1016/j.cell.2020.04.026
