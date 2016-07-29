# Analysis for "Chromatin environment informs CRISPR-Cas9 system activity"

---

Data and code to reproduce the analysis performed in **"Chromatin environment informs CRISPR-Cas9 system activity"** are available in this repository.

---
### Data

Feature tables for the four different datasetes are collected in **Data/**

The following information is collected for each Cas9 target site:

1. Predicted Activity score from sgRNAscorer
2. Predicted Activity score from WU-CRISPR
3. Predicted Activity score from Combined sgRNAscorer and WU-CRISPR scores
4. Experimentally determined score
5. Target site activity class (0=Low, 1=High)
6. Intersect information of the target site with 32 chromatin marks from **Ernst & Kellis, Nature Biotech (2015)** 

---
### Code
The analysis script for construction and testing of the model is outlined in **Results/Workflow.R**

---
### Results
Output of **Results/Workflow.R** are stored in **Results/**