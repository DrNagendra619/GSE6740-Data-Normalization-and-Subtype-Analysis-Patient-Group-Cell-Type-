# GSE6740-Data-Normalization-and-Subtype-Analysis-Patient-Group-Cell-Type-
GSE6740 Data Normalization and Subtype Analysis (Patient Group &amp; Cell Type)
# 📊 Microarray Preprocessing Pipeline: Multivariate Clustering of HIV/AIDS Patients (GSE6740)

This R script automates the essential preprocessing and Quality Control (QC) steps for analyzing the **GSE6740** microarray dataset. This complex study investigates gene expression in T-cells from individuals with different stages of **HIV/AIDS**.

The pipeline is distinguished by its ability to extract and simultaneously visualize **two critical biological factors (Patient Group and Cell Type)** in the Principal Component Analysis (PCA) plot, ensuring the data is properly prepared for multivariate statistical modeling.

## 🚀 Key Features

* **Automated Data Retrieval:** Fetches expression data and metadata directly from the **GEO database (GSE6740)** using `GEOquery`.
* **Quantile Normalization:** Applies **quantile normalization** (`limma::normalizeBetweenArrays`) to correct technical variation and ensure comparable gene distributions across all samples.
* **Dual Factor Extraction:** Extracts and creates factors for **two independent variables** from the metadata:
    1.  **Patient Group:** (Acute, Chronic, Non-Progressor, Uninfected)
    2.  **Cell Type:** ($\text{CD4}^{+}$ and $\text{CD8}^{+}$ T-cells)
* **QC Visualization:** Generates **boxplots** before and after normalization to visually confirm data distribution consistency.
* **Low-Expression Filtering:** Filters out genes with mean expression in the lowest quartile (25th percentile) to reduce noise.
* **Multivariate Clustering:** Generates a **PCA plot** that simultaneously uses **color for Patient Group** and **shape for Cell Type**, providing a powerful visual QC check.
* **Data Persistence:** Saves the fully processed data, metadata, and both factor variables to an `.RData` file for subsequent analysis steps.

---

## 🔬 Analysis Overview

| Component | Method / Test | Purpose |
| :--- | :--- | :--- |
| **Dataset** | GSE6740 | Gene expression from $\text{CD4}^{+}$ and $\text{CD8}^{+}$ T-cells of HIV/AIDS patients. |
| **Normalization** | Quantile Normalization | Standardizes gene expression distributions across all arrays. |
| **Filtering** | Interquartile Range (IQR) Filtering | Removes non-informative, low-expressed genes. |
| **Clustering** | PCA | Assesses global sample similarity and checks if the two experimental factors ($\text{CD4}^{+}/\text{CD8}^{+}$ and disease stage) drive the main sources of variance. |

---

## 🛠️ Prerequisites and Setup

### 📦 Packages

The script loads the following essential packages:
* `GEOquery` (For data retrieval)
* `limma` (For normalization)
* `dplyr` (For metadata manipulation)
* `ggplot2` (For PCA visualization)

### ⚙️ Execution

1.  **Download** the `GSE6740 Data Normalization and Subtype Analysis (Patient Group & Cell Type).R` file.
2.  **Execute** the script in your R environment:
    ```R
    source("GSE6740 Data Normalization and Subtype Analysis (Patient Group & Cell Type).R")
    ```
    *Note: All output files are saved to the current working directory where the script is executed.*

---

## 📁 Output Files (3 Plots + 1 Data File)

| Filename | Type | Description |
| :--- | :--- | :--- |
| `GSE6740_processed_data_step1.RData` | R Binary Data | Contains the final, filtered, and normalized `expression_data`, `metadata`, `patient_group`, and `cell_type` factors. |
| `GSE6740_boxplot_before_normalization.png` | QC | Boxplot illustrating raw data distributions. |
| `GSE6740_boxplot_after_normalization.png` | QC | Boxplot confirming uniform data distributions across all samples after normalization. |
| `GSE6740_pca_plot.png` | Clustering | **Principal Component Analysis (PCA)** plot, showing sample clustering simultaneously colored by **Patient Group** and shaped by **Cell Type**. |
