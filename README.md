# Metallic_Glasses

Welcome to the collaborative repository! This project aims to probe the structure-property relationship in metallic glasses.

---

## Project: Concurrent prediction of metallic glasses’ global energy and internal structural heterogeneity by interpretable machine learning

### 0 Sample Generation using LAMMPS
Inside the "sampleGeneration/" directory, we show an example for generation of sub-Tg (transition temperature) annealing samples.

### 1 Sample Generation
In folder npy/, we convert the LAMMPS output structures into the SOAP descriptors and store them as a matrix.

### 2 Main Part of Machine Learning
Generally, the "xgb/" folder contains the key components of our machine learning process, including dataset creation and training/testing.

#### 2.1 ULE
Inside the "LAEs/" directory, we use the defined threshold to select unique local environments (ULEs).

#### 2.2 Dataset preparation
In folder dataset/, the data is split and shuffled into three parts: training/validation/testing for further XGBoost processing.

#### 2.3 XGBoost
In Folder train/ and test/, we execute the XGBoost algorithm on the prepared dataset and subsequently analyze the obtained results.
