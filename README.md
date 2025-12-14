# Automated Computational Workflow for $S_N2$ Reactions
### Computational Chemistry Project, Fall 2025
**Author:** Sam Sheng (Junior, Chemistry Major)

![Status](https://img.shields.io/badge/Status-Completed-success) ![Python](https://img.shields.io/badge/Python-3.9%2B-blue) ![Gaussian](https://img.shields.io/badge/Gaussian-16-orange) ![License](https://img.shields.io/badge/License-MIT-green)

## 📖 Project Overview
This project implements a **fully automated Python-Gaussian workflow** (AI-assisted) to study the nucleophilic substitution ($S_N2$) reaction between various alkyl amines and methyl iodide ($\mathrm{CH_3I}$).

Beyond traditional Quantum Mechanics (QM) calculations, this project integrates **Machine Learning (SISSO-inspired multivariate regression)** to bridge the gap between calculated activation energies ($E_a$) and molecular descriptors, revealing the subtle interplay between electronic effects and steric hindrance.

**Reaction Scheme:**
$$\mathrm{NR_3} + \mathrm{CH_3I} \rightarrow [\mathrm{R_3N \cdots CH_3 \cdots I}]^\ddagger \rightarrow \mathrm{MeNR_3^+} + \mathrm{I^-}$$

---

## 🚀 Key Features
* **End-to-End Automation**: Python scripts drive the entire Gaussian 16 pipeline (Input Gen $\to$ Job Submission $\to$ Log Parsing $\to$ Error Handling).
* **Robust Error Recovery**: Includes specialized "Rescue Scripts" for difficult cases like $\mathrm{Me_2NH}$ (Smart Scan) and $\mathrm{Me_3N}$ (Rigid Steric Correction).
* **AI-Driven Insight**: Inspired by  **SISSO** (Sure Independence Screening and Sparsifying Operator), utilizes multivariate regression to derive a unified physical law governing reactivity ($R^2 > 0.95$).

---

## 🛠️ Workflow

The project follows a rigorous 7-step computational pipeline:

1.  **Complex Assembly**: Automated generation of reactant complexes using `RDKit`.
2.  **Geometry Optimization**: Pre-optimization of the reactant complex structure.
3.  **Relaxed Scan (PES)**: Driving the C-N bond formation coordinate to locate the approximate energy barrier.
4.  **Transition State (TS) Optimization**:
    * *Zoom-in Strategy*: Refining the scan near the energy maximum.
    * *Pinpoint*: Using `Opt(TS, CalcAll, NoEig)` to locate the exact saddle point.
5.  **Frequency Analysis**: Verifying the TS by confirming a single imaginary frequency (vibrational mode corresponding to the reaction coordinate).
6.  **IRC Analysis**: Tracing the Intrinsic Reaction Coordinate path to ensure continuity between Reactants, TS, and Products.
7.  **🌟 Data Analysis and Feature Discovery (Where my whimsy comes in!):**
    * Extracting descriptors (HOMO, Charge, Entropy, Steric indices, etc.) from logs.
    * Applying SISSO to discover the "formula" of reactivity.

---

## 📂 File Structure

```text
Project_Root/
├── 00_Notebooks_for_Exploration/ # Demonstration of exploration
├── 01_Automated_Workflow/        # Core batch processing scripts (Scan -> TS -> IRC)
├── 02_Case_Studies/              # Troubleshooting scripts for Me2NH and Me3N outliers
├── 03_Data_Analysis/             # Feature engineering and multivariate regression
├── 04_Results_Data/              # Final CSV reports and PNG figures
├── 05_Gaussian_Inputs/           # Archived Gaussian input files
├── 06_Raw_Logs/                  # Archived Gaussian output files
├── 07_Plot_Drawing/              # Scripts for drawing plots.
└── utils/                        # Shared utility library (chem_utils.py)
```

---

## 📐 Methodology: Feature Engineering

Before applying the multivariate regression (SISSO-inspired), we extracted physical descriptors from both the **quantum mechanical outputs** (Gaussian 16) and the **molecular topology** (RDKit). This ensures the model captures electronic, steric, and geometric factors objectively.

The specific descriptors used in the unified model are defined as follows:

| Symbol | Physical Quantity | Physical Meaning | Calculation Source / Method |
| :--- | :--- | :--- | :--- |
| **$Q_N$** | **Nitrogen Partial Charge** | Represents **Electronic Effect**. A more negative charge indicates higher electron density and stronger intrinsic nucleophilicity. | **Gaussian .log File**<br>Extracted from *Mulliken Population Analysis* in the static reactant optimization step. |
| **$N_H$** | **Number of Amine Hydrogens** | Represents **Substitution Level**. <br>($NH_3$=3, $RNH_2$=2, $R_2NH$=1, $R_3N$=0). Acts as a discrete indicator for steric crowding around the N center. | **Molecular Graph**<br>Calculated via `RDKit` by counting H atoms attached to Nitrogen in the SMILES string. |
| **$\kappa_3$** | **Hall-Kier Kappa-3 Index** | Represents **Molecular Shape / Flexibility**. <br>Distinguishes between "globular/rigid" shapes (e.g., $Me_3N$, low $\kappa_3$) and "linear/flexible" chains (e.g., $NEt_3$, high $\kappa_3$). | **Topological Calculation**<br>Computed via `rdkit.Chem.GraphDescriptors.Kappa3()`. It encodes the degree of branching and flexibility based on the molecular skeleton. |

> **Note**: While $Q_N$ reflects the *intrinsic* reactivity calculated by QM, $\kappa_3$ and $N_H$ provide critical *geometric* corrections that purely electronic parameters usually miss.

---

## 📊 Key Findings

**THE FORMULA**：$$E_a = 297.9 + 933.82 \ Q_N + 17.31 \ N_H + 0.02 \ \kappa_3$$ $(R^2 = 0.9544)$

1.  **Electronic Dominance**: Generally, stronger electron-donating alkyl groups lower the activation barrier.
2.  **The "Steric Anomaly"**: $Me_3N$ was identified as a significant outlier. Unlike the flexible $Et_3N$, the rigid methyl groups in $Me_3N$ create a "rigid steric wall," significantly increasing the barrier.
3.  **Unified Model**: A 3-descriptor model was discovered, perfectly predicting gas-phase reactivity by combining **Electronic ($Q_N$)**, **Steric ($N_H$)**, and **Flexibility ($\kappa_3$)** factors.
   
> **Note on Dimensional Consistency**:
> Readers may observe that the input descriptors possess disparate units (e.g., $Q_N$ is in elementary charge, $\kappa_3$ is a dimensionless index). In this multivariate regression model, **dimensional consistency is implicitly maintained by the regression coefficients** ($c_1, c_2, c_3$). Each coefficient effectively carries the necessary physical dimensions to bridge the descriptors to the final energy unit ($\mathrm{kJ/mol}$).

---

## 🔧 Tech Stack

* **Python**: `pandas`, `numpy`, `matplotlib`, `scipy`, `sklearn`
* **Cheminformatics**: `RDKit`
* **QM Software**: Gaussian 16 (PM7 Method)
* **Visualization**: GaussView 6, PyMOL

---

## ✨ Acknowledgement

Special thanks to **Google's Gemini 3 Pro**.

This project was significantly accelerated by Gemini's exceptional capabilities in multi-turn conversational reasoning and multimodal processing. From organizing the workflow to debugging code and supporting my data "whimsy", Gemini acted as an invaluable research co-pilot.

>**Note:** 
>Using LLM to assist in this project is allowed and **highly encouraged** by the lecturer, adhering to all academic integrity and honor codes.

---

*Created by Sam Sheng. Fudan University, Class of 2027.*