This guide provides a comprehensive overview of the **DNA Structural Analyzer**, a professional-grade Streamlit application designed for thermodynamic modeling and engineering of genetic sequences.

---

# DNA Structural Analyzer

A high-fidelity tool for molecular biologists and bioengineers to visualize DNA topology, analyze thermodynamic stability, and perform targeted mutagenesis. The application leverages nearest-neighbor approximations to identify structural weak points and restriction enzyme recognition sites.

##  Features & Functionality

### 1. 3D Structural Analytics

* **Dynamic Topology**: Generates a mathematically accurate double-helix model based on your sequence length, featuring a fixed **3-turn rotation** for visual consistency.
* **Stability Heatmapping**: Rungs are color-coded based on local melting temperatures ($T_m$):
* **Stable (Green)**: $T_m > 20°C$.
* **Partial (Blue)**: $16–20°C$.
* **Unstable (Red)**: $T_m < 16°C$.


* **Interactive Engine**: Includes an **Auto-Rotate** mode and a **Focus System** that snaps the camera to specific indices or enzyme sites.

### 2. Genetic Engineering Suite

* **Mutation Registry**: Automatically identifies "Unstable" windows and suggests specific $A/T \to G/C$ substitutions to increase local structural integrity.
* **Enzyme Mapping**: Detects recognition sites for common endonucleases (EcoRI, BamHI, HindIII, NotI, SmaI, XhoI) and highlights them in gold on the 3D model.
* **Termini Analysis**: Calculates GC-clamp strength at both $5'$ and $3'$ ends and generates the reverse-complement sequence instantly.

### 3. Comparison & Version Control

* **Sequence Comparison**: Visually diffs your current sequence against historical snapshots. It uses a custom engine to highlight deletions (red strike-through) and additions (green bold).
* **History Archive**: Automatically logs manual edits, randomizations, and restorations with timestamps and custom labels.

### 4. Professional PDF Reporting

Generates a comprehensive one-page analytical report including:

* Thermodynamic assessment and grade (Stable/Unstable).
* Nucleotide composition and GC-content interpretation.
* A mapped table of restriction sites and suggested mutations.

---

## Installation & Running the Code

### Prerequisites

Ensure you have Python 3.8+ installed. The application requires the following libraries:

* `streamlit`
* `numpy`
* `plotly`
* `reportlab` (for PDF generation)

### Setup

1. **Install Dependencies**:
```bash
pip install streamlit numpy plotly reportlab

```


2. **Launch the Application**:
Navigate to the directory containing `dna_app_fixed.py` and run:
```bash
streamlit run dna_app_fixed.py

```



---

# How to Use

1. **Input Sequence**: Enter your DNA sequence (A, T, G, C) in the sidebar text area or click **"Generate random"** for a test 60bp strand.
2. **Analyze Stability**: Use the **"Analysis window"** slider to adjust the sensitivity of the $T_m$ calculation (default is 6bp).
3. **Explore the Helix**: In the **3D Analytics** tab, use the "Rotate" button to view the topology. Click **"Focus"** next to any enzyme in the Engineering tab to inspect that specific site.
4. **Optimize**: Review the **Mutation Registry**. Click **"Apply"** on suggested changes to automatically update your sequence with more stable bases.
5. **Compare & Export**: Save a snapshot in the **Comparison** tab before making major changes. Once finished, use the **"Export PDF report"** button to save your findings.

---

### Technical Note

Melting temperatures are calculated using the nearest-neighbor approximation ($4°C$ per $G/C$ and $2°C$ per $A/T$) within the sliding analysis window defined in the sidebar.
