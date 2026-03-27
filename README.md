This GitHub-specific README is tailored to the final version of your code, focusing on the sophisticated UI, thermodynamic calculations, and 3D rendering capabilities.

***

# DNA Structural Analyzer: 3D Topology & Stability Suite

A professional bioengineering tool built with Streamlit for the 3D visualization, thermodynamic modeling, and engineering of DNA sequences. This application utilizes nearest-neighbor approximations to identify structural instabilities and map restriction enzyme recognition sites.

##  Key Features

###  High-Fidelity 3D Rendering
* **Fixed Helical Geometry**: Generates a mathematically accurate double-helix with a locked 2.5-turn rotation for standardized 60bp analysis.
* **Real-time Camera Focus**: Includes a "Focus Engine" that recalculates 3D camera vectors to snap directly to specific nucleotide indices or enzyme cut sites.
* **Stability Isolation**: Toggle views to isolate rungs based on their stability grade (Stable, Partial, or Unstable).

###  Thermodynamic Analysis
* **Nearest-Neighbor Modeling**: Calculates local melting temperatures ($T_m$) using a sliding analysis window (default 6bp) to detect regions prone to "breathing" or denaturation.
* **Automated Mutagenesis**: Suggests $A/T \to G/C$ substitutions for unstable regions and allows for one-click application of these mutations.

###  Engineering Tools
* **Enzyme Registry**: Built-in detection for common endonucleases including EcoRI, BamHI, HindIII, NotI, SmaI, and XhoI.
* **Sequence Comparison**: A visual "Diff" engine that highlights changes between the current sequence and historical snapshots using color-coded HTML markers.
* **Automated PDF Reporting**: Generates a bench-ready summary report including GC-clamp strength, enzyme maps, and stability assessments.

---

##  Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/dna-structural-analyzer.git
   cd dna-structural-analyzer
   ```

2. **Install dependencies**:
   The app requires `streamlit`, `numpy`, `plotly`, and `reportlab`.
   ```bash
   pip install streamlit numpy plotly reportlab
   ```

3. **Run the application**:
   ```bash
   streamlit run dna_app_fixed.py
   ```

---

##  Usage Guide

1. **Input**: Enter a DNA sequence in the sidebar or use the **"Randomize"** button to generate a test sequence.
2. **Analyze**: Navigate to the **3D Structural Analytics** tab to inspect the helix. Select enzymes in the **Engineering** tab to highlight their binding sites in white.
3. **Compare**: Use the **Sequence Comparison** tab to view exactly how your mutations have changed the sequence relative to previous versions in your history.
4. **Export**: Once your design is optimized, click **"Export Narrative PDF Report"** to download your project documentation.

---

##  Built With
* **Streamlit** - The web framework for the UI.
* **Plotly** - Used for the 3D Scatter3d helical rendering.
* **NumPy** - Handles the trigonometric coordinate calculations for the DNA backbone.
* **ReportLab** - Powers the automated PDF generation.

## 📄 License
This project is licensed under the MIT License - see the LICENSE file for details.
