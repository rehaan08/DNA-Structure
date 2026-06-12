import streamlit as st
import random
import numpy as np
import plotly.graph_objects as go
import re
import io
from datetime import datetime
import subprocess, sys

for pkg in ["reportlab", "biopython"]:
    try:
        __import__("reportlab" if pkg == "reportlab" else "Bio")
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg,
                               "--quiet", "--disable-pip-version-check"])

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.MeltingTemp import Tm_Wallace, Tm_NN
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch, EcoRI, BamHI, HindIII, NotI, SmaI, XhoI

st.set_page_config(page_title="DNA Structural Analyzer", layout="wide", initial_sidebar_state="expanded")

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:ital,opsz,wght@0,9..40,300;0,9..40,400;0,9..40,500;0,9..40,600;1,9..40,400&family=DM+Mono:ital,wght@0,300;0,400;0,500;1,400&display=swap');

/* ═══════════════════════════════════════════════════════
   GEMINI × NOTEBOOKLM CROSSOVER THEME
   Light: warm NotebookLM parchment + Gemini gradient accents
   Dark:  warm NotebookLM ink + Gemini luminous gradient accents
═══════════════════════════════════════════════════════ */

:root {
  /* NotebookLM warm parchment surfaces */
  --bg:            #FAF7F2;
  --bg-card:       #FFFEF9;
  --bg-sidebar:    #F3EFE7;
  --bg-hover:      #EDE8DF;
  --bg-input:      #FFFEF9;
  --bg-code:       #F0EBE1;

  /* Borders: warm not cool */
  --border:        #DDD8CE;
  --border-focus:  #8B5CF6;

  /* Typography: ink-on-parchment */
  --text:          #1C1917;
  --text-muted:    #57534E;
  --text-faint:    #A8A29E;

  /* Gemini gradient accent — the signature iridescent sweep */
  --accent:        #7C3AED;          /* deep violet anchor */
  --accent-mid:    #2563EB;          /* electric blue midpoint */
  --accent-teal:   #0891B2;          /* cyan-teal tail */
  --grad:          linear-gradient(115deg, #7C3AED 0%, #2563EB 45%, #0891B2 100%);
  --grad-subtle:   linear-gradient(115deg,
                     rgba(124,58,237,0.12) 0%,
                     rgba(37,99,235,0.10) 45%,
                     rgba(8,145,178,0.08) 100%);
  --grad-glow:     linear-gradient(115deg,
                     rgba(124,58,237,0.25) 0%,
                     rgba(37,99,235,0.20) 45%,
                     rgba(8,145,178,0.15) 100%);

  /* Semantic bio colors — warm analogues */
  --stable:        #059669;   /* emerald */
  --partial:       #2563EB;   /* blue */
  --unstable:      #DC2626;   /* red */
  --enzyme:        #D97706;   /* amber */

  /* Elevation — warm shadows */
  --shadow:        0 1px 3px rgba(28,25,23,0.08), 0 1px 2px rgba(28,25,23,0.06);
  --shadow-md:     0 3px 8px rgba(28,25,23,0.10), 0 1px 4px rgba(28,25,23,0.07);
  --shadow-lg:     0 8px 24px rgba(28,25,23,0.12), 0 2px 8px rgba(28,25,23,0.08);

  --radius:        10px;
  --radius-sm:     6px;
  --radius-pill:   999px;
  --font:          'DM Sans', sans-serif;
  --font-mono:     'DM Mono', 'Roboto Mono', monospace;
}

@media (prefers-color-scheme: dark) {
  :root {
    /* NotebookLM warm ink surfaces — not cool gray, warm dark */
    --bg:            #18161A;
    --bg-card:       #221F25;
    --bg-sidebar:    #1D1A20;
    --bg-hover:      #2E2A33;
    --bg-input:      #221F25;
    --bg-code:       #1D1A20;

    --border:        #3D3843;
    --border-focus:  #A78BFA;

    --text:          #EDE9E3;
    --text-muted:    #B8B0AA;
    --text-faint:    #7A7470;

    /* Gemini luminous gradient — brighter in dark, like glowing stars */
    --accent:        #A78BFA;
    --accent-mid:    #60A5FA;
    --accent-teal:   #22D3EE;
    --grad:          linear-gradient(115deg, #A78BFA 0%, #60A5FA 45%, #22D3EE 100%);
    --grad-subtle:   linear-gradient(115deg,
                       rgba(167,139,250,0.18) 0%,
                       rgba(96,165,250,0.15) 45%,
                       rgba(34,211,238,0.12) 100%);
    --grad-glow:     linear-gradient(115deg,
                       rgba(167,139,250,0.35) 0%,
                       rgba(96,165,250,0.28) 45%,
                       rgba(34,211,238,0.20) 100%);

    --stable:        #34D399;
    --partial:       #60A5FA;
    --unstable:      #F87171;
    --enzyme:        #FBBF24;

    --shadow:        0 1px 4px rgba(0,0,0,0.35), 0 1px 2px rgba(0,0,0,0.25);
    --shadow-md:     0 3px 10px rgba(0,0,0,0.40), 0 1px 4px rgba(0,0,0,0.30);
    --shadow-lg:     0 8px 28px rgba(0,0,0,0.50), 0 2px 8px rgba(0,0,0,0.35);
  }
}

/* ── Base ── */
html, body, .stApp {
  background-color: var(--bg) !important;
  color: var(--text) !important;
  font-family: var(--font) !important;
  font-size: 14px;
  line-height: 1.55;
}

/* ── Global font + color catch-all: override every Streamlit element ── */
.stApp *, .stApp *::before, .stApp *::after {
  font-family: var(--font) !important;
  box-sizing: border-box;
}
/* Restore mono for elements that need it */
.stApp pre, .stApp code, .stApp kbd,
[data-testid="stCodeBlock"] *,
.stTextArea textarea, .stTextInput input {
  font-family: var(--font-mono) !important;
}
/* Kill any remaining white/light text from Streamlit defaults */
.stApp p, .stApp span, .stApp label, .stApp div,
.stApp li, .stApp td, .stApp th, .stApp small,
[data-testid="stMarkdownContainer"] *,
[data-testid="stText"] *,
.stRadio label span, .stCheckbox label span,
[data-baseweb="typo-labelsmall"],
[data-baseweb="typo-labelmedium"],
[data-baseweb="typo-paragraphsmall"],
[data-baseweb="typo-paragraphmedium"] {
  color: var(--text) !important;
}
/* Muted roles */
.stApp label,
[data-testid="stWidgetLabel"] p,
[data-testid="stWidgetLabel"] span {
  color: var(--text-muted) !important;
  font-size: 0.8rem !important;
  font-weight: 500 !important;
}
/* Metric label faint */
[data-testid="stMetricLabel"] p { color: var(--text-faint) !important; }

/* ── Sidebar ── */
[data-testid="stSidebar"] {
  background-color: var(--bg-sidebar) !important;
  border-right: 1px solid var(--border) !important;
}
[data-testid="stSidebarContent"] { padding: 1.25rem 1rem; }
[data-testid="stSidebar"] label,
[data-testid="stSidebar"] p,
[data-testid="stSidebar"] span {
  color: var(--text-muted) !important;
  font-size: 0.8rem !important;
  font-family: var(--font) !important;
}

/* ── Page title — Gemini gradient text ── */
h1 {
  font-family: var(--font) !important;
  font-size: 1.4rem !important;
  font-weight: 600 !important;
  background: var(--grad) !important;
  -webkit-background-clip: text !important;
  -webkit-text-fill-color: transparent !important;
  background-clip: text !important;
  letter-spacing: -0.01em !important;
  margin-bottom: 0 !important;
}

/* ── HR as gradient line ── */
hr {
  border: none !important;
  height: 1px !important;
  background: var(--grad) !important;
  opacity: 0.25 !important;
  margin: 1rem 0 !important;
}

/* ── Tabs ── */
[data-testid="stTabs"] button {
  font-family: var(--font) !important;
  font-size: 0.83rem !important;
  font-weight: 500 !important;
  color: var(--text-muted) !important;
  padding: 0.5rem 1.2rem !important;
  border: none !important;
  background: none !important;
  border-radius: var(--radius) var(--radius) 0 0 !important;
  transition: color 0.15s, background 0.15s !important;
}
[data-testid="stTabs"] button:hover {
  color: var(--accent) !important;
  background: var(--grad-subtle) !important;
}
[data-testid="stTabs"] button[aria-selected="true"] {
  color: var(--accent) !important;
  font-weight: 600 !important;
  background: var(--grad-subtle) !important;
  border-bottom: 2px solid transparent !important;
  border-image: var(--grad) 1 !important;
}
[data-testid="stTabs"] [role="tablist"] {
  border-bottom: 1px solid var(--border) !important;
  gap: 0.2rem !important;
}

/* ── Metric cards ── */
[data-testid="stMetric"] {
  background: var(--bg-card) !important;
  border: 1px solid var(--border) !important;
  border-radius: var(--radius) !important;
  padding: 1rem 1.25rem !important;
  box-shadow: var(--shadow) !important;
  position: relative !important;
  overflow: hidden !important;
  transition: box-shadow 0.2s, transform 0.2s !important;
}
[data-testid="stMetric"]::before {
  content: '' !important;
  position: absolute !important;
  top: 0; left: 0; right: 0 !important;
  height: 2px !important;
  background: var(--grad) !important;
}
[data-testid="stMetric"]:hover {
  box-shadow: var(--shadow-md) !important;
  transform: translateY(-1px) !important;
}
[data-testid="stMetricLabel"] p {
  font-family: var(--font) !important;
  font-size: 0.71rem !important;
  font-weight: 500 !important;
  text-transform: uppercase !important;
  letter-spacing: 0.07em !important;
  color: var(--text-faint) !important;
}
[data-testid="stMetricValue"] {
  font-family: var(--font) !important;
  font-size: 1.55rem !important;
  font-weight: 600 !important;
  background: var(--grad) !important;
  -webkit-background-clip: text !important;
  -webkit-text-fill-color: transparent !important;
  background-clip: text !important;
}

/* ── Buttons ── */
.stButton > button {
  font-family: var(--font) !important;
  font-size: 0.82rem !important;
  font-weight: 500 !important;
  background: var(--bg-card) !important;
  border: 1px solid var(--border) !important;
  color: var(--accent) !important;
  border-radius: var(--radius-pill) !important;
  padding: 0.35rem 1rem !important;
  box-shadow: var(--shadow) !important;
  transition: all 0.18s ease !important;
  text-transform: none !important;
  letter-spacing: 0 !important;
}
.stButton > button:hover {
  background: var(--grad-subtle) !important;
  box-shadow: var(--shadow-md) !important;
  border-color: var(--accent) !important;
}
[data-testid="stDownloadButton"] button {
  width: 100% !important;
  border-radius: var(--radius) !important;
  background: var(--grad) !important;
  color: #fff !important;
  border: none !important;
  font-weight: 600 !important;
  box-shadow: var(--shadow-md) !important;
  letter-spacing: 0.01em !important;
}
[data-testid="stDownloadButton"] button:hover {
  filter: brightness(1.08) saturate(1.1) !important;
  box-shadow: var(--shadow-lg) !important;
}

/* ── Multiselect ── */
[data-testid="stMultiSelect"] { font-size: 0.82rem !important; }

/* ── Inputs — text area, text input ── */
.stTextArea textarea, .stTextInput input {
  font-family: var(--font-mono) !important;
  font-size: 0.82rem !important;
  background-color: var(--bg-input) !important;
  border: 1px solid var(--border) !important;
  color: var(--text) !important;
  border-radius: var(--radius-sm) !important;
  caret-color: var(--accent) !important;
}
.stTextArea textarea::placeholder, .stTextInput input::placeholder {
  color: var(--text-faint) !important;
  opacity: 1 !important;
}
.stTextArea textarea:focus, .stTextInput input:focus {
  border-color: var(--border-focus) !important;
  box-shadow: 0 0 0 3px rgba(124,58,237,0.12) !important;
  outline: none !important;
}
/* Text selection highlight */
.stApp ::selection { background: rgba(124,58,237,0.2); color: var(--text); }
.stApp ::-moz-selection { background: rgba(124,58,237,0.2); color: var(--text); }

/* ── Cursor ── */
.stApp { cursor: default; }
.stApp a, .stApp button, .stApp [role="button"],
.stApp [data-baseweb="tab"], .stApp summary,
.stApp [data-testid="stFileUploadDropzone"] { cursor: pointer !important; }
.stApp textarea, .stApp input[type="text"],
.stApp input[type="search"] { cursor: text !important; }
.stApp [data-baseweb="slider"] { cursor: grab !important; }
.stApp [data-baseweb="slider"]:active { cursor: grabbing !important; }

/* ── Scrollbars ── */
.stApp ::-webkit-scrollbar { width: 6px; height: 6px; }
.stApp ::-webkit-scrollbar-track { background: var(--bg); }
.stApp ::-webkit-scrollbar-thumb {
  background: var(--border);
  border-radius: var(--radius-pill);
}
.stApp ::-webkit-scrollbar-thumb:hover { background: var(--text-faint); }

/* ── Selectbox + Multiselect — every layer ── */

/* Outer container */
[data-testid="stSelectbox"] > div > div,
[data-testid="stMultiSelect"] > div > div {
  background: var(--bg-input) !important;
  border: 1px solid var(--border) !important;
  border-radius: var(--radius-sm) !important;
}
/* The control row itself */
[data-baseweb="select"] > div:first-child {
  background: var(--bg-input) !important;
  border: 1px solid var(--border) !important;
  border-radius: var(--radius-sm) !important;
}
[data-baseweb="select"]:focus-within > div:first-child {
  border-color: var(--border-focus) !important;
  box-shadow: 0 0 0 3px rgba(124,58,237,0.12) !important;
}
/* Text inside select control */
[data-baseweb="select"] [data-baseweb="single-value"],
[data-baseweb="select"] [data-baseweb="placeholder"],
[data-baseweb="select"] input,
[data-baseweb="select"] input::placeholder {
  color: var(--text) !important;
  font-family: var(--font) !important;
  font-size: 0.82rem !important;
  caret-color: var(--accent) !important;
}
[data-baseweb="select"] [data-baseweb="placeholder"] {
  color: var(--text-faint) !important;
}
/* Chevron icon */
[data-baseweb="select"] svg,
[data-baseweb="select"] [data-baseweb="select-arrow"] svg {
  fill: var(--text-faint) !important;
  color: var(--text-faint) !important;
}
/* Tags/chips inside multiselect */
[data-baseweb="tag"] {
  background: var(--grad-subtle) !important;
  border: 1px solid var(--border-focus) !important;
  color: var(--accent) !important;
  border-radius: var(--radius-pill) !important;
  font-family: var(--font) !important;
  font-size: 0.75rem !important;
}
[data-baseweb="tag"] span { color: var(--accent) !important; }
[data-baseweb="tag"] [data-baseweb="tag-action"] svg { fill: var(--accent) !important; }

/* Dropdown popover panel */
[data-baseweb="popover"],
[data-baseweb="popover"] > div,
[data-baseweb="menu"],
[data-baseweb="menu"] ul,
[data-baseweb="select"] [role="listbox"] {
  background: var(--bg-card) !important;
  border: 1px solid var(--border) !important;
  border-radius: var(--radius) !important;
  box-shadow: var(--shadow-lg) !important;
  padding: 0.25rem !important;
}
/* Each option row */
[data-baseweb="option"],
[data-baseweb="menu"] li {
  background: var(--bg-card) !important;
  color: var(--text-muted) !important;
  font-family: var(--font) !important;
  font-size: 0.82rem !important;
  border-radius: var(--radius-sm) !important;
  padding: 0.4rem 0.75rem !important;
  transition: background 0.12s, color 0.12s !important;
}
[data-baseweb="option"]:hover {
  background: var(--grad-subtle) !important;
  color: var(--accent) !important;
}
[data-baseweb="option"][aria-selected="true"] {
  background: var(--grad-subtle) !important;
  color: var(--accent) !important;
  font-weight: 600 !important;
}
/* Disabled option */
[data-baseweb="option"][aria-disabled="true"] {
  color: var(--text-faint) !important;
  cursor: not-allowed !important;
}

/* ── Slider ── */
[data-testid="stSlider"] [data-baseweb="slider"] [role="slider"] {
  background: var(--accent) !important;
  border-color: var(--accent) !important;
  box-shadow: 0 0 0 4px rgba(124,58,237,0.15) !important;
}
/* Track fill */
[data-testid="stSlider"] [data-baseweb="slider"] div[class*="Track"] > div:first-child {
  background: var(--grad) !important;
}
/* Track background */
[data-testid="stSlider"] [data-baseweb="slider"] div[class*="Track"] {
  background: var(--border) !important;
}
/* Tick labels */
[data-testid="stSlider"] [data-testid="stMarkdownContainer"] p,
[data-testid="stSlider"] span {
  color: var(--text-muted) !important;
  font-family: var(--font) !important;
}

/* ── Checkboxes & Radios ── */
[data-baseweb="checkbox"] [data-checked="true"] div,
[data-baseweb="radio"] [data-checked="true"] div {
  background: var(--accent) !important;
  border-color: var(--accent) !important;
}
[data-baseweb="checkbox"] div,
[data-baseweb="radio"] div {
  border-color: var(--border) !important;
  background: var(--bg-input) !important;
}
[data-baseweb="checkbox"] label span,
[data-baseweb="radio"] label span {
  color: var(--text-muted) !important;
  font-family: var(--font) !important;
}

/* ── Number input ── */
[data-baseweb="input"] {
  background: var(--bg-input) !important;
  border: 1px solid var(--border) !important;
  border-radius: var(--radius-sm) !important;
}
[data-baseweb="input"] input {
  background: transparent !important;
  color: var(--text) !important;
  font-family: var(--font-mono) !important;
  caret-color: var(--accent) !important;
}
[data-baseweb="input"] input::placeholder { color: var(--text-faint) !important; }
[data-baseweb="input"]:focus-within {
  border-color: var(--border-focus) !important;
  box-shadow: 0 0 0 3px rgba(124,58,237,0.12) !important;
}

/* ── Tooltip ── */
[data-baseweb="tooltip"] div {
  background: var(--text) !important;
  color: var(--bg) !important;
  font-family: var(--font) !important;
  font-size: 0.75rem !important;
  border-radius: var(--radius-sm) !important;
  padding: 0.3rem 0.6rem !important;
}

/* ── Streamlit alert/info/warning/error boxes ── */
[data-testid="stAlert"] {
  background: var(--bg-card) !important;
  border: 1px solid var(--border) !important;
  border-radius: var(--radius) !important;
  color: var(--text-muted) !important;
}
[data-testid="stAlert"] p, [data-testid="stAlert"] span {
  color: var(--text-muted) !important;
}

/* ── Section labels — warm small caps ── */
.section-label {
  font-family: var(--font) !important;
  font-size: 0.69rem;
  font-weight: 600;
  text-transform: uppercase;
  letter-spacing: 0.09em;
  color: var(--text-faint);
  margin-bottom: 0.65rem;
  padding-bottom: 0.32rem;
  border-bottom: 1px solid var(--border);
}

/* ── Cards ── */
.diff-card {
  background: var(--bg-card);
  border: 1px solid var(--border);
  padding: 1.1rem 1.25rem;
  border-radius: var(--radius);
  font-family: var(--font-mono);
  font-size: 0.8rem;
  line-height: 2.2;
  white-space: pre-wrap;
  word-break: break-all;
  color: var(--text-muted);
  box-shadow: var(--shadow);
}
.diff-del {
  color: var(--unstable);
  text-decoration: line-through;
  background: color-mix(in srgb, var(--unstable) 10%, transparent);
  padding: 0 3px;
  border-radius: 3px;
}
.diff-add {
  color: var(--stable);
  font-weight: 600;
  background: color-mix(in srgb, var(--stable) 10%, transparent);
  padding: 0 3px;
  border-radius: 3px;
}
.diag-card {
  background: var(--bg-card);
  border: 1px solid var(--border);
  padding: 1.1rem 1.25rem;
  border-radius: var(--radius);
  font-family: var(--font-mono);
  font-size: 0.8rem;
  white-space: pre-wrap;
  word-break: break-all;
  box-shadow: var(--shadow);
}

/* ── Primer box — parchment scroll feel ── */
.primer-box {
  background: var(--bg-code);
  border-left: 3px solid transparent;
  border-image: var(--grad) 1;
  padding: 0.85rem 1rem;
  font-family: var(--font-mono);
  font-size: 0.78rem;
  color: var(--text-muted);
  line-height: 1.9;
  word-break: break-all;
  border-radius: 0 var(--radius-sm) var(--radius-sm) 0;
}

/* ── Legend panel ── */
.legend-panel {
  background: var(--bg-card);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  padding: 1rem 1.1rem;
  box-shadow: var(--shadow);
}
.legend-row {
  display: flex;
  align-items: center;
  gap: 0.7rem;
  margin-bottom: 0.55rem;
  font-family: var(--font);
  font-size: 0.78rem;
  color: var(--text-muted);
}
.legend-swatch { width: 24px; height: 4px; border-radius: 2px; flex-shrink: 0; }
.legend-dot    { width: 8px; height: 8px; border-radius: 50%; flex-shrink: 0; }

.mut-row {
  font-family: var(--font);
  font-size: 0.82rem;
  padding: 0.5rem 0;
  border-bottom: 1px solid var(--border);
  color: var(--text-muted);
}
.focus-panel {
  background: var(--bg-card);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  padding: 0.9rem 1.1rem;
  font-family: var(--font-mono);
  font-size: 0.78rem;
  color: var(--text-muted);
  line-height: 2.0;
  box-shadow: var(--shadow);
}

.stAlert {
  background: var(--bg-card) !important;
  border: 1px solid var(--border) !important;
  border-radius: var(--radius) !important;
}

/* ── st.code / code blocks ── */
[data-testid="stCodeBlock"],
[data-testid="stCodeBlock"] pre,
[data-testid="stCodeBlock"] code {
  background: var(--bg-code) !important;
  color: var(--text-muted) !important;
  border: 1px solid var(--border) !important;
  border-radius: var(--radius-sm) !important;
  font-family: var(--font-mono) !important;
  font-size: 0.8rem !important;
}

/* ── File uploader ── */
[data-testid="stFileUploader"] > div,
[data-testid="stFileUploadDropzone"] {
  background: var(--bg-card) !important;
  border: 1px dashed var(--border) !important;
  border-radius: var(--radius) !important;
  color: var(--text-muted) !important;
}
[data-testid="stFileUploadDropzone"]:hover {
  background: var(--grad-subtle) !important;
  border-color: var(--border-focus) !important;
}
[data-testid="stFileUploadDropzone"] span,
[data-testid="stFileUploadDropzone"] p,
[data-testid="stFileUploadDropzone"] small {
  color: var(--text-muted) !important;
}
[data-testid="stFileUploadDropzone"] svg {
  fill: var(--text-faint) !important;
}


#MainMenu, footer, header { visibility: hidden; }

/* ── Restriction Enzyme reference panel ── */
.re-ref-panel {
  position: fixed;
  top: 3.5rem;
  right: 1rem;
  z-index: 9999;
  width: 320px;
  background: var(--bg-card);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  font-family: var(--font);
  font-size: 0.78rem;
  box-shadow: var(--shadow-lg);
  overflow: hidden;
}
.re-ref-header {
  display: flex;
  align-items: center;
  justify-content: space-between;
  padding: 0.6rem 0.9rem;
  background: var(--grad-subtle);
  border-bottom: 1px solid var(--border);
  cursor: pointer;
  user-select: none;
  font-size: 0.76rem;
  font-weight: 600;
  color: var(--accent);
  list-style: none;
  letter-spacing: 0.02em;
}
.re-ref-header::-webkit-details-marker { display: none; }
.re-ref-header::after {
  content: 'be';
  font-size: 0.7rem;
  color: var(--text-faint);
  transition: transform 0.15s ease;
}
details[open] .re-ref-header::after { transform: rotate(180deg); }
.re-ref-header:hover { filter: brightness(1.04); }
.re-ref-body { padding: 0.5rem 0.6rem 0.7rem; }
.re-ref-table { width: 100%; border-collapse: collapse; }
.re-ref-table th {
  color: var(--text-faint);
  text-transform: uppercase;
  letter-spacing: 0.07em;
  font-size: 0.65rem;
  font-weight: 600;
  padding: 0.3rem 0.5rem;
  border-bottom: 1px solid var(--border);
  text-align: left;
}
.re-ref-table td {
  padding: 0.32rem 0.5rem;
  color: var(--text-muted);
  border-bottom: 1px solid var(--border);
  font-size: 0.75rem;
  vertical-align: middle;
}
.re-ref-table tr:last-child td { border-bottom: none; }
.re-ref-table tr:hover td { background: var(--grad-subtle); color: var(--text); }
.re-enz-name   { color: var(--enzyme) !important; font-weight: 600; }
.re-site-seq   { color: var(--stable) !important; font-family: var(--font-mono) !important; letter-spacing: 0.06em; }
.re-cut-sticky { color: var(--partial) !important; font-size: 0.71rem; }
.re-cut-blunt  { color: var(--text-faint) !important; font-size: 0.71rem; }
</style>
""", unsafe_allow_html=True)
# Inject JS to sync CSS vars with system theme changes in real-time
st.markdown("""
<script>
(function() {
  function applyTheme(dark) {
    var root = document.documentElement;
    if (dark) {
      root.style.setProperty('--bg',           '#18161A');
      root.style.setProperty('--bg-card',      '#221F25');
      root.style.setProperty('--bg-sidebar',   '#1D1A20');
      root.style.setProperty('--bg-hover',     '#2E2A33');
      root.style.setProperty('--bg-input',     '#221F25');
      root.style.setProperty('--bg-code',      '#1D1A20');
      root.style.setProperty('--border',       '#3D3843');
      root.style.setProperty('--border-focus', '#A78BFA');
      root.style.setProperty('--text',         '#EDE9E3');
      root.style.setProperty('--text-muted',   '#B8B0AA');
      root.style.setProperty('--text-faint',   '#7A7470');
      root.style.setProperty('--accent',       '#A78BFA');
      root.style.setProperty('--accent-mid',   '#60A5FA');
      root.style.setProperty('--accent-teal',  '#22D3EE');
      root.style.setProperty('--grad',         'linear-gradient(115deg,#A78BFA 0%,#60A5FA 45%,#22D3EE 100%)');
      root.style.setProperty('--grad-subtle',  'linear-gradient(115deg,rgba(167,139,250,0.18) 0%,rgba(96,165,250,0.15) 45%,rgba(34,211,238,0.12) 100%)');
      root.style.setProperty('--grad-glow',    'linear-gradient(115deg,rgba(167,139,250,0.35) 0%,rgba(96,165,250,0.28) 45%,rgba(34,211,238,0.20) 100%)');
      root.style.setProperty('--stable',       '#34D399');
      root.style.setProperty('--partial',      '#60A5FA');
      root.style.setProperty('--unstable',     '#F87171');
      root.style.setProperty('--enzyme',       '#FBBF24');
      root.style.setProperty('--shadow',       '0 1px 4px rgba(0,0,0,0.35),0 1px 2px rgba(0,0,0,0.25)');
      root.style.setProperty('--shadow-md',    '0 3px 10px rgba(0,0,0,0.40),0 1px 4px rgba(0,0,0,0.30)');
      root.style.setProperty('--shadow-lg',    '0 8px 28px rgba(0,0,0,0.50),0 2px 8px rgba(0,0,0,0.35)');
    } else {
      var props = ['--bg','--bg-card','--bg-sidebar','--bg-hover','--bg-input','--bg-code',
                   '--border','--border-focus','--text','--text-muted','--text-faint',
                   '--accent','--accent-mid','--accent-teal','--grad','--grad-subtle','--grad-glow',
                   '--stable','--partial','--unstable','--enzyme',
                   '--shadow','--shadow-md','--shadow-lg'];
      props.forEach(function(p){ root.style.removeProperty(p); });
    }
  }
  var mq = window.matchMedia('(prefers-color-scheme: dark)');
  applyTheme(mq.matches);
  mq.addEventListener('change', function(e) { applyTheme(e.matches); });
})();
</script>
""", unsafe_allow_html=True)


C_STABLE   = "#137333"
C_PARTIAL  = "#1A73E8"
C_UNSTABLE = "#C5221F"
C_ENZYME   = "#E37400"

RB = RestrictionBatch([EcoRI, BamHI, HindIII, NotI, SmaI, XhoI])

GRNA_POS_WEIGHTS = [
    0.000, 0.000, 0.014, 0.000, 0.000, 0.395, 0.317, 0.000, 0.389, 0.079,
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583
]

MISMATCH_PENALTY = {
    frozenset(["A","G"]): 0.4,   # transition
    frozenset(["C","T"]): 0.4,   # transition
    frozenset(["A","C"]): 0.7,   # transversion
    frozenset(["A","T"]): 0.6,   # transversion
    frozenset(["G","C"]): 0.7,   # transversion
    frozenset(["G","T"]): 0.6,   # transversion
}

def find_grna_candidates(seq):
    """
    Scan both strands for NGG PAM motifs and extract the upstream 20 bp protospacer.
    Uses Bio.Seq for strand handling and reverse complement generation.
    Returns list of dicts with: grna, pam_idx, strand, gc_pct.
    Minimum sequence length required: 23 bp (20 bp spacer + 3 bp PAM).
    """
    candidates = []
    bio_seq = Seq(seq)
    strands = [("+", str(bio_seq)), ("-", str(bio_seq.reverse_complement()))]

    for strand, s in strands:
        for i in range(len(s) - 22):
            pam = s[i+20:i+23]
            if pam[1] == "G" and pam[2] == "G":
                grna = s[i:i+20]
                if len(grna) == 20:
                    gc = gc_fraction(Seq(grna)) * 100
                    pam_idx = i + 20 if strand == "+" else len(seq) - (i + 23)
                    candidates.append({
                        "grna":    grna,
                        "pam_idx": pam_idx,
                        "strand":  strand,
                        "gc_pct":  round(gc, 1),
                    })
    return candidates


def score_on_target(grna):
    """
    Simplified on-target efficacy score (0-100).
    Uses Bio.SeqUtils.gc_fraction for GC calculation and Bio.SeqUtils.MeltingTemp.Tm_NN
    for seed region thermodynamic stability assessment.
    """
    gc = gc_fraction(Seq(grna))

    gc_score = max(0.0, 1.0 - abs(gc - 0.55) * 2.5)

    poly_penalty = 0.0
    for base in "ATGC":
        if base * 4 in grna:
            poly_penalty = max(poly_penalty, 0.35)
        elif base * 3 in grna:
            poly_penalty = max(poly_penalty, 0.15)

    if "TTTT" in grna:
        poly_penalty = max(poly_penalty, 0.5)

    seed = grna[8:]
    try:
        seed_tm = Tm_NN(Seq(seed), nn_table=None)
        seed_score = max(0.0, min(1.0, seed_tm / 60.0))
    except Exception:
        seed_gc = gc_fraction(Seq(seed))
        seed_score = max(0.0, 1.0 - abs(seed_gc - 0.5) * 2.0)

    raw = gc_score * 0.45 + seed_score * 0.35 - poly_penalty * 0.20
    return round(min(max(raw, 0.0), 1.0) * 100, 1)


def score_off_target(grna, mock_genome_seqs):
    """
    MIT-style off-target risk score (0–100, lower = safer).
    Compares the gRNA against each sequence in mock_genome_seqs.
    The penalty per off-target hit = product of position-weighted mismatch scores.
    Aggregated with diminishing returns (sqrt sum).
    """
    if not mock_genome_seqs:
        return 0.0

    total_penalty = 0.0
    for candidate in mock_genome_seqs:
        if len(candidate) < 20:
            continue
        aln = candidate[:20]
        hit_score = 1.0
        for pos in range(20):
            g_base = grna[pos]
            t_base = aln[pos]
            if g_base != t_base:
                key = frozenset([g_base, t_base])
                mm_penalty = MISMATCH_PENALTY.get(key, 0.7)
                hit_score *= (1.0 - GRNA_POS_WEIGHTS[pos] * mm_penalty)
        total_penalty += hit_score

    risk = min(total_penalty / len(mock_genome_seqs) * 100, 100.0)
    return round(risk, 1)


def predict_hairpin(grna, window=4):
    """
    Predict if the gRNA will form secondary structure (hairpin) that blocks Cas9.
    Uses Bio.Seq.reverse_complement for complement generation.
    Returns: (has_hairpin: bool, stem_length: int, description: str)
    """
    rc = str(Seq(grna).reverse_complement())
    max_stem = 0
    best_desc = ""
    L = len(grna)

    for i in range(L - window * 2):
        for j in range(i + window, L - window + 1):
            stem = 0
            while (i + stem < j) and (j + stem < L):
                if grna[i + stem] == rc[L - 1 - (j + stem)]:
                    stem += 1
                else:
                    break
            if stem >= window and stem > max_stem:
                max_stem = stem
                best_desc = f"stem {i}-{i+stem-1}  <->  {j}-{j+stem-1} ({stem} bp)"

    if max_stem >= 6:
        return True, max_stem, best_desc + "  strong hairpin"
    elif max_stem >= 4:
        return True, max_stem, best_desc + "  moderate hairpin"
    else:
        return False, max_stem, "No significant hairpin predicted"


def generate_mock_genome(seed_seq, n_decoys=30):
    """
    Generate a mock 'genome' for off-target scoring.
    Creates n_decoys sequences by applying 1–4 random mutations to the seed gRNA.
    This simulates plausible off-target loci with varying similarity.
    """
    bases = list("ATGC")
    decoys = []
    rng = random.Random(42)   # deterministic for reproducibility
    for _ in range(n_decoys):
        d = list(seed_seq)
        n_muts = rng.randint(1, 4)
        positions = rng.sample(range(20), n_muts)
        for p in positions:
            alts = [b for b in bases if b != d[p]]
            d[p] = rng.choice(alts)
        decoys.append("".join(d))
    return decoys


def parse_genome_file(uploaded_file):
    """
    Parse an uploaded genomic context file using Bio.SeqIO.
    Supports FASTA (.fa, .fasta, .fna) and plain text (.txt).
    Returns (cleaned_sequence, source_label).
    Only ATGC bases are retained; IUPAC ambiguity codes are stripped.
    """
    raw_bytes = uploaded_file.read()
    raw_text  = raw_bytes.decode("utf-8", errors="ignore")
    handle    = io.StringIO(raw_text)

    is_fasta = any(line.startswith(">") for line in raw_text.splitlines())

    if is_fasta:
        records = list(SeqIO.parse(handle, "fasta"))
        seq = "".join(str(r.seq).upper() for r in records)
        label = f"FASTA · {len(records)} record(s)"
    else:
        handle.seek(0)
        seq = "".join(line.strip().upper() for line in handle)
        label = "Plain text"

    seq = re.sub(r"[^ATGC]", "", seq)
    label += f" · {len(seq):,} bp loaded"
    return seq, label


def genome_to_decoys(genome_seq: str, grna: str, max_decoys: int = 200) -> list[str]:
    """
    Extract real off-target candidate windows from a genomic context sequence.
    Slides a 20-bp window across the genome and collects windows that sit
    upstream of an NGG PAM and have ≤4 mismatches with the guide.
    Falls back to the mock generator if no NGG windows are found.
    Returns a list of 20-bp candidate strings (may be empty).
    """
    candidates = []
    L = len(genome_seq)
    for i in range(L - 22):
        pam = genome_seq[i+20:i+23]
        if len(pam) == 3 and pam[1] == "G" and pam[2] == "G":
            window = genome_seq[i:i+20]
            if len(window) == 20:
                mm = sum(a != b for a, b in zip(grna, window))
                if mm <= 4:
                    candidates.append(window)
                    if len(candidates) >= max_decoys:
                        break
    return candidates if candidates else generate_mock_genome(grna, n_decoys=30)


def rank_grnas(candidates, n_decoys=30, genome_seq: str | None = None):
    """
    Score all candidates and return the top 3 by composite score.
    If genome_seq is provided, real off-target windows are extracted from it.
    Otherwise falls back to the mock genome generator.
    Composite = 0.55 * on_target - 0.30 * off_target_risk - 0.15 * hairpin_penalty
    """
    scored = []
    for c in candidates:
        grna = c["grna"]
        on   = score_on_target(grna)
        if genome_seq:
            decoys = genome_to_decoys(genome_seq, grna, max_decoys=n_decoys)
        else:
            decoys = generate_mock_genome(grna, n_decoys)
        off  = score_off_target(grna, decoys)
        hp, stem_len, hp_desc = predict_hairpin(grna)
        hp_penalty = stem_len * 5.0

        composite = round(on * 0.55 - off * 0.30 - hp_penalty * 0.15, 2)
        scored.append({**c,
            "on_target":  on,
            "off_target": off,
            "hairpin":    hp,
            "stem_len":   stem_len,
            "hp_desc":    hp_desc,
            "composite":  composite,
        })

    scored.sort(key=lambda x: x["composite"], reverse=True)
    return scored[:3]


def generate_crispr_pdf(top_grnas, gene_seq):
    """
    Generate a 'Ready-to-Order' PDF sheet for the top gRNA candidates.
    """
    import io
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.lib.enums import TA_LEFT, TA_JUSTIFY
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable
    from reportlab.lib import colors

    buf = io.BytesIO()
    W, H   = A4
    margin = 20 * mm

    doc = SimpleDocTemplate(buf, pagesize=A4,
                            leftMargin=margin, rightMargin=margin,
                            topMargin=16*mm, bottomMargin=16*mm)

    FG      = colors.HexColor("#1A2530")
    MUTED   = colors.HexColor("#4A6070")
    ACCENT  = colors.HexColor("#1A5A7A")
    RULE    = colors.HexColor("#C0CDD4")
    TBL_HDR = colors.HexColor("#E4EEF2")
    TBL_ROW = colors.HexColor("#F7FAFB")
    TBL_ALT = colors.HexColor("#EBF2F5")
    SAFE    = colors.HexColor("#1A5A4A")
    RISK    = colors.HexColor("#7A2A2A")
    MONO    = "Courier"
    SANS    = "Helvetica"

    def sty(name, font=SANS, size=9, color=FG, leading=14, align=TA_LEFT, **kw):
        return ParagraphStyle(name, fontName=font, fontSize=size,
                              textColor=color, leading=leading, alignment=align, **kw)

    S_title = sty("t", MONO, 13, ACCENT, 18, spaceAfter=2, letterSpacing=2)
    S_sub   = sty("s", MONO,  7, MUTED,   9, spaceAfter=0, letterSpacing=2)
    S_head  = sty("h", MONO,  8, ACCENT, 12, spaceBefore=6, spaceAfter=2, letterSpacing=1.5)
    S_body  = sty("b", SANS,  9, FG,     13, TA_JUSTIFY, spaceAfter=3)
    S_mono  = sty("m", MONO,  7, MUTED,  10, spaceAfter=2)
    S_cap   = sty("c", MONO,  7, MUTED,   9)

    def HR(): return HRFlowable(width="100%", thickness=0.5, color=RULE, spaceAfter=6, spaceBefore=4)
    def SP(h=3): return Spacer(1, h*mm)

    story = []
    story.append(Paragraph("CRISPR-Cas9 gRNA READY-TO-ORDER SHEET", S_title))
    story.append(Paragraph(
        f"Generated {datetime.now().strftime('%Y-%m-%d  %H:%M:%S')}  |  "
        f"Input length: {len(gene_seq)} bp  |  Top 3 candidates ranked by composite score",
        S_sub
    ))
    story.append(SP(2)); story.append(HR())

    story.append(Paragraph("SCORING METHODOLOGY", S_head))
    story.append(Paragraph(
        "On-target efficacy is estimated using GC content optimisation (40-70%), seed region "
        "thermodynamics, and poly-run penalties. Off-target risk is computed via a simplified "
        "MIT-style position-weighted mismatch scoring against a mock genome ensemble. "
        "Secondary structure risk uses reverse-complement self-alignment to detect hairpin stems "
        "that would block Cas9 loading. Composite = 0.55 x On-target - 0.30 x Off-risk - 0.15 x Hairpin.",
        S_body
    ))
    story.append(SP(1)); story.append(HR())

    for rank, g in enumerate(top_grnas, 1):
        story.append(Paragraph(f"CANDIDATE #{rank}  /  Composite: {g['composite']:.1f}", S_head))

        story.append(Paragraph("5'-NGG spacer sequence (order this 20-mer):", S_cap))
        story.append(Paragraph(g["grna"], S_mono))
        story.append(SP(1))

        hp_label = f"Yes ({g['stem_len']} bp stem)" if g["hairpin"] else "No"
        hp_color = RISK if g["hairpin"] else SAFE
        rows = [
            ["Metric", "Value", "Interpretation"],
            ["On-target score",  f"{g['on_target']:.1f} / 100",
             "Higher = better cleavage efficiency"],
            ["Off-target risk",  f"{g['off_target']:.1f} / 100",
             "Lower = fewer predicted off-target cuts"],
            ["GC content",       f"{g['gc_pct']:.1f}%",
             "Optimal range 40-70%"],
            ["Strand",           g["strand"],
             "+ forward  /  - reverse complement"],
            ["PAM index",        str(g["pam_idx"]),
             "Position in input sequence"],
            ["Hairpin risk",     hp_label,
             g["hp_desc"]],
        ]
        t = Table(rows, colWidths=[(W - 2*margin)/3]*3)
        t.setStyle(TableStyle([
            ("BACKGROUND",  (0,0), (-1,0),  TBL_HDR),
            ("BACKGROUND",  (0,1), (-1,-1), TBL_ROW),
            ("ROWBACKGROUNDS", (0,1), (-1,-1), [TBL_ROW, TBL_ALT]),
            ("TEXTCOLOR",   (0,0), (-1,0),  MUTED),
            ("TEXTCOLOR",   (0,1), (-1,-1), FG),
            ("TEXTCOLOR",   (1,5), (1,5),   hp_color),
            ("FONTNAME",    (0,0), (-1,-1), MONO),
            ("FONTSIZE",    (0,0), (-1,-1), 7.5),
            ("PADDING",     (0,0), (-1,-1), 4),
            ("GRID",        (0,0), (-1,-1), 0.4, RULE),
            ("BOX",         (0,0), (-1,-1), 0.8, MUTED),
        ]))
        story.append(t); story.append(SP(3))

    story.append(HR())
    story.append(Paragraph(
        "Scores are computational predictions only. Validate experimentally. "
        "Off-target analysis uses a stochastic mock-genome model and does not "
        "replace whole-genome in-silico alignment (e.g. Cas-OFFinder). "
        "Generated by DNA Structural Analyzer — research use only.",
        sty("foot", SANS, 7, MUTED, 10, TA_JUSTIFY)
    ))

    doc.build(story)
    buf.seek(0)
    return buf.read()

for k, v in [("dna",""), ("history",[]), ("selected_enzymes",[]), ("focus_idx",None)]:
    if k not in st.session_state: st.session_state[k] = v

def add_history(seq, label):
    if seq and (not st.session_state.history or st.session_state.history[0]["seq"] != seq):
        st.session_state.history.insert(0, {"time": datetime.now().strftime("%H:%M:%S"), "seq": seq, "label": label})

def force_snapshot(seq, label):
    if seq:
        st.session_state.history.insert(0, {"time": datetime.now().strftime("%H:%M:%S"), "seq": seq, "label": label})

def get_rev_comp(seq):
    return str(Seq(seq).reverse_complement())

def apply_mutation(seq, idx, base):
    s = list(seq); s[idx] = base; return "".join(s)

def generate_pdf_report(dna_seq, n, mean_tm, gc_pct, found_enz, mutations, rev_comp, gc5, gc3):
    """Generate a styled one-page A4 PDF report using reportlab."""
    import io
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.lib.enums import TA_LEFT, TA_JUSTIFY
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable
    from reportlab.lib import colors

    buf = io.BytesIO()
    W, H = A4
    margin = 20 * mm

    doc = SimpleDocTemplate(
        buf, pagesize=A4,
        leftMargin=margin, rightMargin=margin,
        topMargin=16*mm, bottomMargin=16*mm
    )

    FG      = colors.HexColor("#1A2530")    # dark slate - body text
    MUTED   = colors.HexColor("#4A6070")    # medium - labels
    ACCENT  = colors.HexColor("#1A5A7A")    # dark teal - headings/values
    STABLE  = colors.HexColor("#1A5A4A")
    UNSTAB  = colors.HexColor("#7A2A2A")
    ENZ_CLR = colors.HexColor("#7A5A10")
    RULE    = colors.HexColor("#C0CDD4")    # light rule lines
    TBL_HDR = colors.HexColor("#E4EEF2")    # table header bg
    TBL_ROW = colors.HexColor("#F7FAFB")    # table row bg
    TBL_ALT = colors.HexColor("#EBF2F5")    # alternate row tint
    MONO = "Courier"
    SANS = "Helvetica"

    def safe(s):
        """Sanitise string for ReportLab WinAnsiEncoding (Latin-1 safe)."""
        s = str(s)
        s = s.replace('°C', ' degC').replace('°', ' deg')
        s = s.replace('→', '->').replace('←', '<-')
        s = s.replace('–', '-').replace('—', ' - ')
        s = s.replace('‘', "'").replace('’', "'")
        s = s.replace('“', '"').replace('”', '"')
        s = s.replace('·', '.').replace('•', '*')
        s = s.replace('×', 'x').replace('≥', '>=').replace('≤', '<=')
        s = s.replace('é', 'e').replace('ö', 'o').replace('ü', 'u')
        s = s.replace('─', '-').replace('│', '|')
        return s.encode('latin-1', errors='replace').decode('latin-1')
    def sty(name, font=SANS, size=9, color=FG, leading=14, align=TA_LEFT, **kw):
        return ParagraphStyle(name, fontName=font, fontSize=size, textColor=color,
                              leading=leading, alignment=align, **kw)

    S_title    = sty("title",   MONO, 13, ACCENT,  18, spaceAfter=2, letterSpacing=2)
    S_sub      = sty("sub",     MONO,  7, MUTED,    9, spaceAfter=0, letterSpacing=2)
    S_head     = sty("head",    MONO,  8, ACCENT,  12, spaceBefore=6, spaceAfter=2, letterSpacing=1.5)
    S_body     = sty("body",    SANS,  9, FG,      13, TA_JUSTIFY, spaceAfter=3)
    S_mono     = sty("mono",    MONO,  7, MUTED,   10, spaceAfter=2)
    S_caption  = sty("caption", MONO,  7, MUTED,    9)

    def HR(): return HRFlowable(width="100%", thickness=0.5, color=RULE, spaceAfter=6, spaceBefore=4)
    def SP(h=3): return Spacer(1, h*mm)

    grade = "STABLE" if mean_tm > 18 else "UNSTABLE"
    grade_color = "#1A6A5A" if mean_tm > 18 else "#8A1A1A"

    if gc_pct < 40:
        gc_interp = (f"The GC content of {gc_pct:.1f}% is below the optimal 40–60% range, "
                     "suggesting AT-richness which may reduce thermostability and increase "
                     "susceptibility to denaturation under physiological conditions.")
    elif gc_pct > 60:
        gc_interp = (f"The GC content of {gc_pct:.1f}% exceeds the optimal 40–60% range. "
                     "This elevated GC fraction increases melting temperature but may also "
                     "promote secondary structure formation such as hairpin loops and G-quadruplexes.")
    else:
        gc_interp = (f"The GC content of {gc_pct:.1f}% falls within the optimal 40–60% range, "
                     "indicating a balanced nucleotide composition conducive to stable hybridisation "
                     "and reliable PCR amplification.")

    if len(mutations) == 0:
        mut_interp = ("No unstable windows were identified across the sequence. "
                      "All local melting temperatures meet the minimum threshold, "
                      "indicating uniform structural integrity along the strand.")
    else:
        mut_interp = (f"{len(mutations)} unstable region(s) were identified where the local melting "
                      f"temperature falls below 16 °C. These sites represent potential weak points "
                      f"under thermal or chemical stress and are candidates for targeted mutagenesis "
                      f"— specifically substituting A/T bases with G/C to increase local Tm.")

    if len(found_enz) == 0:
        enz_interp = ("No recognition sites for the six screened restriction endonucleases (EcoRI, "
                      "BamHI, HindIII, NotI, SmaI, XhoI) were detected. The sequence may be "
                      "introduced into restriction-based cloning vectors without risk of internal digestion.")
    else:
        enz_names = list({f['name'] for f in found_enz})
        enz_interp = (f"{len(found_enz)} restriction site(s) were mapped, involving "
                      f"{', '.join(enz_names)}. These sites must be considered when selecting a "
                      f"cloning strategy, as they represent cut points that could fragment the "
                      f"insert during vector preparation.")

    if mean_tm < 16:
        tm_interp = (f"The mean melting temperature of {mean_tm:.1f} °C is critically low, "
                     "indicating widespread thermodynamic instability. Significant sequence "
                     "re-engineering is recommended before proceeding with downstream applications.")
    elif mean_tm < 20:
        tm_interp = (f"The mean melting temperature of {mean_tm:.1f} °C is marginal. "
                     "While portions of the sequence are stable, AT-rich regions reduce overall "
                     "thermodynamic robustness. Selective G/C substitutions at flagged indices "
                     "are advisable.")
    else:
        tm_interp = (f"The mean melting temperature of {mean_tm:.1f} °C is within an acceptable "
                     "range for standard molecular biology applications including PCR, cloning, "
                     "and hybridisation-based assays.")

    story = []

    story.append(Paragraph("DNA STRUCTURAL ANALYSIS REPORT", S_title))
    story.append(Paragraph(
        f"Generated {datetime.now().strftime('%Y-%m-%d  %H:%M:%S')}  |  "
        f"Sequence length {n} bp  |  Grade: "
        f'<font color="{grade_color}"><b>{grade}</b></font>',
        S_sub
    ))
    story.append(SP(2)); story.append(HR())

    story.append(Paragraph("SEQUENCE METRICS", S_head))
    metrics_data = [
        ["Length", f"{n} bp",          "Mean Tm",    f"{mean_tm:.1f} degC"],
        ["GC Content", f"{gc_pct:.1f}%", "RE Sites",   str(len(found_enz))],
        ["Unstable Regions", str(len(mutations)), "GC Clamp 5p", f"{gc5}/5"],
        ["Rev-Comp Length", f"{n} bp",  "GC Clamp 3p", f"{gc3}/5"],
    ]
    t = Table(metrics_data, colWidths=[(W - 2*margin)/4]*4)
    t.setStyle(TableStyle([
        ("BACKGROUND",  (0,0), (-1,-1), TBL_ROW),
        ("BACKGROUND",  (0,0), (0,-1),  TBL_ALT),
        ("BACKGROUND",  (2,0), (2,-1),  TBL_ALT),
        ("TEXTCOLOR",   (0,0), (0,-1),  MUTED),
        ("TEXTCOLOR",   (2,0), (2,-1),  MUTED),
        ("TEXTCOLOR",   (1,0), (1,-1),  ACCENT),
        ("TEXTCOLOR",   (3,0), (3,-1),  ACCENT),
        ("FONTNAME",    (0,0), (-1,-1), MONO),
        ("FONTSIZE",    (0,0), (-1,-1), 8),
        ("PADDING",     (0,0), (-1,-1), 5),
        ("GRID",        (0,0), (-1,-1), 0.4, RULE),
        ("BOX",         (0,0), (-1,-1), 0.8, MUTED),
    ]))
    story.append(t); story.append(SP(3)); story.append(HR())

    story.append(Paragraph("ANALYTICAL COMMENTARY", S_head))
    story.append(Paragraph(
        "<b>Thermodynamic Assessment.</b>  " + safe(tm_interp), S_body))
    story.append(Paragraph(
        "<b>Nucleotide Composition.</b>  " + safe(gc_interp), S_body))
    story.append(Paragraph(
        "<b>Structural Stability.</b>  " + safe(mut_interp), S_body))
    story.append(Paragraph(
        "<b>Restriction Enzyme Compatibility.</b>  " + safe(enz_interp), S_body))
    story.append(SP(1)); story.append(HR())

    if mutations:
        story.append(Paragraph("FLAGGED UNSTABLE SITES", S_head))
        mut_rows = [["Index", "From", "Suggested", "Local Tm"]]
        for m in mutations:
            mut_rows.append([str(m["Idx"]), m["From"], m["To"], f"{m['Tm']} degC"])
        mt = Table(mut_rows, colWidths=[(W - 2*margin)/4]*4)
        mt.setStyle(TableStyle([
            ("BACKGROUND",  (0,0), (-1,0),  TBL_HDR),
            ("BACKGROUND",  (0,1), (-1,-1), TBL_ROW),
            ("ROWBACKGROUNDS", (0,1), (-1,-1), [TBL_ROW, TBL_ALT]),
            ("TEXTCOLOR",   (0,0), (-1,0),  MUTED),
            ("TEXTCOLOR",   (0,1), (-1,-1), FG),
            ("TEXTCOLOR",   (1,1), (1,-1),  colors.HexColor("#8A2020")),
            ("TEXTCOLOR",   (2,1), (2,-1),  colors.HexColor("#1A6040")),
            ("FONTNAME",    (0,0), (-1,-1), MONO),
            ("FONTSIZE",    (0,0), (-1,-1), 7.5),
            ("PADDING",     (0,0), (-1,-1), 4),
            ("GRID",        (0,0), (-1,-1), 0.4, RULE),
            ("BOX",         (0,0), (-1,-1), 0.8, MUTED),
        ]))
        story.append(mt); story.append(SP(2)); story.append(HR())

    if found_enz:
        story.append(Paragraph("RESTRICTION ENZYME MAP", S_head))
        enz_site_map = {str(e): e.site for e in RB}
        enz_rows = [["Enzyme", "Index", "Recognition Site"]]
        for fe in found_enz:
            enz_rows.append([fe["name"], str(fe["idx"]), enz_site_map.get(fe["name"], "")])
        et = Table(enz_rows, colWidths=[(W - 2*margin)/3]*3)
        et.setStyle(TableStyle([
            ("BACKGROUND",  (0,0), (-1,0),  TBL_HDR),
            ("BACKGROUND",  (0,1), (-1,-1), TBL_ROW),
            ("ROWBACKGROUNDS", (0,1), (-1,-1), [TBL_ROW, TBL_ALT]),
            ("TEXTCOLOR",   (0,0), (-1,0),  MUTED),
            ("TEXTCOLOR",   (0,1), (-1,-1), FG),
            ("TEXTCOLOR",   (2,1), (2,-1),  ENZ_CLR),
            ("FONTNAME",    (0,0), (-1,-1), MONO),
            ("FONTSIZE",    (0,0), (-1,-1), 7.5),
            ("PADDING",     (0,0), (-1,-1), 4),
            ("GRID",        (0,0), (-1,-1), 0.4, RULE),
            ("BOX",         (0,0), (-1,-1), 0.8, MUTED),
        ]))
        story.append(et); story.append(SP(2)); story.append(HR())

    story.append(Paragraph("RAW SEQUENCE DATA", S_head))
    def wrap60(s):
        return "  ".join(s[i:i+60] for i in range(0, len(s), 60))
    story.append(Paragraph("5' -> 3'  (forward)", S_caption))
    story.append(Paragraph(wrap60(dna_seq), S_mono))
    story.append(SP(1))
    story.append(Paragraph("3' -> 5'  (reverse complement)", S_caption))
    story.append(Paragraph(wrap60(rev_comp), S_mono))
    story.append(SP(2)); story.append(HR())

    story.append(Paragraph(
        "This report was generated by DNA Structural Analyzer. "
        "Melting temperatures are calculated using the nearest-neighbour approximation "
        "(4 degC per G/C, 2 degC per A/T) within the selected analysis window. "
        "Results are intended for research use only.",
        sty("footer", SANS, 7, MUTED, 10, TA_JUSTIFY)
    ))

    doc.build(story)
    buf.seek(0)
    return buf.read()

def _note(msg):
    return (f'<div style="font-family:var(--font);font-size:0.82rem;'
            f'color:var(--text-faint);padding:0.4rem 0">{msg}</div>')

st.sidebar.markdown('<div class="section-label">Sequence</div>', unsafe_allow_html=True)

if st.sidebar.button("Generate random · 60 bp"):
    st.session_state.dna = "".join(random.choice("ATGC") for _ in range(60))
    add_history(st.session_state.dna, "Random")
    st.session_state.focus_idx = None
    st.rerun()

dna_in = st.sidebar.text_area(
    "seq", value=st.session_state.dna, placeholder="5′ ATGC... 3′",
    height=90, label_visibility="collapsed"
).upper().strip()

if dna_in != st.session_state.dna:
    st.session_state.dna = dna_in
    add_history(dna_in, "Manual")

st.sidebar.markdown('<div class="section-label" style="margin-top:1.2rem">Parameters</div>', unsafe_allow_html=True)
window_size = st.sidebar.slider("Analysis window (bp)", 3, 10, 6)

st.sidebar.markdown('<div class="section-label" style="margin-top:1.2rem">Enzyme overlay</div>', unsafe_allow_html=True)

invalid_idx = [i for i, c in enumerate(dna_in) if c not in "ATGC"]
is_valid = len(dna_in) >= window_size and not invalid_idx

st.markdown("<h1>DNA Structural Analyzer</h1>", unsafe_allow_html=True)
st.markdown("<hr>", unsafe_allow_html=True)

if not is_valid and dna_in:
    highlighted = "".join(
        f'<span style="color:var(--unstable);font-weight:600;text-decoration:underline">{c}</span>'
        if i in invalid_idx else f'<span style="color:var(--text-muted)">{c}</span>'
        for i, c in enumerate(dna_in)
    )
    st.markdown(f'<div class="diag-card">{highlighted}</div>', unsafe_allow_html=True)
    if st.button("Clean sequence"):
        cleaned = re.sub(r"[^ATGC]", "", dna_in)
        st.session_state.dna = cleaned
        add_history(cleaned, "Cleaned")
        st.rerun()

elif is_valid:
    N = len(dna_in)
    bio_seq = Seq(dna_in)
    tm_list, mutations = [], []
    for i in range(N - window_size + 1):
        chunk = dna_in[i:i+window_size]
        tm = round(Tm_Wallace(Seq(chunk)), 1)
        status = "Stable" if tm > 20 else ("Partial" if tm >= 16 else "Unstable")
        tm_list.append({"idx": i, "tm": tm, "status": status})
        if tm < 16:
            mutations.append({"Idx": i, "From": dna_in[i], "To": "G" if dna_in[i] in "AT" else "C", "Tm": tm})

    rb_analysis = RB.search(bio_seq, linear=True)
    found_enz = []
    enz_name_map = {
        EcoRI: "EcoRI", BamHI: "BamHI", HindIII: "HindIII",
        NotI: "NotI", SmaI: "SmaI", XhoI: "XhoI",
    }
    for enzyme, positions in rb_analysis.items():
        name = enz_name_map.get(enzyme, str(enzyme))
        site_len = len(enzyme.site)
        for pos in positions:
            found_enz.append({"name": name, "idx": pos - 1, "len": site_len})

    unique_enz = sorted({f["name"] for f in found_enz})
    if unique_enz:
        valid_prev = [e for e in st.session_state.selected_enzymes if e in unique_enz]
        st.session_state.selected_enzymes = st.sidebar.multiselect(
            "sites", unique_enz, default=valid_prev, label_visibility="collapsed"
        )
    else:
        st.session_state.selected_enzymes = []
        st.sidebar.caption("No sites detected")

    enz_bits = set()
    for f in found_enz:
        if f["name"] in st.session_state.selected_enzymes:
            for b in range(f["idx"], f["idx"]+f["len"]): enz_bits.add(b)

    mean_tm = np.mean([t["tm"] for t in tm_list])
    gc_pct  = gc_fraction(bio_seq) * 100

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Length",     f"{N} bp")
    c2.metric("Mean Tm",    f"{mean_tm:.1f} degC")
    c3.metric("GC Content", f"{gc_pct:.1f}%")
    c4.metric("RE Sites",   len(found_enz))

    st.markdown("<div style='height:0.5rem'></div>", unsafe_allow_html=True)

    tabs = st.tabs(["Helix · 3D", "Engineering", "CRISPR · gRNA", "Comparison", "History"])

    with tabs[0]:
        st.markdown("""
<div class="re-ref-panel">
  <details>
    <summary class="re-ref-header">Restriction Enzymes</summary>
    <div class="re-ref-body">
      <table class="re-ref-table">
        <thead>
          <tr>
            <th>Enzyme</th>
            <th>Organism</th>
            <th>Recognition</th>
            <th>Cut Type</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td class="re-enz-name">EcoRI</td>
            <td><i>E. coli</i></td>
            <td class="re-site-seq">GAATTC</td>
            <td class="re-cut-sticky">Sticky 4-base</td>
          </tr>
          <tr>
            <td class="re-enz-name">BamHI</td>
            <td><i>B. amyloliquefaciens</i></td>
            <td class="re-site-seq">GGATCC</td>
            <td class="re-cut-sticky">Sticky 4-base</td>
          </tr>
          <tr>
            <td class="re-enz-name">HindIII</td>
            <td><i>H. influenzae</i></td>
            <td class="re-site-seq">AAGCTT</td>
            <td class="re-cut-sticky">Sticky 4-base</td>
          </tr>
          <tr>
            <td class="re-enz-name">HaeIII</td>
            <td><i>H. aegyptius</i></td>
            <td class="re-site-seq">GGCC</td>
            <td class="re-cut-blunt">Blunt</td>
          </tr>
          <tr>
            <td class="re-enz-name">NotI</td>
            <td><i>N. otitidiscaviarum</i></td>
            <td class="re-site-seq">GCGGCCGC</td>
            <td class="re-cut-sticky">Sticky 4-base</td>
          </tr>
        </tbody>
      </table>
    </div>
  </details>
</div>
""", unsafe_allow_html=True)
        col_v, col_l = st.columns([5, 2], gap="medium")

        with col_l:
            st.markdown("<div style='height:0.6rem'></div>", unsafe_allow_html=True)
            st.markdown('<div class="section-label">Stability key</div>', unsafe_allow_html=True)
            st.markdown(f"""
            <div class="legend-panel">
              <div class="legend-row">
                <div class="legend-swatch" style="background:{C_STABLE}"></div>
                <span>Stable &nbsp;<span style="color:var(--stable)">Tm &gt; 20 °C</span></span>
              </div>
              <div class="legend-row">
                <div class="legend-swatch" style="background:{C_PARTIAL}"></div>
                <span>Partial &nbsp;<span style="color:var(--partial)">16 – 20 °C</span></span>
              </div>
              <div class="legend-row">
                <div class="legend-swatch" style="background:{C_UNSTABLE}"></div>
                <span>Unstable &nbsp;<span style="color:var(--unstable)">Tm &lt; 16 °C</span></span>
              </div>
              <div class="legend-row">
                <div class="legend-swatch" style="background:{C_ENZYME};height:3px;border-radius:0"></div>
                <span>Enzyme site</span>
              </div>
              <div class="legend-row" style="margin-top:0.5rem">
                <div class="legend-swatch" style="background:rgba(50,100,120,0.7)"></div>
                <span>Strand 1</span>
              </div>
              <div class="legend-row">
                <div class="legend-swatch" style="background:rgba(40,70,100,0.7)"></div>
                <span>Strand 2</span>
              </div>
            </div>
            """, unsafe_allow_html=True)

            st.markdown("<div style='height:1rem'></div>", unsafe_allow_html=True)
            st.markdown('<div class="section-label">Filter</div>', unsafe_allow_html=True)
            iso_filter = st.multiselect(
                "filter", ["Stable","Partial","Unstable"],
                default=["Stable","Partial","Unstable"],
                label_visibility="collapsed"
            )

            if st.session_state.focus_idx is not None:
                fi = st.session_state.focus_idx
                st.markdown("<div style='height:0.75rem'></div>", unsafe_allow_html=True)
                st.markdown(f'<div class="section-label">Focus · index {fi}</div>', unsafe_allow_html=True)
                if fi < len(tm_list):
                    st.markdown(f"""
                    <div class="focus-panel">
                    Base &nbsp;&nbsp;&nbsp;&nbsp;{dna_in[fi]}<br>
                    Tm &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;{tm_list[fi]["tm"]} °C<br>
                    Status &nbsp;&nbsp;{tm_list[fi]["status"]}
                    </div>""", unsafe_allow_html=True)
                if st.button("Clear focus"):
                    st.session_state.focus_idx = None
                    st.rerun()

        with col_v:
            turns = 3
            z     = np.linspace(-1, 1, N)
            theta = np.linspace(0, turns * 2 * np.pi, N)
            r     = 1.0
            x1 = r * np.cos(theta);          y1 = r * np.sin(theta)
            x2 = r * np.cos(theta + np.pi);  y2 = r * np.sin(theta + np.pi)

            color_map = {
                ("Stable",   False): C_STABLE,
                ("Partial",  False): C_PARTIAL,
                ("Unstable", False): C_UNSTABLE,
                ("Stable",   True):  C_ENZYME,
                ("Partial",  True):  C_ENZYME,
                ("Unstable", True):  C_ENZYME,
            }

            groups = {}
            for i in range(N):
                status = tm_list[i]["status"] if i < len(tm_list) else tm_list[-1]["status"]
                tv     = tm_list[i]["tm"]     if i < len(tm_list) else tm_list[-1]["tm"]
                if status not in iso_filter: continue
                is_enz = i in enz_bits
                key = (status, is_enz)
                if key not in groups: groups[key] = {"xs":[],"ys":[],"zs":[],"txt":[]}
                groups[key]["xs"].extend([x1[i], x2[i], None])
                groups[key]["ys"].extend([y1[i], y2[i], None])
                groups[key]["zs"].extend([z[i],  z[i],  None])
                groups[key]["txt"].extend([
                    f"{dna_in[i]}  ·  idx {i}  ·  {tv} °C",
                    f"{dna_in[i]}  ·  idx {i}  ·  {tv} °C", ""
                ])

            fig = go.Figure()

            fig.add_trace(go.Scatter3d(
                x=x1, y=y1, z=z, mode="lines+markers",
                line=dict(color="rgba(50,100,120,0.55)", width=3),
                marker=dict(size=2.2, color="rgba(55,105,125,0.75)"),
                showlegend=False, hoverinfo="skip"
            ))
            fig.add_trace(go.Scatter3d(
                x=x2, y=y2, z=z, mode="lines+markers",
                line=dict(color="rgba(40,70,100,0.55)", width=3),
                marker=dict(size=2.2, color="rgba(45,75,105,0.75)"),
                showlegend=False, hoverinfo="skip"
            ))

            for (status, is_enz), g in groups.items():
                fig.add_trace(go.Scatter3d(
                    x=g["xs"], y=g["ys"], z=g["zs"],
                    mode="lines",
                    line=dict(color=color_map[(status,is_enz)], width=12 if is_enz else 4),
                    showlegend=False,
                    hovertext=g["txt"], hoverinfo="text"
                ))

            if st.session_state.focus_idx is not None:
                fi_clamped = min(st.session_state.focus_idx, N-1)
                fz = float(z[fi_clamped])
                init_cam = dict(eye=dict(x=3.2, y=0.0, z=0.0),
                                center=dict(x=0, y=0, z=fz), up=dict(x=0, y=0, z=1))
            else:
                init_cam = dict(eye=dict(x=0.0, y=3.2, z=0.0),
                                center=dict(x=0, y=0, z=0), up=dict(x=0, y=0, z=1))

            fig.update_layout(
                scene=dict(
                    xaxis_visible=False, yaxis_visible=False, zaxis_visible=False,
                    aspectmode="manual", aspectratio=dict(x=1, y=1, z=3.5),
                    camera=init_cam, bgcolor="rgba(0,0,0,0)"
                ),
                paper_bgcolor="rgba(0,0,0,0)",
                height=500,
                margin=dict(l=0, r=0, t=0, b=0),
                showlegend=False,
                hoverlabel=dict(
                    bgcolor="rgba(255,255,255,0.97)", bordercolor="#DADCE0",
                    font=dict(family="Inter", size=11, color="#202124")
                )
            )

            st.plotly_chart(fig, use_container_width=True,
                            config=dict(displayModeBar=False),
                            key="helix_chart")

            st.markdown("""
<script>
(function() {
  // Wait for Plotly to be ready, then wire up the rotate/pause buttons
  function initHelixControls() {
    var gd = document.querySelector('[data-testid="stPlotlyChart"] .js-plotly-plot');
    if (!gd) { setTimeout(initHelixControls, 200); return; }

    var rotating = false;
    var rafId    = null;
    var stepRad  = (2 * Math.PI) / (72 * 2); // ~2.5 deg per frame at 60fps

    function getEye() {
      try { return gd._fullLayout.scene.camera.eye; } catch(e) { return {x:0,y:3.2,z:0}; }
    }

    function step() {
      if (!rotating) return;
      var eye = getEye();
      var r   = Math.sqrt(eye.x * eye.x + eye.y * eye.y);
      if (r < 0.01) r = 3.2;
      var angle = Math.atan2(eye.y, eye.x) + stepRad;
      Plotly.relayout(gd, {
        'scene.camera.eye.x': r * Math.cos(angle),
        'scene.camera.eye.y': r * Math.sin(angle)
      });
      rafId = requestAnimationFrame(step);
    }

    // Build the custom button bar
    var wrap = document.querySelector('[data-testid="stPlotlyChart"]');
    if (!wrap) return;

    // Avoid double-injecting on re-renders
    if (wrap.querySelector('.helix-btn-bar')) return;

    var bar = document.createElement('div');
    bar.className = 'helix-btn-bar';
    bar.style.cssText = [
      'position:absolute', 'bottom:12px', 'left:12px', 'z-index:100',
      'display:flex', 'gap:6px'
    ].join(';');

    function makeBtn(label) {
      var b = document.createElement('button');
      b.textContent = label;
      b.style.cssText = [
        'font-family:DM Sans,sans-serif', 'font-size:11px', 'font-weight:500',
        'padding:4px 12px', 'border-radius:16px', 'cursor:pointer',
        'border:1px solid #DADCE0', 'background:rgba(255,255,255,0.95)',
        'color:#5F6368', 'transition:all 0.15s'
      ].join(';');
      b.onmouseenter = function() { b.style.background='#E8F0FE'; b.style.color='#1A73E8'; b.style.borderColor='#1A73E8'; };
      b.onmouseleave = function() { b.style.background='rgba(255,255,255,0.95)'; b.style.color='#5F6368'; b.style.borderColor='#DADCE0'; };
      return b;
    }

    var btnRotate = makeBtn('▶  Rotate');
    var btnPause  = makeBtn('⏸  Pause');

    btnRotate.onclick = function() {
      if (rotating) return;
      rotating = true;
      rafId = requestAnimationFrame(step);
    };
    btnPause.onclick = function() {
      rotating = false;
      if (rafId) { cancelAnimationFrame(rafId); rafId = null; }
    };

    bar.appendChild(btnRotate);
    bar.appendChild(btnPause);

    // The chart wrapper needs relative positioning for absolute children
    wrap.style.position = 'relative';
    wrap.appendChild(bar);
  }

  // Dark mode aware button styling
  var mq = window.matchMedia('(prefers-color-scheme: dark)');
  function patchDarkButtons() {
    var bar = document.querySelector('.helix-btn-bar');
    if (!bar) return;
    var btns = bar.querySelectorAll('button');
    btns.forEach(function(b) {
      if (mq.matches) {
        b.style.background = 'rgba(43,41,48,0.95)';
        b.style.color      = '#CAC4D0';
        b.style.borderColor= '#49454F';
        b.onmouseenter = function() { b.style.background='#3B3944'; b.style.color='#A8C7FA'; b.style.borderColor='#A8C7FA'; };
        b.onmouseleave = function() { b.style.background='rgba(43,41,48,0.95)'; b.style.color='#CAC4D0'; b.style.borderColor='#49454F'; };
      } else {
        b.style.background = 'rgba(255,255,255,0.95)';
        b.style.color      = '#5F6368';
        b.style.borderColor= '#DADCE0';
        b.onmouseenter = function() { b.style.background='#E8F0FE'; b.style.color='#1A73E8'; b.style.borderColor='#1A73E8'; };
        b.onmouseleave = function() { b.style.background='rgba(255,255,255,0.95)'; b.style.color='#5F6368'; b.style.borderColor='#DADCE0'; };
      }
    });
  }
  mq.addEventListener('change', patchDarkButtons);

  initHelixControls();
})();
</script>
""", unsafe_allow_html=True)

    with tabs[1]:
        el, er = st.columns(2, gap="large")

        with el:
            st.markdown('<div class="section-label">Mutation registry</div>', unsafe_allow_html=True)
            if mutations:
                if st.button("Apply all"):
                    new_seq = dna_in
                    for m in mutations: new_seq = apply_mutation(new_seq, m["Idx"], m["To"])
                    st.session_state.dna = new_seq; add_history(new_seq, "Auto-fixed"); st.rerun()
                for m in mutations:
                    ca, cb = st.columns([4, 1])
                    ca.markdown(
                        f'<div class="mut-row">idx <b>{m["Idx"]}</b> &nbsp;'
                        f'<span style="color:var(--unstable)">{m["From"]}</span> → '
                        f'<span style="color:var(--stable)">{m["To"]}</span> &nbsp;'
                        f'<span style="color:var(--text-faint)">{m["Tm"]} °C</span></div>',
                        unsafe_allow_html=True
                    )
                    if cb.button("Apply", key=f"mut_{m['Idx']}"):
                        st.session_state.dna = apply_mutation(dna_in, m["Idx"], m["To"])
                        add_history(st.session_state.dna, f"Fix {m['Idx']}"); st.rerun()
            else:
                st.markdown('<div style="font-family:var(--font);font-size:0.73rem;'
                            'color:#1E4A3A;padding:0.4rem 0">No unstable regions.</div>', unsafe_allow_html=True)

            st.markdown("<div style='height:1rem'></div>", unsafe_allow_html=True)
            st.markdown('<div class="section-label">Primer</div>', unsafe_allow_html=True)
            gc5 = round(gc_fraction(Seq(dna_in[:5])) * 5)
            gc3 = round(gc_fraction(Seq(dna_in[-5:])) * 5)
            st.markdown(
                f'<div class="primer-box">'
                f'<span style="color:var(--text-faint)">5′ → 3′ &nbsp;</span>{dna_in}<br>'
                f'<span style="color:var(--text-faint)">3′ → 5′ &nbsp;</span>{get_rev_comp(dna_in)}'
                f'</div>'
                f'<div style="font-family:var(--font);font-size:0.68rem;'
                f'color:var(--text-faint);margin-top:0.5rem">GC clamp &nbsp; 5′: {gc5}/5 · 3′: {gc3}/5</div>',
                unsafe_allow_html=True
            )

        with er:
            st.markdown('<div class="section-label">Restriction enzyme sites</div>', unsafe_allow_html=True)
            if found_enz:
                for f in found_enz:
                    ca, cb = st.columns([4, 1])
                    active = f["name"] in st.session_state.selected_enzymes
                    badge = f'<span style="color:{C_ENZYME}">◆</span>' if active else '<span style="color:var(--border)">◇</span>'
                    ca.markdown(
                        f'<div class="mut-row">{badge} <b>{f["name"]}</b> &nbsp; idx {f["idx"]}</div>',
                        unsafe_allow_html=True
                    )
                    if cb.button("Focus", key=f"foc_{f['name']}_{f['idx']}"):
                        st.session_state.focus_idx = f["idx"]; st.rerun()
            else:
                st.markdown('<div style="font-family:var(--font);font-size:0.73rem;'
                            'color:var(--text-faint);padding:0.4rem 0">No restriction sites found.</div>', unsafe_allow_html=True)

            st.markdown("<div style='height:1rem'></div>", unsafe_allow_html=True)
            gc5_val = round(gc_fraction(Seq(dna_in[:5])) * 5)
            gc3_val = round(gc_fraction(Seq(dna_in[-5:])) * 5)
            try:
                pdf_bytes = generate_pdf_report(
                    dna_seq=dna_in, n=N, mean_tm=mean_tm, gc_pct=gc_pct,
                    found_enz=found_enz, mutations=mutations,
                    rev_comp=get_rev_comp(dna_in), gc5=gc5_val, gc3=gc3_val
                )
                st.download_button(
                    "Export PDF report", pdf_bytes,
                    file_name="dna_analysis.pdf",
                    mime="application/pdf",
                    use_container_width=True
                )
            except Exception as _pdf_err:
                _msg = (
                    '<div style="font-family: monospace; font-size:0.72rem; '
                    'color:var(--unstable); border:1px solid #2A1010; padding:0.75rem; border-radius:2px;">'
                    'PDF generation requires <b>reportlab</b>.<br>'
                    'Run in your terminal, then restart the app:<br><br>'
                    '<code style="color:#C05050; background:var(--bg-code); padding:3px 6px;">pip install reportlab</code>'
                    '</div>'
                )
                st.markdown(_msg, unsafe_allow_html=True)

    with tabs[2]:
        st.markdown('<div class="section-label">CRISPR-Cas9 Guide RNA Optimizer</div>', unsafe_allow_html=True)

        if N < 23:
            st.markdown(
                _note("Sequence must be at least 23 bp to identify gRNA candidates (20 bp spacer + NGG PAM)."),
                unsafe_allow_html=True
            )
        else:
            cr_left, cr_right = st.columns([3, 2], gap="large")

            with cr_left:
                st.markdown('<div class="section-label">PAM scanner · NGG</div>', unsafe_allow_html=True)

                with st.spinner("Scanning for gRNA candidates…"):
                    candidates = find_grna_candidates(dna_in)

                if not candidates:
                    st.markdown(_note("No NGG PAM sites found in this sequence."), unsafe_allow_html=True)
                else:
                    st.markdown(
                        '<div class="section-label" style="margin-top:0.25rem">Genomic context for off-target scan</div>',
                        unsafe_allow_html=True
                    )
                    uploaded_genome = st.file_uploader(
                        "genome_upload",
                        type=["fa", "fasta", "fna", "txt"],
                        label_visibility="collapsed",
                        help="Upload a FASTA or plain-text genomic sequence. "
                             "Off-target windows will be extracted from real NGG loci in your file."
                    )

                    genome_seq   = None
                    genome_label = None

                    if uploaded_genome is not None:
                        try:
                            genome_seq, genome_label = parse_genome_file(uploaded_genome)
                            if len(genome_seq) < 23:
                                st.markdown(
                                    _note("Uploaded file contains fewer than 23 usable bases after cleaning. "
                                          "Falling back to mock genome."),
                                    unsafe_allow_html=True
                                )
                                genome_seq = None
                            else:
                                st.markdown(
                                    f'<div style="font-family:var(--font);font-size:0.65rem;'
                                    f'color:var(--stable);margin-bottom:0.6rem;padding:0.4rem 0.7rem;'
                                    f'border-left:2px solid var(--stable);background:var(--bg-code)">'
                                    f'✓ {genome_label}</div>',
                                    unsafe_allow_html=True
                                )
                        except Exception as e:
                            st.markdown(
                                _note(f"Could not parse uploaded file: {e}. Using mock genome."),
                                unsafe_allow_html=True
                            )
                            genome_seq = None
                    else:
                        st.markdown(
                            '<div style="font-family:var(--font);font-size:0.63rem;'
                            'color:var(--text-faint);margin-bottom:0.6rem;padding:0.35rem 0.7rem;'
                            'border-left:2px solid var(--border);background:var(--bg-code)">'
                            'No file uploaded — using mock genome for off-target context. '
                            'Upload a FASTA for real-locus scanning.'
                            '</div>',
                            unsafe_allow_html=True
                        )

                    n_decoys = st.sidebar.slider("Off-target windows", 10, 200, 30,
                                                  key="crispr_decoys")
                    top3 = rank_grnas(candidates, n_decoys=n_decoys, genome_seq=genome_seq)

                    st.markdown(
                        f'<div style="font-family:var(--font);font-size:0.68rem;'
                        f'color:var(--text-faint);margin-bottom:0.75rem">'
                        f'{len(candidates)} candidate sites found · showing top 3 by composite score'
                        f'</div>',
                        unsafe_allow_html=True
                    )

                    for rank, g in enumerate(top3, 1):
                        hp_color  = "#8A3030" if g["hairpin"] else "#2A6A4A"
                        hp_icon   = "⚠" if g["hairpin"] else "✓"
                        off_color = "#8A3030" if g["off_target"] > 40 else "#2A6A4A"
                        gc_ok     = 40 <= g["gc_pct"] <= 70
                        gc_color  = "#2A6A4A" if gc_ok else "#8A5020"

                        st.markdown(f"""
                        <div style="background:var(--bg-card);border:1px solid var(--border);border-left:3px solid var(--partial);
                                    border-radius:2px;padding:1rem 1.1rem;margin-bottom:0.85rem">
                          <div style="font-family:'var(--font);font-size:0.74rem;
                                      text-transform:uppercase;letter-spacing:0.18em;color:var(--partial);margin-bottom:0.5rem">
                            Rank #{rank} &nbsp;·&nbsp; Composite {g['composite']:.1f}
                          </div>
                          <div style="font-family:'var(--font-mono);font-size:0.9rem;
                                      color:var(--accent);word-break:break-all;letter-spacing:0.06em;margin-bottom:0.65rem">
                            5'–{g['grna']}–NGG–3'
                          </div>
                          <div style="display:flex;flex-wrap:wrap;gap:0.6rem">
                            <span style="font-family:'var(--font);font-size:0.76rem;
                                         background:var(--bg-code);border:1px solid var(--border);
                                         padding:2px 8px;border-radius:2px;color:#4A7A6A">
                              On-target: {g['on_target']:.1f}
                            </span>
                            <span style="font-family:'var(--font);font-size:0.76rem;
                                         background:var(--bg-code);border:1px solid var(--border);
                                         padding:2px 8px;border-radius:2px;color:{off_color}">
                              Off-risk: {g['off_target']:.1f}
                            </span>
                            <span style="font-family:'var(--font);font-size:0.76rem;
                                         background:var(--bg-code);border:1px solid var(--border);
                                         padding:2px 8px;border-radius:2px;color:{gc_color}">
                              GC: {g['gc_pct']:.0f}%
                            </span>
                            <span style="font-family:'var(--font);font-size:0.76rem;
                                         background:var(--bg-code);border:1px solid var(--border);
                                         padding:2px 8px;border-radius:2px;color:{hp_color}">
                              {hp_icon} Hairpin {'risk' if g['hairpin'] else 'clear'}
                            </span>
                            <span style="font-family:'var(--font);font-size:0.76rem;
                                         background:var(--bg-code);border:1px solid var(--border);
                                         padding:2px 8px;border-radius:2px;color:#3A6070">
                              Strand {g['strand']} · PAM@{g['pam_idx']}
                            </span>
                          </div>
                          <div style="font-family:'var(--font);font-size:0.74rem;
                                      color:var(--text-faint);margin-top:0.45rem">{g['hp_desc']}</div>
                        </div>
                        """, unsafe_allow_html=True)

                    st.markdown('<div class="section-label" style="margin-top:0.5rem">Score comparison</div>',
                                unsafe_allow_html=True)

                    fig_c = go.Figure()
                    labels   = [f"#{i+1}  {g['grna'][:8]}…" for i, g in enumerate(top3)]
                    on_vals  = [g["on_target"]  for g in top3]
                    off_vals = [g["off_target"] for g in top3]
                    comp_vals= [max(g["composite"], 0) for g in top3]

                    for vals, name, color in [
                        (on_vals,   "On-target",   "#2A6A8A"),
                        (off_vals,  "Off-risk",    "#8A3030"),
                        (comp_vals, "Composite",   "#4A7A5A"),
                    ]:
                        fig_c.add_trace(go.Bar(
                            name=name, x=labels, y=vals,
                            marker_color=color, marker_line_width=0,
                            text=[f"{v:.1f}" for v in vals],
                            textposition="outside",
                            textfont=dict(family="Inter", size=9, color="#5F6368")
                        ))

                    fig_c.update_layout(
                        barmode="group", bargap=0.2, bargroupgap=0.05,
                        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                        height=220,
                        margin=dict(l=10, r=10, t=10, b=10),
                        font=dict(family="Inter", size=9, color="#5F6368"),
                        legend=dict(bgcolor="rgba(255,255,255,0.95)", bordercolor="#DADCE0", borderwidth=1,
                                    font=dict(size=9)),
                        yaxis=dict(range=[0, 110], showgrid=True, gridcolor="#E8EAED",
                                   zeroline=False, tickfont=dict(size=8)),
                        xaxis=dict(showgrid=False, tickfont=dict(size=8)),
                    )
                    st.plotly_chart(fig_c, use_container_width=True, config=dict(displayModeBar=False))

                    try:
                        crispr_pdf = generate_crispr_pdf(top3, dna_in)
                        st.download_button(
                            "Export Ready-to-Order PDF", crispr_pdf,
                            file_name="crispr_grna_order.pdf",
                            mime="application/pdf",
                            use_container_width=True
                        )
                    except Exception as _pe:
                        st.markdown(_note(f"PDF export failed: {_pe}"), unsafe_allow_html=True)

            with cr_right:
                st.markdown('<div class="section-label">Scoring guide</div>', unsafe_allow_html=True)
                st.markdown("""
                <div class="legend-panel" style="line-height:1.9">
                  <div class="legend-row">
                    <div class="legend-dot" style="background:var(--partial)"></div>
                    <span><b style="color:var(--partial)">On-target</b>&nbsp; Efficacy estimate. Higher = better cleavage. Rewards 40–70% GC, optimal seed region thermodynamics, no poly-runs.</span>
                  </div>
                  <div class="legend-row" style="margin-top:0.5rem">
                    <div class="legend-dot" style="background:#8A3030"></div>
                    <span><b style="color:var(--unstable)">Off-target risk</b>&nbsp; MIT-style position-weighted score vs. mock-genome ensemble. Lower = safer. Seed region mismatches penalised most.</span>
                  </div>
                  <div class="legend-row" style="margin-top:0.5rem">
                    <div class="legend-dot" style="background:var(--enzyme)"></div>
                    <span><b style="color:var(--enzyme)">Hairpin</b>&nbsp; Self-complementarity check. Stems &ge; 4 bp risk blocking Cas9 loading. Flagged with stem position map.</span>
                  </div>
                  <div class="legend-row" style="margin-top:0.5rem">
                    <div class="legend-dot" style="background:var(--stable)"></div>
                    <span><b style="color:var(--stable)">Composite</b>&nbsp; Weighted ranking score = 0.55 × On − 0.30 × Off-risk − 0.15 × Hairpin. Highest composite = recommended top pick.</span>
                  </div>
                </div>
                <div style="font-family:'var(--font);font-size:0.74rem;color:var(--text-faint);
                            margin-top:1rem;line-height:1.8;padding:0.7rem;border:1px solid var(--border);border-radius:2px">
                  PAM · NGG (SpCas9 default)<br>
                  Spacer · 20 bp upstream of PAM<br>
                  Both strands scanned<br>
                  Off-target model · {n_decoys if len(candidates) > 0 else 30} mock-genome decoys<br>
                  Note: Validate with Cas-OFFinder or CRISPOR before ordering
                </div>
                """, unsafe_allow_html=True)

                if len(candidates) > 0:
                    st.markdown('<div class="section-label" style="margin-top:1.2rem">All candidates</div>',
                                unsafe_allow_html=True)
                    st.markdown(
                        f'<div style="font-family:var(--font);font-size:0.65rem;'
                        f'color:var(--text-faint)">{len(candidates)} NGG sites total — '
                        f'{"only top 3 scored for performance" if len(candidates) > 3 else "all scored"}</div>',
                        unsafe_allow_html=True
                    )
                    for i, c in enumerate(candidates[:8]):
                        st.markdown(
                            f'<div style="font-family:var(--font);font-size:0.65rem;'
                            f'color:var(--text-faint);padding:0.25rem 0;border-bottom:1px solid var(--border)">'
                            f'{i+1}. {c["grna"]} &nbsp; <span style="color:var(--text-faint)">PAM@{c["pam_idx"]} {c["strand"]}</span>'
                            f'</div>',
                            unsafe_allow_html=True
                        )
                    if len(candidates) > 8:
                        st.markdown(
                            _note(f"… and {len(candidates)-8} more sites"),
                            unsafe_allow_html=True
                        )

    with tabs[3]:
        sc1, sc2 = st.columns([3, 1])
        snap_label = sc1.text_input("label", placeholder="Snapshot label...", label_visibility="collapsed")
        if sc2.button("Save snapshot"):
            force_snapshot(dna_in, snap_label or "Snapshot"); st.rerun()

        st.markdown("<hr>", unsafe_allow_html=True)

        cmp_candidates = []
        for idx_h, h in enumerate(st.session_state.history):
            skip = (idx_h == 0 and h["seq"] == dna_in and h["label"] == "Manual")
            if not skip:
                cmp_candidates.append(h)

        if cmp_candidates:
            opts   = [f"{h['time']}  {h['label']}" for h in cmp_candidates]
            choice = st.selectbox("Compare against snapshot", opts,
                                  index=0, label_visibility="collapsed")
            old_seq = cmp_candidates[opts.index(choice)]["seq"]

            diff_html = ""; changed = 0
            maxL = max(len(dna_in), len(old_seq))
            for i in range(maxL):
                curr = dna_in[i]  if i < len(dna_in)  else "-"
                old  = old_seq[i] if i < len(old_seq) else "-"
                if curr == old:
                    diff_html += f'<span style="color:var(--text-faint)">{curr}</span>'
                else:
                    changed += 1
                    diff_html += f'<span class="diff-del">{old}</span><span class="diff-add">{curr}</span>'

            st.markdown("""
            <div style="display:flex;gap:1.5rem;margin-bottom:0.65rem;flex-wrap:wrap">
              <div style="display:flex;align-items:center;gap:0.5rem;font-family:'var(--font);font-size:0.76rem;color:var(--text-muted)">
                <span style="color:var(--text-faint);font-size:0.85rem;letter-spacing:0.05em">ATGC</span>
                <span>Unchanged</span>
              </div>
              <div style="display:flex;align-items:center;gap:0.5rem;font-family:'var(--font);font-size:0.76rem;color:var(--text-muted)">
                <span class="diff-del" style="font-size:0.85rem">X</span>
                <span>Removed / old base</span>
              </div>
              <div style="display:flex;align-items:center;gap:0.5rem;font-family:'var(--font);font-size:0.76rem;color:var(--text-muted)">
                <span class="diff-add" style="font-size:0.85rem">X</span>
                <span>Added / new base</span>
              </div>
              <div style="display:flex;align-items:center;gap:0.5rem;font-family:'var(--font);font-size:0.76rem;color:var(--text-muted)">
                <span style="color:var(--text-faint);font-size:0.85rem">-</span>
                <span>Position absent in shorter sequence</span>
              </div>
            </div>
            """, unsafe_allow_html=True)
            st.markdown(f'<div class="diff-card">{diff_html}</div>', unsafe_allow_html=True)
            identity = (maxL - changed) / maxL * 100 if maxL else 100
            d1, d2, d3 = st.columns(3)
            d1.metric("Changed bases",     changed)
            d2.metric("Sequence identity", f"{identity:.1f}%")
            d3.metric("Length delta",      f"{len(dna_in)-len(old_seq):+d} bp")
        else:
            st.markdown(_note("No saved snapshots differ from the current sequence. "
                              "Save a snapshot before editing, then return here to compare."),
                        unsafe_allow_html=True)

    with tabs[4]:
        hc1, hc2 = st.columns([3, 1])
        q = hc1.text_input("q", placeholder="Search label or sequence...", label_visibility="collapsed").strip().upper()
        if hc2.button("Clear all"):
            st.session_state.history = []; st.rerun()

        if not st.session_state.history:
            st.markdown(_note("No history yet."), unsafe_allow_html=True)
        else:
            filt = [(i, h) for i, h in enumerate(st.session_state.history)
                     if not q or q in h["seq"] or q in h["label"].upper()]
            if not filt:
                st.markdown(_note("No matches."), unsafe_allow_html=True)
            for orig_idx, h in filt:
                ca, cb, cc = st.columns([1, 5, 1])
                ca.markdown(
                    f'<div style="font-family:var(--font);font-size:0.62rem;'
                    f'color:var(--text-faint);line-height:1.8">{h["time"]}<br>'
                    f'<span style="color:var(--accent)">{h["label"]}</span></div>',
                    unsafe_allow_html=True
                )
                cb.markdown(
                    f'<div style="background:var(--bg-code);border:1px solid var(--border);'
                    f'border-radius:var(--radius-sm);padding:0.5rem 0.75rem;'
                    f'font-family:var(--font-mono);font-size:0.78rem;color:var(--text-muted);'
                    f'word-break:break-all;line-height:1.6">{h["seq"]}</div>',
                    unsafe_allow_html=True
                )
                if cc.button("Restore", key=f"res_{orig_idx}"):
                    st.session_state.dna = h["seq"]; st.session_state.focus_idx = None
                    add_history(h["seq"], "Restored"); st.rerun()

elif not dna_in:
    st.markdown(
        '<div style="font-family:var(--font);font-size:0.78rem;'
        'color:var(--text-faint);margin-top:2rem;letter-spacing:0.05em">'
        "Enter a nucleotide sequence in the sidebar to begin."
        "</div>",
        unsafe_allow_html=True
    )
