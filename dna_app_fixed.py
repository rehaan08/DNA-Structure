import streamlit as st
import random
import numpy as np
import plotly.graph_objects as go
import re
from datetime import datetime
import subprocess, sys

# Auto-install reportlab if not present (user may not have it installed)
try:
    import reportlab
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "reportlab",
                           "--quiet", "--disable-pip-version-check"])

# ── PAGE CONFIG ────────────────────────────────────────────────────────────────
st.set_page_config(page_title="DNA Structural Analyzer", layout="wide", initial_sidebar_state="expanded")

# ── DESIGN SYSTEM ──────────────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@300;400;500&family=IBM+Plex+Sans:wght@300;400;500&display=swap');

html, body, .stApp { background-color: #080C10 !important; color: #C8D0D8; font-family: 'IBM Plex Sans', sans-serif; }

[data-testid="stSidebar"] { background-color: #0C1018 !important; border-right: 1px solid #1C2430; }
[data-testid="stSidebar"] * { font-family: 'IBM Plex Mono', monospace !important; font-size: 0.78rem !important; }
[data-testid="stSidebarContent"] { padding: 1.5rem 1rem; }

h1 { font-family: 'IBM Plex Mono', monospace !important; font-size: 1.05rem !important; font-weight: 500 !important;
     letter-spacing: 0.18em !important; color: #6A8A9A !important; text-transform: uppercase; margin-bottom: 0 !important; }

[data-testid="stTabs"] button { font-family: 'IBM Plex Mono', monospace !important; font-size: 0.7rem !important;
    letter-spacing: 0.12em; text-transform: uppercase; color: #3A5060 !important;
    padding: 0.4rem 1.1rem; border: none; background: none !important; }
[data-testid="stTabs"] button[aria-selected="true"] { color: #8AAABA !important; border-bottom: 1px solid #3A6A8A !important; }
[data-testid="stTabs"] [role="tablist"] { border-bottom: 1px solid #131D27; gap: 0; }

[data-testid="stMetric"] { background: #0C1018; border: 1px solid #131D27; border-radius: 2px; padding: 0.75rem 1rem; }
[data-testid="stMetricLabel"] p { font-family: 'IBM Plex Mono', monospace !important; font-size: 0.62rem !important;
    text-transform: uppercase; letter-spacing: 0.14em; color: #3A5060 !important; }
[data-testid="stMetricValue"] { font-family: 'IBM Plex Mono', monospace !important; font-size: 1.25rem !important;
    color: #6A9AAA !important; font-weight: 400 !important; }

.stButton > button { font-family: 'IBM Plex Mono', monospace !important; font-size: 0.7rem !important;
    text-transform: uppercase; letter-spacing: 0.1em; background: transparent !important;
    border: 1px solid #1C2A36 !important; color: #4A6A7A !important; border-radius: 2px !important;
    padding: 0.3rem 0.85rem !important; transition: all 0.12s ease; }
.stButton > button:hover { border-color: #3A6A8A !important; color: #8AAABA !important; background: #0A141E !important; }

[data-testid="stMultiSelect"] { font-size: 0.75rem !important; }
.stMultiSelect [data-baseweb="tag"] { background-color: #0C1828 !important; border: 1px solid #1E3A52 !important;
    font-family: 'IBM Plex Mono', monospace !important; font-size: 0.68rem !important; border-radius: 2px !important; }

.stTextArea textarea, .stTextInput input, .stSelectbox select {
    font-family: 'IBM Plex Mono', monospace !important; font-size: 0.75rem !important;
    background-color: #0C1018 !important; border: 1px solid #1C2A36 !important;
    color: #8AAABA !important; border-radius: 2px !important; }
.stTextArea textarea:focus, .stTextInput input:focus { border-color: #2A4A6A !important; box-shadow: none !important; }

.stSlider [data-testid="stMarkdownContainer"] p { font-family: 'IBM Plex Mono', monospace !important; font-size: 0.7rem !important; }

.section-label { font-family: 'IBM Plex Mono', monospace; font-size: 0.62rem; text-transform: uppercase;
    letter-spacing: 0.18em; color: #2A4050; margin-bottom: 0.65rem; padding-bottom: 0.35rem;
    border-bottom: 1px solid #131D27; }

.diff-card { background: #0C1018; border: 1px solid #131D27; padding: 1.1rem 1.25rem; border-radius: 2px;
    font-family: 'IBM Plex Mono', monospace; font-size: 0.78rem; line-height: 2.3;
    white-space: pre-wrap; word-break: break-all; color: #4A6A7A; }
.diff-del { color: #8A3030; text-decoration: line-through; background: #150C0C; padding: 0 3px; border-radius: 1px; }
.diff-add { color: #308060; font-weight: 500; background: #0A1612; padding: 0 3px; border-radius: 1px; }

.diag-card { background: #0C1018; border: 1px solid #1E1010; padding: 1.1rem 1.25rem; border-radius: 2px;
    font-family: 'IBM Plex Mono', monospace; font-size: 0.78rem; white-space: pre-wrap; word-break: break-all; }

.primer-box { background: #0C1018; border-left: 2px solid #1E3A52; padding: 0.7rem 1rem;
    font-family: 'IBM Plex Mono', monospace; font-size: 0.7rem; color: #4A6A7A;
    line-height: 1.8; word-break: break-all; border-radius: 0 2px 2px 0; }

.legend-panel { background: #0C1018; border: 1px solid #131D27; border-radius: 2px; padding: 1rem 1.1rem; }
.legend-row { display: flex; align-items: center; gap: 0.7rem; margin-bottom: 0.55rem;
    font-family: 'IBM Plex Mono', monospace; font-size: 0.68rem; color: #3A5060; }
.legend-swatch { width: 24px; height: 3px; border-radius: 1px; flex-shrink: 0; }
.legend-dot { width: 7px; height: 7px; border-radius: 50%; flex-shrink: 0; }

.mut-row { font-family: 'IBM Plex Mono', monospace; font-size: 0.75rem; padding: 0.38rem 0;
    border-bottom: 1px solid #0E1820; color: #4A6A7A; }

.focus-panel { background: #0C1018; border: 1px solid #131D27; border-radius: 2px;
    padding: 0.9rem 1.1rem; font-family: 'IBM Plex Mono', monospace; font-size: 0.7rem;
    color: #4A6A7A; line-height: 2.0; }

hr { border: none; border-top: 1px solid #131D27 !important; margin: 1rem 0 !important; }
.stAlert { background: #0C1018 !important; border: 1px solid #131D27 !important; border-radius: 2px !important; }
[data-testid="stDownloadButton"] button { width: 100%; }
#MainMenu, footer, header { visibility: hidden; }
</style>
""", unsafe_allow_html=True)

# ── CONSTANTS ──────────────────────────────────────────────────────────────────
C_STABLE   = "#2A7A6A"
C_PARTIAL  = "#2A4A6A"
C_UNSTABLE = "#7A2A2A"
C_ENZYME   = "#9A7A20"

RE_SITES = {
    "EcoRI":   "GAATTC",   "BamHI":   "GGATCC",
    "HindIII": "AAGCTT",   "NotI":    "GCGGCCGC",
    "SmaI":    "CCCGGG",   "XhoI":    "CTCGAG",
}

# ── SESSION STATE ──────────────────────────────────────────────────────────────
for k, v in [("dna",""), ("history",[]), ("selected_enzymes",[]), ("focus_idx",None)]:
    if k not in st.session_state: st.session_state[k] = v

def add_history(seq, label):
    """Auto-save: skips if seq matches most recent entry (avoids keystroke spam)."""
    if seq and (not st.session_state.history or st.session_state.history[0]["seq"] != seq):
        st.session_state.history.insert(0, {"time": datetime.now().strftime("%H:%M:%S"), "seq": seq, "label": label})

def force_snapshot(seq, label):
    """Always save: used when user explicitly clicks Save Snapshot."""
    if seq:
        st.session_state.history.insert(0, {"time": datetime.now().strftime("%H:%M:%S"), "seq": seq, "label": label})

def get_rev_comp(seq):
    return "".join({"A":"T","C":"G","G":"C","T":"A"}.get(b,b) for b in reversed(seq))

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

    # -- Colour palette: dark text on white for print readability ---------------
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
    # ── Styles ─────────────────────────────────────────────────────────────────
    MONO = "Courier"
    SANS = "Helvetica"

    def safe(s):
        """Sanitise string for ReportLab WinAnsiEncoding (Latin-1 safe)."""
        s = str(s)
        # Handle degree+C together first to avoid double-C artifact
        s = s.replace('°C', ' degC').replace('°', ' deg')
        s = s.replace('→', '->').replace('←', '<-')
        s = s.replace('–', '-').replace('—', ' - ')
        s = s.replace('‘', "'").replace('’', "'")
        s = s.replace('“', '"').replace('”', '"')
        s = s.replace('·', '.').replace('•', '*')
        s = s.replace('×', 'x').replace('≥', '>=').replace('≤', '<=')
        s = s.replace('é', 'e').replace('ö', 'o').replace('ü', 'u')
        s = s.replace('─', '-').replace('│', '|')
        # Final safety net: drop any remaining non-Latin-1
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

    # ── Analytical write-up ────────────────────────────────────────────────────
    grade = "STABLE" if mean_tm > 18 else "UNSTABLE"
    grade_color = "#1A6A5A" if mean_tm > 18 else "#8A1A1A"

    # GC interpretation
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

    # Mutation interpretation
    if len(mutations) == 0:
        mut_interp = ("No unstable windows were identified across the sequence. "
                      "All local melting temperatures meet the minimum threshold, "
                      "indicating uniform structural integrity along the strand.")
    else:
        mut_interp = (f"{len(mutations)} unstable region(s) were identified where the local melting "
                      f"temperature falls below 16 °C. These sites represent potential weak points "
                      f"under thermal or chemical stress and are candidates for targeted mutagenesis "
                      f"— specifically substituting A/T bases with G/C to increase local Tm.")

    # Enzyme interpretation
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

    # Tm interpretation
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

    # ── Build story ────────────────────────────────────────────────────────────
    story = []

    # Header block
    story.append(Paragraph("DNA STRUCTURAL ANALYSIS REPORT", S_title))
    story.append(Paragraph(
        f"Generated {datetime.now().strftime('%Y-%m-%d  %H:%M:%S')}  |  "
        f"Sequence length {n} bp  |  Grade: "
        f'<font color="{grade_color}"><b>{grade}</b></font>',
        S_sub
    ))
    story.append(SP(2)); story.append(HR())

    # Metrics table
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

    # Analytical commentary
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

    # Unstable sites table (if any)
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

    # Enzyme sites table (if any)
    if found_enz:
        story.append(Paragraph("RESTRICTION ENZYME MAP", S_head))
        enz_rows = [["Enzyme", "Index", "Recognition Site"]]
        for fe in found_enz:
            enz_rows.append([fe["name"], str(fe["idx"]), RE_SITES[fe["name"]]])
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

    # Raw sequences
    story.append(Paragraph("RAW SEQUENCE DATA", S_head))
    # Wrap sequence in 60-char lines for readability
    def wrap60(s):
        return "  ".join(s[i:i+60] for i in range(0, len(s), 60))
    story.append(Paragraph("5' -> 3'  (forward)", S_caption))
    story.append(Paragraph(wrap60(dna_seq), S_mono))
    story.append(SP(1))
    story.append(Paragraph("3' -> 5'  (reverse complement)", S_caption))
    story.append(Paragraph(wrap60(rev_comp), S_mono))
    story.append(SP(2)); story.append(HR())

    # Footer note
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
    return (f'<div style="font-family:\'IBM Plex Mono\',monospace;font-size:0.73rem;'
            f'color:#1E3A4A;padding:0.4rem 0">{msg}</div>')

# ── SIDEBAR ────────────────────────────────────────────────────────────────────
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

# ── VALIDATION ─────────────────────────────────────────────────────────────────
invalid_idx = [i for i, c in enumerate(dna_in) if c not in "ATGC"]
is_valid = len(dna_in) >= window_size and not invalid_idx

# ── PAGE HEADER ────────────────────────────────────────────────────────────────
st.markdown("<h1>DNA Structural Analyzer</h1>", unsafe_allow_html=True)
st.markdown("<hr>", unsafe_allow_html=True)

# ── INVALID ────────────────────────────────────────────────────────────────────
if not is_valid and dna_in:
    highlighted = "".join(
        f'<span style="color:#8A3030;font-weight:600;text-decoration:underline">{c}</span>'
        if i in invalid_idx else f'<span style="color:#3A5060">{c}</span>'
        for i, c in enumerate(dna_in)
    )
    st.markdown(f'<div class="diag-card">{highlighted}</div>', unsafe_allow_html=True)
    if st.button("Clean sequence"):
        cleaned = re.sub(r"[^ATGC]", "", dna_in)
        st.session_state.dna = cleaned
        add_history(cleaned, "Cleaned")
        st.rerun()

elif is_valid:
    # ── COMPUTE ────────────────────────────────────────────────────────────────
    N = len(dna_in)
    tm_list, mutations = [], []
    for i in range(N - window_size + 1):
        chunk = dna_in[i:i+window_size]
        tm = 4*(chunk.count("G")+chunk.count("C")) + 2*(chunk.count("A")+chunk.count("T"))
        status = "Stable" if tm > 20 else ("Partial" if tm >= 16 else "Unstable")
        tm_list.append({"idx": i, "tm": tm, "status": status})
        if tm < 16:
            mutations.append({"Idx": i, "From": dna_in[i], "To": "G" if dna_in[i] in "AT" else "C", "Tm": tm})

    found_enz = []
    for name, site in RE_SITES.items():
        for m in re.finditer(f"(?={site})", dna_in):
            found_enz.append({"name": name, "idx": m.start(), "len": len(site)})

    # Enzyme sidebar selector — must be before helix render
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

    # ── METRICS ────────────────────────────────────────────────────────────────
    mean_tm = np.mean([t["tm"] for t in tm_list])
    gc_pct  = (dna_in.count("G") + dna_in.count("C")) / N * 100

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Length",     f"{N} bp")
    c2.metric("Mean Tm",    f"{mean_tm:.1f} degC")
    c3.metric("GC Content", f"{gc_pct:.1f}%")
    c4.metric("RE Sites",   len(found_enz))

    st.markdown("<div style='height:0.5rem'></div>", unsafe_allow_html=True)

    tabs = st.tabs(["Helix · 3D", "Engineering", "Comparison", "History"])

    # ────────────────────────────────────────────────────────────────────────────
    # TAB 1 — HELIX
    # ────────────────────────────────────────────────────────────────────────────
    with tabs[0]:
        col_v, col_l = st.columns([5, 2], gap="medium")

        # ── LEGEND & CONTROLS (right column, renders first logically) ──────────
        with col_l:
            st.markdown("<div style='height:0.6rem'></div>", unsafe_allow_html=True)
            st.markdown('<div class="section-label">Stability key</div>', unsafe_allow_html=True)
            st.markdown(f"""
            <div class="legend-panel">
              <div class="legend-row">
                <div class="legend-swatch" style="background:{C_STABLE}"></div>
                <span>Stable &nbsp;<span style="color:#1A3028">Tm &gt; 20 °C</span></span>
              </div>
              <div class="legend-row">
                <div class="legend-swatch" style="background:{C_PARTIAL}"></div>
                <span>Partial &nbsp;<span style="color:#1A2838">16 – 20 °C</span></span>
              </div>
              <div class="legend-row">
                <div class="legend-swatch" style="background:{C_UNSTABLE}"></div>
                <span>Unstable &nbsp;<span style="color:#381818">Tm &lt; 16 °C</span></span>
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

        # ── HELIX VIEWER (left column) ─────────────────────────────────────────
        with col_v:
            turns = 3
            # z centred on 0 so camera center=(0,0,0) sits exactly mid-helix
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

            # Backbones
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

            # Rungs
            for (status, is_enz), g in groups.items():
                fig.add_trace(go.Scatter3d(
                    x=g["xs"], y=g["ys"], z=g["zs"],
                    mode="lines",
                    line=dict(color=color_map[(status,is_enz)], width=12 if is_enz else 4),
                    showlegend=False,
                    hovertext=g["txt"], hoverinfo="text"
                ))

            # ── ROTATION FRAMES ────────────────────────────────────────────
            # z is centred on 0, so all camera centers use z=0 (true midpoint).
            # eye.z=0 keeps the viewpoint level with the helix centre.
            n_frames = 72
            cam_r    = 3.2
            cam_elev = 0.0   # eye level with z=0 centre — no vertical offset

            frames = [
                go.Frame(
                    layout=dict(scene=dict(camera=dict(
                        eye=dict(x=cam_r*np.cos(2*np.pi*k/n_frames),
                                 y=cam_r*np.sin(2*np.pi*k/n_frames),
                                 z=cam_elev),
                        center=dict(x=0, y=0, z=0),
                        up=dict(x=0, y=0, z=1)
                    ))),
                    name=str(k)
                )
                for k in range(n_frames)
            ]
            fig.frames = frames

            # Focus or default camera
            if st.session_state.focus_idx is not None:
                fi_clamped = min(st.session_state.focus_idx, N-1)
                fz = float(z[fi_clamped])
                init_cam = dict(eye=dict(x=cam_r, y=0.0, z=0.0),
                                center=dict(x=0, y=0, z=fz), up=dict(x=0, y=0, z=1))
            else:
                init_cam = dict(eye=dict(x=0.0, y=cam_r, z=cam_elev),
                                center=dict(x=0, y=0, z=0), up=dict(x=0, y=0, z=1))

            fig.update_layout(
                scene=dict(
                    xaxis_visible=False, yaxis_visible=False, zaxis_visible=False,
                    aspectmode="manual", aspectratio=dict(x=1, y=1, z=3.5),
                    camera=init_cam, bgcolor="#080C10"
                ),
                paper_bgcolor="#080C10",
                height=500,
                margin=dict(l=0, r=0, t=0, b=0),
                showlegend=False,
                updatemenus=[dict(
                    type="buttons", showactive=False,
                    x=0.01, y=0.01, xanchor="left", yanchor="bottom",
                    bgcolor="#0C1018", bordercolor="#1C2A36",
                    font=dict(family="IBM Plex Mono", size=10, color="#4A6A7A"),
                    buttons=[
                        dict(label="▶  Rotate", method="animate",
                             args=[None, dict(frame=dict(duration=40, redraw=True),
                                              fromcurrent=True, transition=dict(duration=0),
                                              mode="immediate")]),
                        dict(label="⏸  Pause", method="animate",
                             args=[[None], dict(frame=dict(duration=0, redraw=False),
                                                mode="immediate")]),
                    ]
                )],
                hoverlabel=dict(
                    bgcolor="#0C1018", bordercolor="#1E3A52",
                    font=dict(family="IBM Plex Mono", size=11, color="#7AAABA")
                )
            )

            st.plotly_chart(fig, use_container_width=True, config=dict(displayModeBar=False))

    # ────────────────────────────────────────────────────────────────────────────
    # TAB 2 — ENGINEERING
    # ────────────────────────────────────────────────────────────────────────────
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
                        f'<span style="color:#6A2828">{m["From"]}</span> → '
                        f'<span style="color:#286050">{m["To"]}</span> &nbsp;'
                        f'<span style="color:#1E3A3A">{m["Tm"]} °C</span></div>',
                        unsafe_allow_html=True
                    )
                    if cb.button("Apply", key=f"mut_{m['Idx']}"):
                        st.session_state.dna = apply_mutation(dna_in, m["Idx"], m["To"])
                        add_history(st.session_state.dna, f"Fix {m['Idx']}"); st.rerun()
            else:
                st.markdown('<div style="font-family:\'IBM Plex Mono\',monospace;font-size:0.73rem;'
                            'color:#1E4A3A;padding:0.4rem 0">No unstable regions.</div>', unsafe_allow_html=True)

            st.markdown("<div style='height:1rem'></div>", unsafe_allow_html=True)
            st.markdown('<div class="section-label">Primer</div>', unsafe_allow_html=True)
            gc5 = dna_in[:5].count("G")+dna_in[:5].count("C")
            gc3 = dna_in[-5:].count("G")+dna_in[-5:].count("C")
            st.markdown(
                f'<div class="primer-box">'
                f'<span style="color:#1E3A52">5′ → 3′ &nbsp;</span>{dna_in}<br>'
                f'<span style="color:#1E3A52">3′ → 5′ &nbsp;</span>{get_rev_comp(dna_in)}'
                f'</div>'
                f'<div style="font-family:\'IBM Plex Mono\',monospace;font-size:0.68rem;'
                f'color:#2A4050;margin-top:0.5rem">GC clamp &nbsp; 5′: {gc5}/5 · 3′: {gc3}/5</div>',
                unsafe_allow_html=True
            )

        with er:
            st.markdown('<div class="section-label">Restriction enzyme sites</div>', unsafe_allow_html=True)
            if found_enz:
                for f in found_enz:
                    ca, cb = st.columns([4, 1])
                    active = f["name"] in st.session_state.selected_enzymes
                    badge = f'<span style="color:{C_ENZYME}">◆</span>' if active else '<span style="color:#1E2E3A">◇</span>'
                    ca.markdown(
                        f'<div class="mut-row">{badge} <b>{f["name"]}</b> &nbsp; idx {f["idx"]}</div>',
                        unsafe_allow_html=True
                    )
                    if cb.button("Focus", key=f"foc_{f['name']}_{f['idx']}"):
                        st.session_state.focus_idx = f["idx"]; st.rerun()
            else:
                st.markdown('<div style="font-family:\'IBM Plex Mono\',monospace;font-size:0.73rem;'
                            'color:#1E3A4A;padding:0.4rem 0">No restriction sites found.</div>', unsafe_allow_html=True)

            st.markdown("<div style='height:1rem'></div>", unsafe_allow_html=True)
            gc5_val = dna_in[:5].count("G")+dna_in[:5].count("C")
            gc3_val = dna_in[-5:].count("G")+dna_in[-5:].count("C")
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
                    'color:#8A3030; border:1px solid #2A1010; padding:0.75rem; border-radius:2px;">'
                    'PDF generation requires <b>reportlab</b>.<br>'
                    'Run in your terminal, then restart the app:<br><br>'
                    '<code style="color:#C05050; background:#1A0808; padding:3px 6px;">pip install reportlab</code>'
                    '</div>'
                )
                st.markdown(_msg, unsafe_allow_html=True)

    # ────────────────────────────────────────────────────────────────────────────
    # TAB 3 — COMPARISON
    # ────────────────────────────────────────────────────────────────────────────
    with tabs[2]:
        sc1, sc2 = st.columns([3, 1])
        snap_label = sc1.text_input("label", placeholder="Snapshot label...", label_visibility="collapsed")
        if sc2.button("Save snapshot"):
            force_snapshot(dna_in, snap_label or "Snapshot"); st.rerun()

        st.markdown("<hr>", unsafe_allow_html=True)

        # Only offer snapshots that differ from the current sequence,
        # so the user never accidentally compares a sequence to itself.
        # Exclude history[0] only if it's the same seq AND label is "Manual"
        # (auto-save). Force-saved snapshots with same seq ARE useful for labelling.
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
                    diff_html += f'<span style="color:#2A4050">{curr}</span>'
                else:
                    changed += 1
                    diff_html += f'<span class="diff-del">{old}</span><span class="diff-add">{curr}</span>'

            # ── Diff legend ────────────────────────────────────────────────
            st.markdown("""
            <div style="display:flex;gap:1.5rem;margin-bottom:0.65rem;flex-wrap:wrap">
              <div style="display:flex;align-items:center;gap:0.5rem;font-family:'IBM Plex Mono',monospace;font-size:0.65rem;color:#3A5060">
                <span style="color:#2A4050;font-size:0.85rem;letter-spacing:0.05em">ATGC</span>
                <span>Unchanged</span>
              </div>
              <div style="display:flex;align-items:center;gap:0.5rem;font-family:'IBM Plex Mono',monospace;font-size:0.65rem;color:#3A5060">
                <span class="diff-del" style="font-size:0.85rem">X</span>
                <span>Removed / old base</span>
              </div>
              <div style="display:flex;align-items:center;gap:0.5rem;font-family:'IBM Plex Mono',monospace;font-size:0.65rem;color:#3A5060">
                <span class="diff-add" style="font-size:0.85rem">X</span>
                <span>Added / new base</span>
              </div>
              <div style="display:flex;align-items:center;gap:0.5rem;font-family:'IBM Plex Mono',monospace;font-size:0.65rem;color:#3A5060">
                <span style="color:#2A4050;font-size:0.85rem">-</span>
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

    # ────────────────────────────────────────────────────────────────────────────
    # TAB 4 — HISTORY
    # ────────────────────────────────────────────────────────────────────────────
    with tabs[3]:
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
                    f'<div style="font-family:\'IBM Plex Mono\',monospace;font-size:0.62rem;'
                    f'color:#2A4050;line-height:1.8">{h["time"]}<br>'
                    f'<span style="color:#2A5060">{h["label"]}</span></div>',
                    unsafe_allow_html=True
                )
                cb.code(h["seq"], language=None)
                # Use the stable original-list index as the widget key
                if cc.button("Restore", key=f"res_{orig_idx}"):
                    st.session_state.dna = h["seq"]; st.session_state.focus_idx = None
                    add_history(h["seq"], "Restored"); st.rerun()

elif not dna_in:
    st.markdown(
        '<div style="font-family:\'IBM Plex Mono\',monospace;font-size:0.78rem;'
        'color:#1E3A4A;margin-top:2rem;letter-spacing:0.05em">'
        "Enter a nucleotide sequence in the sidebar to begin."
        "</div>",
        unsafe_allow_html=True
    )