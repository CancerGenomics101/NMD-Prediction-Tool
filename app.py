import streamlit as st
import re
import pandas as pd
import matplotlib.pyplot as plt

# === Configure page layout (full width) ==================================================
st.set_page_config(page_title="NMD Predictor V1.0", layout="wide")

# === TRANSCRIPT DATABASE ==================================================================
TRANSCRIPTS = {
    "ASXL1_NM_015338.5": ("NM_015338.5", 1669, 4623, 1541),
    "TET2_NM_001127208.2": ("NM_001127208.2", 4487, 6006, 2002),
}

def get_params(gene_tx_key):
    tx, nmd_cut, mrna_len, prot_len = TRANSCRIPTS[gene_tx_key]
    return {
        "gene_tx_key": gene_tx_key,
        "gene": gene_tx_key.split("_")[0],
        "transcript": tx,
        "nmd_cutoff_cdna": nmd_cut,
        "reference_mrna_len": mrna_len,
        "protein_length_aa": prot_len,
    }

# === HGVS PARSING ========================================================================
def extract_c_pos_from_c_hgvs(hgvs_c):
    m = re.search(r"c\.(\d+)", hgvs_c)
    return int(m.group(1)) if m else None

def parse_p_ptc_position(hgvs_p):
    hgvs_p = hgvs_p.strip()
    m = re.search(r"p\.[A-Z][a-z]{2}(\d+)(?:\*|Ter)", hgvs_p, re.IGNORECASE)
    if m: return int(m.group(1))
    m = re.search(r"p\.[A-Z][a-z]{2}?(\d+)([A-Z][a-z]{2})?fs(?:Ter|\*)?(\d+)", hgvs_p, re.IGNORECASE)
    if m: return int(m.group(1)) + int(m.group(3)) - 1
    return None

def hgvs_to_ptc_c_pos(hgvs_str):
    hgvs_str = hgvs_str.strip()
    c_match = re.search(r"c\.\d+.*", hgvs_str, re.IGNORECASE)
    if not c_match: return None, "No c. part found"
    c_start = extract_c_pos_from_c_hgvs(c_match.group())
    p_match = re.search(r"p\.[^ ]*", hgvs_str, re.IGNORECASE)
    if not p_match: return None, "No p. part found"
    ptc_codon = parse_p_ptc_position(p_match.group())
    return (3 * ptc_codon), None

# === STREAMLIT LAYOUT ====================================================================
st.markdown("<h1>NMD Predictor V1.0</h1>", unsafe_allow_html=True)

selected = st.selectbox(
    "Select gene and transcript:",
    options=list(TRANSCRIPTS.keys()),
    index=None,
    placeholder="Select a transcript...",
    format_func=lambda x: x.replace("_", " / "),
    key="gene_tx_key",
)

if selected is None:
    st.info("Please select a transcript from the dropdown to start.")
else:
    # EVERYTHING BELOW IS INDENTED UNDER 'else'
    current = get_params(selected)
    INPUT_DATA = []
    tab_input, tab_report = st.tabs(["Input (Tool)", "Report (Structured Output)"])

    with tab_input:
        st.markdown("""
        **Ownership and use notice:**  
        This NMD Predictor is an in‑house tool developed by **Ashley Sunderland**.  
        You may use it for internal educational and analytical purposes, but reproduction, redistribution, or commercial use without prior written permission is not permitted.  
        This tool is intended for **education and research only** and is **NOT intended for diagnostic purposes**.
        """, unsafe_allow_html=True)

        st.markdown(f"You are currently using: **Gene:** {current['gene']} | **Transcript:** {current['transcript']}")

        input_text = st.text_area("Paste HGVS variant (c. and p.):", height=60)

        if input_text.strip():
            hgvs_lines = [line.strip() for line in input_text.split("\n") if line.strip()]
            for i, line in enumerate(hgvs_lines):
                st.divider()
                st.markdown(f"**Variant #{i+1}:** `{line}`")
                ptc_c_pos, err = hgvs_to_ptc_c_pos(line)
                if err:
                    st.error(f"Parsing error: {err}")
                    continue
                # ... (rest of your logic goes here, all indented by 4 spaces)
                INPUT_DATA.append({"Variant": line, "PTC cDNA": ptc_c_pos})

    with tab_report:
        if not INPUT_DATA: st.info("Paste variants in 'Input' tab.")
        else: st.dataframe(pd.DataFrame(INPUT_DATA))

    # Gene track logic inside the else block
    if INPUT_DATA:
        # (Add your gene track code here, all indented)
        pass

# === COPYRIGHT FOOTER (Always visible - NO INDENTATION) ===
st.markdown(
    "<p style='font-size:12px; color:#888; text-align:center; margin-top:20px;'>"
    "© 2026 Ashley Sunderland • NMD Predictor (educational use only, no reproduction without permission)."
    "</p>",
    unsafe_allow_html=True,
)
