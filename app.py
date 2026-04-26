import streamlit as st
import re
import pandas as pd
import matplotlib.pyplot as plt

# === Configure page layout (full width) ==================================================
st.set_page_config(
    page_title="NMD Predictor V1.0",
    layout="wide",
)

# === TRANSCRIPT DATABASE ==================================================================
TRANSCRIPTS = {
    "ASXL1_NM_015338.5": (
        "NM_015338.5",
        1669, # c.1669 → NMD cutoff
        4623, # CDS ≈ 3 * 1541
        1541,
    ),
    "TET2_NM_001127208.2": (
        "NM_001127208.2",
        4487, # c.4487 → NMD cutoff
        6006, # CDS ≈ 3 * 2002
        2002,
    ),
}

# === PROTEIN DOMAINS (curated from UniProt) =============================================
PROTEIN_DOMAINS = {
    "ASXL1_NM_015338.5": [
        {"name": "HARE-HTH", "start": 1,   "end": 150,  "color": "#4CAF50"},
        {"name": "DEUBAD",   "start": 300, "end": 550,  "color": "#2196F3"},
        {"name": "PHD",      "start": 1400,"end": 1541, "color": "#FF9800"},
    ],
    "TET2_NM_001127208.2": [
        {"name": "Cys-rich", "start": 1130, "end": 1300, "color": "#9C27B0"},
        {"name": "Tet_JBP",  "start": 1300, "end": 2002, "color": "#FF5722"},
    ],
}

def get_params(gene_tx_key):
    tx, nmd_cut, mrna_len, prot_len = TRANSCRIPTS[gene_tx_key]
    return {
        "gene_tx_key": gene_tx_key,
        "gene": gene_tx_key.split("_")[0],
        "transcript": tx,
        "nmd_cutoff_cdna": nmd_cut,
        "reference_mrna_len": mrna_len, # CDS for % protein
        "protein_length_aa": prot_len,
    }

def get_domains(gene_tx_key):
    return PROTEIN_DOMAINS.get(gene_tx_key, [])

# === HGVS PARSING ========================================================================
def extract_c_pos_from_c_hgvs(hgvs_c):
    """
    Extract the first cDNA position from c. strings:
      - c.1234A>G
      - c.245_247delinsTGA
      - c.234_235del
      - c.234_235insA
      - c.234_235dup
    """
    m = re.search(r"c\.(\d+)", hgvs_c)
    if m:
        return int(m.group(1))
    return None

def parse_p_ptc_position(hgvs_p):
    """
    Given p. text, return the PTC codon number.
    Supports:
      - p.Ser156* / p.Ser156Ter
      - p.Trp343Alafs*39, p.Trp343fs*39, p.Arg78fs*50, etc.
    """
    hgvs_p = hgvs_p.strip()
    # Case 1: simple nonsense (p.Ser156* or p.Ser156Ter)
    m = re.search(r"p\.[A-Z][a-z]{2}(\d+)(?:\*|Ter)", hgvs_p, re.IGNORECASE)
    if m:
        aa_stop = int(m.group(1))
        return aa_stop
    # Case 2: frameshift (p.Trp343Alafs*39, p.Trp343fs*39, p.Arg78fs*50, etc.)
    m = re.search(
        r"p\.[A-Z][a-z]{2}?(\d+)([A-Z][a-z]{2})?fs(?:Ter|\*)?(\d+)",
        hgvs_p,
        re.IGNORECASE
    )
    if m:
        aa_start = int(m.group(1))
        n_aa_new = int(m.group(3)) # the *39 or *60 number
        ptc_codon = aa_start + n_aa_new - 1
        return ptc_codon
    return None

def hgvs_to_ptc_c_pos(hgvs_str):
    """
    From HGVS‑style text (c. + p.), return the PTC cDNA position.
    """
    hgvs_str = hgvs_str.strip()
    if not hgvs_str:
        return None, "Empty string"
    # Extract c. part (flexible: del/delins/ins/dup)
    c_match = re.search(r"c\.\d+.*", hgvs_str, re.IGNORECASE)
    if not c_match:
        return None, "No c. part found"
    c_str = c_match.group()
    c_start = extract_c_pos_from_c_hgvs(c_str)
    if c_start is None:
        return None, "Failed to parse c. position"
    # Extract p. part
    p_match = re.search(r"p\.[^ ]*", hgvs_str, re.IGNORECASE)
    if not p_match:
        return None, "No p. part found"
    p_str = p_match.group()
    ptc_codon = parse_p_ptc_position(p_str)
    if ptc_codon is None:
        return None, "Failed to parse p. frameshift/stop"
    # codon → cDNA (first base of codon)
    ptc_c_pos = 3 * ptc_codon
    return ptc_c_pos, None

# === STREAMLIT LAYOUT ====================================================================
st.markdown("<h1>NMD Predictor V1.0</h1>", unsafe_allow_html=True)

# Create two tabs: Input and Report
tab_input, tab_report = st.tabs(["Input (Tool)", "Report (Structured Output)"])

INPUT_DATA = []

with tab_input:
    st.markdown("""
    **Ownership and use notice:**
    This NMD Predictor is an in‑house tool developed by **Ashley Sunderland**.
    You may use it for internal educational and analytical purposes, but reproduction, redistribution, or commercial use without prior written permission is not permitted.
    This tool is intended for **education and research only** and is **NOT intended for diagnostic purposes**.
    """, unsafe_allow_html=True)

    # Simple format: ASXL1_NM_015338.5 → ASXL1 / NM / 015338.5
    gene_tx_key = st.selectbox(
        "Select gene and transcript:",
        options=list(TRANSCRIPTS.keys()),
        format_func=lambda x: x.replace("_", " / "),
        placeholder="Please select a gene and transcript",
        index=None,                    # Makes dropdown blank on load
        key="gene_tx_key",
    )

    if not gene_tx_key:
        st.info("⚠️ Please select a gene and transcript to continue.")
    else:
        current = get_params(gene_tx_key)
        st.markdown(f"""
        You are currently using:
        - **Gene:** {current['gene']}
        - **Transcript:** {current['transcript']}
        - **NMD cutoff:** cDNA position ≤ {current['nmd_cutoff_cdna']} → NMD predicted
        - **Protein length:** {current['protein_length_aa']} aa
        (CDS ≈ {current['reference_mrna_len']} bp, no premature stop.)
        """)

        # HGVS input: smaller box, no help text, just example line
        st.markdown(
            "<small style='color:#888; font-style:italic; display:block; margin-bottom:4px;'>"
            "Example: c.424_425del p.Arg143Thrfs*110"
            "</small>",
            unsafe_allow_html=True,
        )

        input_text = st.text_area(
            "Paste HGVS variant (c. and p.):",
            height=60, # small height
        )

        if not input_text.strip():
            st.info("Paste a variant in the box above (e.g., `c.424_425del p.Arg143Thrfs*110`).")
        else:
            hgvs_lines = [line.strip() for line in input_text.split("\n") if line.strip()]
            for i, line in enumerate(hgvs_lines):
                st.divider()
                st.markdown(f"**Variant #{i+1}:** `{line}`")
                ptc_c_pos, err = hgvs_to_ptc_c_pos(line)
                if err:
                    st.error(f"Parsing error: {err}")
                    continue
                st.write(f"**PTC codon:** {ptc_c_pos // 3}")
                st.write(f"**PTC cDNA position:** `c.{ptc_c_pos}`")

                # First 100 bp warning
                if ptc_c_pos <= 100:
                    st.warning(
                        "⚠️ **WARNING:** This variant occurs within the first 100 bp. "
                        "Currently this tool does not acknowledge the use of alternative transcripts. "
                        "Use scientific judgement."
                    )

                # NMD vs truncated vs UTR‑extension
                cds_len = current["reference_mrna_len"] # ≈ CDS length
                prot_len = current["protein_length_aa"]
                cds_end = 3 * prot_len

                # Detect frameshift start codon if possible
                frameshift_start_codon = None
                if "fs" in line:
                    m_p = re.search(
                        r"p\.[A-Z][a-z]{2}?(\d+)([A-Z][a-z]{2})?fs(?:Ter|\*)?(\d+)",
                        line,
                        re.IGNORECASE
                    )
                    if m_p:
                        frameshift_start_codon = int(m_p.group(1))

                if ptc_c_pos <= current["nmd_cutoff_cdna"]:
                    # NMD → 100% lost
                    nmd = "YES"
                    nmd_label = "NMD predicted"
                    impact = "Full loss (NMD) – 100% of protein lost"
                    fraction_lost = 1.0
                    perc_lost = 100.0
                    extra = "<span style='color:green; font-weight:bold'>DRIVER</span>"
                elif ptc_c_pos > cds_end:
                    # 3′ UTR PTC → chimera‑like
                    nmd = "NO"
                    nmd_label = "Extended / chimera‑like"
                    impact = "Extended protein (3′ UTR PTC)"
                    if frameshift_start_codon is not None:
                        fraction_corrupted = 1.0 - (frameshift_start_codon - 1) / prot_len
                        fraction_corrupted = max(0.0, fraction_corrupted)
                        perc_corrupted = fraction_corrupted * 100
                        perc_text = f"{perc_corrupted:.1f}%"
                        impact += f" – {perc_text} of canonical protein corrupted by frameshift"
                    else:
                        ptc_codon = ptc_c_pos // 3
                        fraction_corrupted = 1.0 - (ptc_codon - 1) / prot_len
                        fraction_corrupted = max(0.0, fraction_corrupted)
                        perc_corrupted = fraction_corrupted * 100
                        perc_text = f"{perc_corrupted:.1f}%"
                        impact += f" – {perc_text} of canonical protein corrupted by frameshift"
                    extra = (
                        "<span style='color:orange; font-weight:bold'>"
                        "Chimera‑like construct generated. Possible driver – "
                        "please consider % of the canonical protein compromised and downstream loss of function (LOF) variants."
                        "</span>"
                    )
                else:
                    # Truncated protein
                    nmd = "NO"
                    nmd_label = "Truncated protein"
                    impact = "Truncated protein"
                    ratio = ptc_c_pos / cds_len
                    fraction_lost = 1.0 - ratio
                    fraction_lost = max(0.0, fraction_lost)
                    perc_lost = fraction_lost * 100
                    extra = (
                        "<span style='color:orange; font-weight:bold'>"
                        "Possible driver variant – requires assessment of the % of the canonical "
                        "protein compromised and database evidence, including downstream loss of "
                        "function (LOF) variants"
                        "</span>"
                    )

                # === COLLECT DATA FOR REPORT ===================================================
                assessment = "Possible driver"
                mechanism = "Fits mechanism and spectrum of variants"
                if "DRIVER" in extra:
                    assessment = "DRIVER"
                elif "Chimera‑like construct" in extra:
                    mechanism = "Chimera‑like construct"

                INPUT_DATA.append({
                    "Variant": line,
                    "Gene": current["gene"],
                    "Transcript": current["transcript"],
                    "PTC codon": ptc_c_pos // 3,
                    "PTC cDNA": f"c.{ptc_c_pos}",
                    "NMD": nmd_label,
                    "Assessment": assessment,
                    "Mechanism / driver note": mechanism,
                    "Protein impact": impact,
                    "Extra": extra,
                })

                # === SHOW IN INTERACTIVE TAB (immediately visible, no expander) =============
                st.markdown(f"**NMD?:** {nmd} ({nmd_label})")
                st.markdown(f"**Protein impact:** {impact}")
                if "corrupted by frameshift" in impact:
                    # % already in impact text
                    pass
                else:
                    st.markdown(f"**Approx. protein lost:** {fraction_lost:.2f} ({perc_lost:.1f}%)")
                st.markdown(extra, unsafe_allow_html=True)

with tab_report:
    if not INPUT_DATA:
        st.info("Paste and run variants in the 'Input' tab to generate a structured report.")
    else:
        df = pd.DataFrame(INPUT_DATA)
        st.markdown("### Structured report (click a row to see details)")
        # Clean view for the table
        display_df = df[[
            "Variant",
            "Gene",
            "Transcript",
            "NMD",
            "Assessment",
            "Mechanism / driver note",
            "Protein impact",
        ]].copy()
        st.dataframe(display_df, use_container_width=True, hide_index=True)
        st.divider()
        # Simple CSV‑style block you can copy‑paste into Excel / paper
        st.markdown("**CSV‑style summary (for copy‑paste):**")
        st.text(df[["Variant", "Gene", "Transcript", "NMD", "Assessment", "Mechanism / driver note", "Protein impact"]].to_csv(index=False))

# --- === GENE TRACK WITH NMD CUTOFF + DOMAINS (prettier version) === ---
if INPUT_DATA:
    current = get_params(st.session_state.gene_tx_key)
    prot_len = current["protein_length_aa"]
    cds_end = 3 * prot_len
    nmd_cutoff = current["nmd_cutoff_cdna"]
    
    # Use last variant processed to drive the track
    last = INPUT_DATA[-1]
    variant_label = last["Variant"].split()[-1] # e.g. p.Ser21*
   
    # Calculate positions
    codon_ptc = parse_p_ptc_position(last["Variant"].strip())
    ptc_c_pos = 3 * codon_ptc
   
    # Frameshift start
    frameshift_start_codon = None
    if "fs" in last["Variant"]:
        m_p = re.search(r"p\.[A-Z][a-z]{2}?(\d+)", last["Variant"], re.IGNORECASE)
        if m_p: frameshift_start_codon = int(m_p.group(1))
   
    var_origin = 3 * frameshift_start_codon if frameshift_start_codon else ptc_c_pos

    # Domain summary table
    domains = get_domains(st.session_state.gene_tx_key)
    if domains:
        st.markdown("**Protein Domains (AA positions from UniProt):**")
        domain_df = pd.DataFrame(domains)
        st.dataframe(domain_df[["name", "start", "end"]], hide_index=True, use_container_width=True)

    # Plot (prettier)
    fig, ax = plt.subplots(figsize=(14, 3.5), tight_layout=True)
    y = 0
    height = 1.0

    # Light background full-length bar
    ax.barh(y, cds_end, height=height, color="#eeeeee", edgecolor="#444444", alpha=0.7)

    # Draw domains
    for d in domains:
        start_bp = d["start"] * 3
        width_bp = (d["end"] - d["start"] + 1) * 3
        ax.barh(y, width_bp, left=start_bp, height=height, 
                color=d["color"], edgecolor="black", alpha=0.9)
        ax.text(start_bp + width_bp/2, y + 0.55, d["name"], 
                ha="center", va="center", fontsize=9, fontweight="bold", color="white")

    # Variant effect
    ax.barh(y, var_origin, height=height*0.65, color="cornflowerblue", edgecolor="black", label="Intact")
    if var_origin < cds_end:
        ax.barh(y, cds_end - var_origin, left=var_origin, height=height*0.65, 
                color="salmon", edgecolor="black", label="Affected")

    # Formatting
    ax.set_xlim(1, max(cds_end, ptc_c_pos + 300))
    ax.set_ylim(-4, 4)
    ax.set_yticks([])
    ax.set_xlabel("cDNA position along transcript (bp)", fontsize=11)

    ax.text(5, 2.6, "Start", ha="left", va="bottom", fontsize=10)
    ax.text(cds_end, 2.6, "CDS End", ha="right", va="bottom", fontsize=10)

    # NMD Cutoff
    if nmd_cutoff <= cds_end:
        ax.axvline(nmd_cutoff, color="purple", linestyle=":", linewidth=2.8)
        ax.text(nmd_cutoff, -2.2, "NMD cutoff", ha="center", va="top", 
                fontsize=10, color="purple", fontweight="bold")

    # Variant annotations
    ax.annotate(variant_label, xy=(var_origin, 0.4), xytext=(var_origin, 2.9),
                arrowprops=dict(arrowstyle="->", color="black", lw=1.8),
                ha="center", fontsize=11, fontweight="bold")

    ax.annotate("PTC", xy=(ptc_c_pos, -0.6), xytext=(ptc_c_pos, -3.0),
                arrowprops=dict(arrowstyle="->", color="black", lw=1.8), 
                fontsize=11, ha="center", fontweight="bold")

    if ptc_c_pos > cds_end:
        ax.text(cds_end + 30, -0.6, "3′ UTR", fontsize=10, color="#D2691E", ha="left")

    ax.legend(loc="upper right", fontsize=10)
    st.pyplot(fig, use_container_width=True)

# Footer (copyright notice)
st.markdown(
    "<p style='font-size:12px; color:#888; text-align:center; margin-top:20px;'>"
    "© 2026 Ashley Sunderland • NMD Predictor (educational use only, no reproduction without permission)."
    "</p>",
    unsafe_allow_html=True,
)
