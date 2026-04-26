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
        1669,
        4623,
        1541,
    ),
    "TET2_NM_001127208.2": (
        "NM_001127208.2",
        4487,
        6006,
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

# === EXON BOUNDARIES (AA) - Accurate from NCBI RefSeq for exact transcripts =============
EXONS = {
    "ASXL1_NM_015338.5": [
        (1,   19,  "ex1"),  (20,  47,  "ex2"),  (48,  48,  "ex3"),
        (49,  84,  "ex4"),  (85, 124,  "ex5"),  (125, 157, "ex6"),
        (158, 188, "ex7"),  (189, 239, "ex8"),  (240, 294, "ex9"),
        (295, 327, "ex10"), (328, 362, "ex11"), (363, 573, "ex12"),
        (574,1541, "ex13"),
    ],
    "TET2_NM_001127208.2": [
        (1,   60,  "ex1"),  (61, 117, "ex2"),  (118, 173, "ex3"),
        (174, 240, "ex4"),  (241, 317, "ex5"),  (318, 393, "ex6"),
        (394, 473, "ex7"),  (474, 550, "ex8"),  (551, 607, "ex9"),
        (608,1496, "ex10"), (1497,2002,"ex11"),
    ],
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

def get_domains(gene_tx_key):
    return PROTEIN_DOMAINS.get(gene_tx_key, [])

def get_exons(gene_tx_key):
    return EXONS.get(gene_tx_key, [])

# === HGVS PARSING ========================================================================
def extract_c_pos_from_c_hgvs(hgvs_c):
    m = re.search(r"c\.(\d+)", hgvs_c)
    if m:
        return int(m.group(1))
    return None

def parse_p_ptc_position(hgvs_p):
    hgvs_p = hgvs_p.strip()
    m = re.search(r"p\.[A-Z][a-z]{2}(\d+)(?:\*|Ter)", hgvs_p, re.IGNORECASE)
    if m:
        return int(m.group(1))
    m = re.search(
        r"p\.[A-Z][a-z]{2}?(\d+)([A-Z][a-z]{2})?fs(?:Ter|\*)?(\d+)",
        hgvs_p,
        re.IGNORECASE
    )
    if m:
        aa_start = int(m.group(1))
        n_aa_new = int(m.group(3))
        return aa_start + n_aa_new - 1
    return None

def hgvs_to_ptc_c_pos(hgvs_str):
    hgvs_str = hgvs_str.strip()
    if not hgvs_str:
        return None, "Empty string"
    c_match = re.search(r"c\.\d+.*", hgvs_str, re.IGNORECASE)
    if not c_match:
        return None, "No c. part found"
    c_str = c_match.group()
    c_start = extract_c_pos_from_c_hgvs(c_str)
    if c_start is None:
        return None, "Failed to parse c. position"
    p_match = re.search(r"p\.[^ ]*", hgvs_str, re.IGNORECASE)
    if not p_match:
        return None, "No p. part found"
    p_str = p_match.group()
    ptc_codon = parse_p_ptc_position(p_str)
    if ptc_codon is None:
        return None, "Failed to parse p. frameshift/stop"
    ptc_c_pos = 3 * ptc_codon
    return ptc_c_pos, None

# === STREAMLIT LAYOUT ====================================================================
st.markdown("<h1>NMD Predictor V1.0</h1>", unsafe_allow_html=True)

tab_input, tab_report = st.tabs(["Input (Tool)", "Report (Structured Output)"])

INPUT_DATA = []

with tab_input:
    st.markdown("""
    **Ownership and use notice:**
    This NMD Predictor is an in‑house tool developed by **Ashley Sunderland**.
    You may use it for internal educational and analytical purposes, but reproduction, redistribution, or commercial use without prior written permission is not permitted.
    This tool is intended for **education and research only** and is **NOT intended for diagnostic purposes**.
    """, unsafe_allow_html=True)

    gene_tx_key = st.selectbox(
        "Select gene and transcript:",
        options=list(TRANSCRIPTS.keys()),
        format_func=lambda x: x.replace("_", " / "),
        placeholder="Please select a gene and transcript",
        index=None,
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

        st.markdown(
            "<small style='color:#888; font-style:italic; display:block; margin-bottom:4px;'>"
            "Example: c.424_425del p.Arg143Thrfs*110"
            "</small>",
            unsafe_allow_html=True,
        )

        input_text = st.text_area(
            "Paste HGVS variant (c. and p.):",
            height=60,
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

                if ptc_c_pos <= 100:
                    st.warning(
                        "⚠️ **WARNING:** This variant occurs within the first 100 bp. "
                        "Currently this tool does not acknowledge the use of alternative transcripts. "
                        "Use scientific judgement."
                    )

                cds_len = current["reference_mrna_len"]
                prot_len = current["protein_length_aa"]
                cds_end = 3 * prot_len

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
                    nmd = "YES"
                    nmd_label = "NMD predicted"
                    impact = "Full loss (NMD) – 100% of protein lost"
                    fraction_lost = 1.0
                    perc_lost = 100.0
                    extra = "<span style='color:green; font-weight:bold'>DRIVER</span>"
                elif ptc_c_pos > cds_end:
                    nmd = "NO"
                    nmd_label = "Extended / chimera‑like"
                    impact = "Extended protein (3′ UTR PTC)"
                    if frameshift_start_codon is not None:
                        fraction_corrupted = max(0.0, 1.0 - (frameshift_start_codon - 1) / prot_len)
                        perc_text = f"{fraction_corrupted*100:.1f}%"
                        impact += f" – {perc_text} of canonical protein corrupted by frameshift"
                    else:
                        ptc_codon = ptc_c_pos // 3
                        fraction_corrupted = max(0.0, 1.0 - (ptc_codon - 1) / prot_len)
                        perc_text = f"{fraction_corrupted*100:.1f}%"
                        impact += f" – {perc_text} of canonical protein corrupted by frameshift"
                    extra = (
                        "<span style='color:orange; font-weight:bold'>"
                        "Chimera‑like construct generated. Possible driver – "
                        "please consider % of the canonical protein compromised and downstream loss of function (LOF) variants."
                        "</span>"
                    )
                else:
                    nmd = "NO"
                    nmd_label = "Truncated protein"
                    impact = "Truncated protein"
                    ratio = ptc_c_pos / cds_len
                    fraction_lost = max(0.0, 1.0 - ratio)
                    perc_lost = fraction_lost * 100
                    extra = (
                        "<span style='color:orange; font-weight:bold'>"
                        "Possible driver variant – requires assessment of the % of the canonical "
                        "protein compromised and database evidence, including downstream loss of "
                        "function (LOF) variants"
                        "</span>"
                    )

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

                st.markdown(f"**NMD?:** {nmd} ({nmd_label})")
                st.markdown(f"**Protein impact:** {impact}")
                if "corrupted by frameshift" not in impact:
                    st.markdown(f"**Approx. protein lost:** {fraction_lost:.2f} ({perc_lost:.1f}%)")
                st.markdown(extra, unsafe_allow_html=True)

with tab_report:
    if not INPUT_DATA:
        st.info("Paste and run variants in the 'Input' tab to generate a structured report.")
    else:
        df = pd.DataFrame(INPUT_DATA)
        st.markdown("### Structured report (click a row to see details)")
        display_df = df[[
            "Variant", "Gene", "Transcript", "NMD", "Assessment",
            "Mechanism / driver note", "Protein impact"
        ]].copy()
        st.dataframe(display_df, use_container_width=True, hide_index=True)
        st.divider()
        st.markdown("**CSV‑style summary (for copy‑paste):**")
        st.text(df[["Variant", "Gene", "Transcript", "NMD", "Assessment", 
                    "Mechanism / driver note", "Protein impact"]].to_csv(index=False))

# --- === GENE TRACK - AMINO ACID SCALE WITH EXONS === ---
if INPUT_DATA:
    current = get_params(st.session_state.gene_tx_key)
    prot_len = current["protein_length_aa"]
    nmd_cutoff_aa = current["nmd_cutoff_cdna"] // 3
    
    last = INPUT_DATA[-1]
    variant_label = last["Variant"].split()[-1]
   
    codon_ptc = parse_p_ptc_position(last["Variant"].strip())
    ptc_aa = codon_ptc if codon_ptc else 0
   
    frameshift_start_codon = None
    if "fs" in last["Variant"]:
        m_p = re.search(r"p\.[A-Z][a-z]{2}?(\d+)", last["Variant"], re.IGNORECASE)
        if m_p: 
            frameshift_start_codon = int(m_p.group(1))
   
    var_origin_aa = frameshift_start_codon if frameshift_start_codon else ptc_aa

    domains = get_domains(st.session_state.gene_tx_key)
    if domains:
        st.markdown("**Protein Domains (AA positions from UniProt):**")
        st.dataframe(pd.DataFrame(domains)[["name", "start", "end"]], hide_index=True, use_container_width=True)

    st.caption("⚠️ **Note:** Exon boundaries are based on NCBI RefSeq for the exact transcripts. "
               "They are approximate AA mappings for visualisation only.")

    # Plot
    fig, ax = plt.subplots(figsize=(14, 5.5), tight_layout=True)
    y = 0
    height = 1.0

    ax.barh(y, prot_len, height=height, color="#f5f5f5", edgecolor="#666", alpha=0.8)

    # Exons
    exons = get_exons(st.session_state.gene_tx_key)
    for start, end, label in exons:
        ax.axvline(start, color="#444444", linestyle="--", linewidth=1.0, alpha=0.6)
        ax.text(start + (end - start)/2, y + 0.5, label, 
                ha="center", va="center", fontsize=9, fontweight="bold", 
                color="black", bbox=dict(facecolor="white", alpha=0.85, pad=1))

    # Domains
    for d in domains:
        width = d["end"] - d["start"] + 1
        ax.barh(y, width, left=d["start"], height=height, 
                color=d["color"], edgecolor="black", alpha=0.9)
        ax.text(d["start"] + width/2, y + 1.25, d["name"], 
                ha="center", va="bottom", fontsize=10, fontweight="bold", 
                color="white", bbox=dict(facecolor='black', alpha=0.75, pad=2))

    # Variant effect
    ax.barh(y, var_origin_aa, height=height*0.75, color="cornflowerblue", edgecolor="black", label="Intact")
    if var_origin_aa < prot_len:
        ax.barh(y, prot_len - var_origin_aa, left=var_origin_aa, height=height*0.75, 
                color="salmon", edgecolor="black", label="Affected")

    ax.set_xlim(1, max(prot_len, ptc_aa + 100))
    ax.set_ylim(-5.5, 6.5)
    ax.set_yticks([])
    ax.set_xlabel("Amino Acid Position", fontsize=12, fontweight="bold")

    ax.text(5, 3.4, "Start", ha="left", va="bottom", fontsize=11)
    ax.text(prot_len, 3.4, "End", ha="right", va="bottom", fontsize=11)

    if nmd_cutoff_aa <= prot_len:
        ax.axvline(nmd_cutoff_aa, color="purple", linestyle=":", linewidth=3)
        ax.text(nmd_cutoff_aa, -3.8, f"NMD cutoff (AA {nmd_cutoff_aa})", 
                ha="center", va="top", fontsize=10, color="purple", fontweight="bold")

    ax.annotate(variant_label, xy=(var_origin_aa, 0.7), xytext=(var_origin_aa, 4.2),
                arrowprops=dict(arrowstyle="->", color="black", lw=2), 
                ha="center", fontsize=11, fontweight="bold")

    ax.annotate("PTC", xy=(ptc_aa, -1.2), xytext=(ptc_aa, -4.8),
                arrowprops=dict(arrowstyle="->", color="black", lw=2), 
                fontsize=11, ha="center", fontweight="bold")

    ax.legend(loc="upper right", fontsize=10)
    st.pyplot(fig, use_container_width=True)

# Footer
st.markdown(
    "<p style='font-size:12px; color:#888; text-align:center; margin-top:20px;'>"
    "© 2026 Ashley Sunderland • NMD Predictor (educational use only, no reproduction without permission)."
    "</p>",
    unsafe_allow_html=True,
)
