import streamlit as st
import re
import pandas as pd
import matplotlib.pyplot as plt
import random
import sys
from pathlib import Path

# === Configure page layout (full width) ==================================================
st.set_page_config(
    page_title="NMD Predictor v1.0",
    layout="wide",
)

# Add current directory to path
sys.path.append(str(Path(__file__).parent))

# Import data from data.py
from data import TRANSCRIPTS, PROTEIN_DOMAINS, EXONS, EDUCATIONAL_FACTS

# === Helper Functions ================================================================
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

# === SPLICE SITE DETECTION ===========================================================
def is_canonical_splice_site_variant(hgvs_c):
    if not hgvs_c:
        return False
    splice_patterns = [
        r"[cC]\.\d+[+-][12]",
        r"[cC]\.\d+_\d+[+-][12]",
        r"[cC]\.\d+[+-]\d+_\d+",
        r"[+-][12][delinsdup]",
    ]
    hgvs_lower = hgvs_c.lower()
    for pattern in splice_patterns:
        if re.search(pattern, hgvs_lower):
            return True
    return False

# === S-VIG O2 STRENGTH SUGGESTION (Improved with full vs partial domain logic) ========
def get_svig_o2_suggestion(ptc_c_pos: int, prot_len: int, nmd_cutoff: int, 
                           frameshift_start_codon: int = None, gene_tx_key: str = None):
    """
    Returns (svig_code, explanation, caveat) following SVIG-UK O2 guidelines.
    """
    corruption_aa = frameshift_start_codon if frameshift_start_codon is not None else (ptc_c_pos // 3)
    percent_lost = max(0.0, 1.0 - corruption_aa / prot_len) * 100

    domains = get_domains(gene_tx_key)
    
    full_domains_lost = []
    partial_domains_lost = []

    for d in domains:
        domain_start = d["start"]
        domain_end = d["end"]
        
        if corruption_aa <= domain_start:
            full_domains_lost.append(d["name"])
        elif corruption_aa < domain_end and corruption_aa > domain_start:
            partial_domains_lost.append(d["name"])

    if ptc_c_pos <= nmd_cutoff:
        # NMD predicted
        code = "O2_VSTR"
        expl = "NMD predicted → Very Strong"
        caveat = ""
    else:
        # NMD evaded
        if full_domains_lost or percent_lost >= 15:
            code = "O2_STR"
            expl = f"NMD evaded but significant impact ({percent_lost:.1f}% lost"
            if full_domains_lost:
                expl += f", full loss of: {', '.join(full_domains_lost)})"
            else:
                expl += ")"
            caveat = ""
        elif partial_domains_lost or percent_lost >= 10:
            code = "O2_STR"
            expl = f"NMD evaded, {percent_lost:.1f}% protein lost"
            if partial_domains_lost:
                caveat = (f"⚠️ Partial loss of domain(s): {', '.join(partial_domains_lost)}. "
                          "Please manually review whether the remaining portion of the domain is still functional. "
                          "Consider downgrading to O2_Mod if the partial domain is not critical.")
            else:
                caveat = ""
        elif percent_lost >= 5:
            code = "O2_Mod"
            expl = f"NMD evaded, moderate C-terminal truncation ({percent_lost:.1f}%)"
            caveat = ""
        else:
            code = "O2_Supp"
            expl = f"NMD evaded, small C-terminal truncation ({percent_lost:.1f}%)"
            caveat = ""

    return code, expl, caveat

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
st.markdown("<h1>NMD Predictor v1.0</h1>", unsafe_allow_html=True)

tab_input, tab_report = st.tabs(["Input (Tool)", "Report (Structured Output)"])

INPUT_DATA = []

with tab_input:
    st.markdown("""
    **Ownership and use notice:**  
    This NMD Predictor is an in‑house tool developed by **Ashley Sunderland**.  
    You may use it for internal educational and analytical purposes, but reproduction, redistribution, or commercial use without prior written permission is not permitted.  
    This tool is intended for **education and research only** and is **NOT intended for diagnostic purposes**.
    """, unsafe_allow_html=True)

    sorted_genes = sorted(TRANSCRIPTS.keys())
    gene_tx_key = st.selectbox(
        "Select gene and transcript:",
        options=sorted_genes,
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
        """)

        st.markdown(
            "<small style='color:#888; font-style:italic; display:block; margin-bottom:4px;'>"
            "Example: c.424_425del p.Arg143Thrfs*110"
            "</small>",
            unsafe_allow_html=True,
        )

        input_text = st.text_area(
            "Paste HGVS variant (c. and p.):",
            height=100,
        )

        if input_text.strip():
            hgvs_lines = [line.strip() for line in input_text.split("\n") if line.strip()]

            for i, line in enumerate(hgvs_lines):
                st.divider()
                st.markdown(f"**Variant #{i+1}:** `{line}`")

                # === SPLICE SITE CHECK ===
                c_part = re.search(r"c\.[^ ]*", line, re.IGNORECASE)
                c_str = c_part.group() if c_part else ""

                if is_canonical_splice_site_variant(c_str):
                    st.markdown("""
                    <div style="background-color:#ffe6e6; padding:20px; border-radius:10px; 
                                border:3px solid #ff4444; margin:10px 0;">
                        <h4 style="color:#cc0000; margin:0 0 8px 0;">⚠️ Warning, further analysis not possible.</h4>
                        <strong>This variant impacts a canonical splice site.</strong><br>
                        Please analyse this variant as a splice consequence.
                    </div>
                    """, unsafe_allow_html=True)
                    continue

                # ==================== NORMAL PROCESSING ====================
                ptc_c_pos, err = hgvs_to_ptc_c_pos(line)
                if err:
                    st.error(f"Parsing error: {err}")
                    continue

                st.write(f"**PTC codon:** {ptc_c_pos // 3}")
                st.write(f"**PTC cDNA position:** `c.{ptc_c_pos}`")

                if ptc_c_pos <= 100:
                    st.warning(
                        "⚠️ **WARNING:** This variant generates a PTC within the first 100 bp. "
                        "Currently, this tool does not acknowledge the use of alternative transcripts. "
                        "Use scientific judgement."
                    )

                prot_len = current["protein_length_aa"]
                cds_end = 3 * prot_len

                frameshift_start_codon = None
                if "fs" in line.lower():
                    m_p = re.search(
                        r"p\.[A-Z][a-z]{2}?(\d+)([A-Z][a-z]{2})?fs(?:Ter|\*)?(\d+)",
                        line,
                        re.IGNORECASE
                    )
                    if m_p:
                        frameshift_start_codon = int(m_p.group(1))

                # ==================== CORE LOGIC ====================
                if ptc_c_pos <= current["nmd_cutoff_cdna"]:
                    nmd_label = "NMD predicted"
                    impact = "Full loss (NMD) – 100% of protein lost"
                    fraction_lost = 1.0
                    perc_lost = 100.0
                    extra = "<span style='color:green; font-weight:bold'>DRIVER</span>"
                    card_color = "🟢"
                elif ptc_c_pos > cds_end:
                    nmd_label = "Extended / chimera‑like"
                    impact = "Extended protein (3′ UTR PTC)"
                    fraction_lost = 0.0
                    perc_lost = 0.0
                    extra = "<span style='color:orange; font-weight:bold'>Chimera‑like construct generated. Possible driver.</span>"
                    card_color = "🧬"
                else:
                    nmd_label = "Truncated protein"
                    impact = "Truncated protein"
                    corruption_position = frameshift_start_codon if frameshift_start_codon is not None else (ptc_c_pos // 3)
                    fraction_lost = max(0.0, 1.0 - corruption_position / prot_len)
                    perc_lost = fraction_lost * 100
                    extra = "<span style='color:orange; font-weight:bold'>Possible driver variant – requires further assessment.</span>"
                    card_color = "🟠"

                # Colored Result Card
                st.markdown(f"""
                <div style="padding:15px; border-radius:8px; border-left:6px solid {'#28a745' if card_color=='🟢' else '#17a2b8' if card_color=='🧬' else '#ffc107'}; background-color:#f8f9fa; margin:10px 0;">
                    <h4>{card_color} {nmd_label}</h4>
                    <p><strong>Impact:</strong> {impact}</p>
                </div>
                """, unsafe_allow_html=True)

                if perc_lost > 0:
                    st.progress(fraction_lost)
                    st.caption(f"**Approx. protein lost:** {perc_lost:.1f}%")

                st.markdown(extra, unsafe_allow_html=True)

                # ====================== S-VIG O2 SUGGESTION ======================
                svig_code, svig_expl, svig_caveat = get_svig_o2_suggestion(
                    ptc_c_pos, prot_len, current["nmd_cutoff_cdna"],
                    frameshift_start_codon, current["gene_tx_key"]
                )

                st.markdown(f"""
                <div style="background-color:#e6f7ff; padding:15px; border-radius:8px; border-left:5px solid #1890ff; margin:12px 0;">
                    <strong>S-VIG O2 suggestion:</strong> <span style="font-size:1.15em; font-weight:bold">{svig_code}</span><br>
                    {svig_expl}
                </div>
                """, unsafe_allow_html=True)

                if svig_caveat:
                    st.warning(svig_caveat)

                # Report data collection
                INPUT_DATA.append({
                    "Variant": line,
                    "Gene": current["gene"],
                    "Transcript": current["transcript"],
                    "PTC codon": ptc_c_pos // 3,
                    "PTC cDNA": f"c.{ptc_c_pos}",
                    "NMD": nmd_label,
                    "Assessment": "DRIVER" if card_color == "🟢" else "Possible driver",
                    "Mechanism / driver note": "NMD" if card_color == "🟢" else "Truncation/Chimera",
                    "Protein impact": impact,
                    "S-VIG O2": svig_code,
                })

with tab_report:
    if not INPUT_DATA:
        st.info("Paste and run variants in the 'Input' tab to generate a structured report.")
    else:
        df = pd.DataFrame(INPUT_DATA)
        st.markdown("### Structured report (click a row to see details)")
        display_df = df[[
            "Variant", "Gene", "Transcript", "NMD", "Assessment",
            "Mechanism / driver note", "Protein impact", "S-VIG O2"
        ]].copy()
        st.dataframe(display_df, use_container_width=True, hide_index=True)
        
        st.divider()
        st.markdown("**CSV‑style summary (for copy‑paste):**")
        st.text(df[["Variant", "Gene", "Transcript", "NMD", "Assessment",
                    "Mechanism / driver note", "Protein impact", "S-VIG O2"]].to_csv(index=False))
        
        csv = df.to_csv(index=False)
        st.download_button(
            label="📥 Download Full Report as CSV",
            data=csv,
            file_name="NMD_Predictor_Report.csv",
            mime="text/csv"
        )

# === Gene Track ===
if INPUT_DATA:
    current = get_params(st.session_state.gene_tx_key)
    prot_len = current["protein_length_aa"]
    nmd_cutoff_aa = current["nmd_cutoff_cdna"] // 3
    
    last = INPUT_DATA[-1]
    variant_label = last["Variant"].split()[-1] if " " in last["Variant"] else last["Variant"]

    codon_ptc = parse_p_ptc_position(last["Variant"].strip())
    ptc_aa = codon_ptc if codon_ptc else 0

    frameshift_start_codon = None
    if "fs" in last["Variant"].lower():
        m_p = re.search(r"p\.[A-Z][a-z]{2}?(\d+)", last["Variant"], re.IGNORECASE)
        if m_p:
            frameshift_start_codon = int(m_p.group(1))

    var_origin_aa = frameshift_start_codon if frameshift_start_codon else ptc_aa
    domains = get_domains(st.session_state.gene_tx_key)

    if domains:
        st.markdown("**Protein Domains (AA positions from UniProt):**")
        st.dataframe(pd.DataFrame(domains)[["name", "start", "end"]], hide_index=True, use_container_width=True)

    st.caption("⚠️ Exon boundaries are approximate visual aids based on NCBI RefSeq for these specific transcripts.")

    fig, ax = plt.subplots(figsize=(14, 5.5), tight_layout=True)
    y = 0
    height = 1.0
    ax.barh(y, prot_len, height=height, color="#f5f5f5", edgecolor="#666", alpha=0.8)

    exons = get_exons(st.session_state.gene_tx_key)
    for start, end, label in exons:
        ax.axvline(start, color="#444444", linestyle="--", linewidth=1.0, alpha=0.6)
        ax.text(start + (end - start)/2, y + 0.5, label,
                ha="center", va="center", fontsize=9, fontweight="bold",
                color="black", bbox=dict(facecolor="white", alpha=0.85, pad=1))

    for d in domains:
        width = d["end"] - d["start"] + 1
        ax.barh(y, width, left=d["start"], height=height,
                color=d.get("color", "#1f77b4"), edgecolor="black", alpha=0.9)
        ax.text(d["start"] + width/2, y + 1.25, d["name"],
                ha="center", va="bottom", fontsize=10, fontweight="bold",
                color="white", bbox=dict(facecolor='black', alpha=0.75, pad=2))

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

# === EDUCATIONAL FACT ===
st.markdown("---")
fact = random.choice(EDUCATIONAL_FACTS)
st.info(f"💡 **Did you know?** {fact}")

# Footer
st.markdown(
    "<p style='font-size:12px; color:#888; text-align:center; margin-top:20px;'>"
    "© 2026 Ashley Sunderland • NMD Predictor (educational use only, no reproduction without permission)."
    "</p>",
    unsafe_allow_html=True,
)
