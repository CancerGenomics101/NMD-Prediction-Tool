import streamlit as st
import re


# === TRANSCRIPT DATABASE ==================================================================

TRANSCRIPTS = {
    "ASXL1_NM_015338.5": (
        "NM_015338.5",
        1669,          # c.1669 → NMD cutoff
        4623,          # CDS ≈ 3 * 1541
        1541,
    ),
    "TET2_NM_001127208.2": (
        "NM_001127208.2",
        4487,          # c.4487 → NMD cutoff
        6006,          # CDS ≈ 3 * 2002
        2002,
    ),
}


def get_params(gene_tx_key):
    tx, nmd_cut, mrna_len, prot_len = TRANSCRIPTS[gene_tx_key]
    return {
        "gene_tx_key": gene_tx_key,
        "gene": gene_tx_key.split("_")[0],
        "transcript": tx,
        "nmd_cutoff_cdna": nmd_cut,
        "reference_mrna_len": mrna_len,  # CDS for % protein
        "protein_length_aa": prot_len,
    }


# === HGVS PARSING (now more flexible) =====================================================

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
    hgvs_p = hgvs_p.strip()

    # p.Ser156* / p.Ser156Ter
    m = re.search(r"p\.[A-Z][a-z]{2}(\d+)(?:\*|Ter)", hgvs_p, re.IGNORECASE)
    if m:
        aa_stop = int(m.group(1))
        return aa_stop

    # p.Trp343Alafs*39
    m = re.search(r"p\.[A-Z][a-z]{2}(\d+)[A-Z][a-z]{2}?fs(?:Ter|\*)?(\d+)", hgvs_p, re.IGNORECASE)
    if m:
        aa_start = int(m.group(1))
        n_aa_new = int(m.group(2))
        ptc_codon = aa_start + n_aa_new - 1
        return ptc_codon

    return None


def hgvs_to_ptc_c_pos(hgvs_str):
    hgvs_str = hgvs_str.strip()
    if not hgvs_str:
        return None, "Empty string"

    # c. part (now handles del/delins/ins/dup)
    c_match = re.search(r"c\.\d+.*", hgvs_str, re.IGNORECASE)
    if not c_match:
        return None, "No c. part found"
    c_str = c_match.group()
    c_start = extract_c_pos_from_c_hgvs(c_str)
    if c_start is None:
        return None, "Failed to parse c. position"

    # p. part
    p_match = re.search(r"p\.[^ ]*", hgvs_str, re.IGNORECASE)
    if not p_match:
        return None, "No p. part found"
    p_str = p_match.group()

    ptc_codon = parse_p_ptc_position(p_str)
    if ptc_codon is None:
        return None, "Failed to parse p. frameshift/stop"

    # codon → cDNA (first base of the codon)
    ptc_c_pos = 3 * ptc_codon
    return ptc_c_pos, None


# === STREAMLIT LAYOUT ====================================================================

st.title("NMD Predictor (HGVS input, GRCh37)")

st.selectbox(
    "Select gene and transcript:",
    key="gene_tx_key",
    options=list(TRANSCRIPTS.keys()),
    format_func=lambda x: x.replace("_", " / "),
)

current = get_params(st.session_state.gene_tx_key)

st.markdown(f"""
You are currently using:
- **Gene:** {current['gene']}  
- **Transcript:** {current['transcript']}  
- **NMD cutoff:** cDNA position ≤ {current['nmd_cutoff_cdna']} → NMD predicted  
- **Protein length:** {current['protein_length_aa']} aa  
(CDS ≈ {current['reference_mrna_len']} bp, no premature stop.)
""")

input_text = st.text_area(
    "Paste HGVS variant (c. and p.):",
    value="c.994A>G p.Trp343Alafs*39",
    height=100,
)

if not input_text.strip():
    st.info("Paste a variant in the box above.")
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
        cds_len = current["reference_mrna_len"]
        prot_len = current["protein_length_aa"]

        cds_end = 3 * prot_len

        if ptc_c_pos <= current["nmd_cutoff_cdna"]:
            nmd = "YES (NMD predicted)"
            impact = "Full loss (NMD) – 100% of protein lost"
            fraction_lost = 1.0
            perc_lost = 100.0
            extra = "<span style='color:green; font-weight:bold'>DRIVER</span>"
        elif ptc_c_pos > cds_end:
            # PTC in 3′ UTR → extension / chimera‑like
            nmd = "NO (extended / chimera‑like)"
            impact = "Extended protein (3′ UTR PTC)"
            ratio = ptc_c_pos / cds_len
            fraction_lost = 1.0 - ratio
            fraction_lost = max(0.0, fraction_lost)
            perc_lost = fraction_lost * 100
            extra = (
                "<span style='color:orange; font-weight:bold'>"
                "Chimera‑like construct generated"
                "</span>"
            )
        else:
            # Truncated protein (within CDS but after NMD cutoff)
            nmd = "NO (truncated protein)"
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

        st.markdown(f"**NMD?:** {nmd}")
        st.markdown(f"**Protein impact:** {impact}")
        st.markdown(f"**Approx. protein lost:** {fraction_lost:.2f} ({perc_lost:.1f}%)")
        st.markdown(extra, unsafe_allow_html=True)
