import streamlit as st
import re


GENE = "ASXL1"
TRANSCRIPT = "NM_015338.5"
NMD_CUTOFF_CDNA = 1669             # cDNA ≤ 1669 → NMD
REFERENCE_MRNA_LEN = 4623          # cDNA length
PROTEIN_LENGTH_AA = 1541           # canonical protein length


def extract_c_pos_from_c_hgvs(hgvs_c):
    """Given c.994A>G, return 994."""
    m = re.search(r"c\.(\d+)", hgvs_c)
    if m:
        return int(m.group(1))
    return None


def parse_p_ptc_position(hgvs_p):
    """
    From p. text, return the PTC codon number (e.g., 156).
    Supports:
      - p.Ser156*  / p.Ser156Ter
      - p.Trp343Alafs*39 / p.Trp343AlafsTer40
    """
    hgvs_p = hgvs_p.strip()

    # Case 1: simple nonsense (p.Ser156* or p.Ser156Ter)
    m = re.search(r"p\.[A-Z][a-z]{2}(\d+)(?:\*|Ter)", hgvs_p, re.IGNORECASE)
    if m:
        aa_stop = int(m.group(1))
        return aa_stop

    # Case 2: frameshift with new stop (p.Trp343Alafs*39)
    m = re.search(r"p\.[A-Z][a-z]{2}(\d+)[A-Z][a-z]{2}?fs(?:Ter|\*)?(\d+)", hgvs_p, re.IGNORECASE)
    if m:
        aa_start = int(m.group(1))
        n_aa_new = int(m.group(2))
        ptc_codon = aa_start + n_aa_new - 1
        return ptc_codon

    return None


def hgvs_to_ptc_c_pos(hgvs_str):
    """
    From a string like "c.994A>G p.Trp343Alafs*39" or "c.466T>G p.Ser156*",
    return the PTC cDNA position (or None if parsing fails).
    """
    hgvs_str = hgvs_str.strip()
    if not hgvs_str:
        return None, "Empty string"

    # Extract c. part
    c_match = re.search(r"c\.\d+[ACGT>][a-zA-Z0-9]*", hgvs_str, re.IGNORECASE)
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

    # codon → cDNA (first base of the codon)
    ptc_c_pos = 3 * ptc_codon
    return ptc_c_pos, None


st.title("ASXL1 NMD Predictor (HGVS input, NM_015338.5)")

st.markdown(f"""
For **ASXL1, transcript {TRANSCRIPT}**, paste an HGVS‑style variant.

- **NMD cutoff:** cDNA position ≤ {NMD_CUTOFF_CDNA} → NMD.
- **Canonical protein length:** {PROTEIN_LENGTH_AA} aa.

Examples:
- `c.994A>G p.Trp343Alafs*39`  → frameshift, PTC at codon 381
- `c.466T>G p.Ser156*`         → nonsense, PTC at codon 156
- `c.1234A>T p.Trp412Ter`      → nonsense, PTC at codon 412
- `c.1800C>T p.Arg600fs*60`    → frameshift, PTC at codon 659
""")

input_text = st.text_area(
    "Paste HGVS variant (c. and p.):",
    value="c.467C>A p.Ser156*",
    height=100,
)

if not input_text.strip():
    st.info("Paste a variant in the box above.")
else:
    # Split into lines (allow multiple variants)
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

        # NMD decision
        if ptc_c_pos <= NMD_CUTOFF_CDNA:
            nmd = "YES (NMD predicted)"
            fraction_lost = 1.0
            impact = "Full loss (NMD)"
        else:
            nmd = "NO (truncated protein)"
            ratio = ptc_c_pos / REFERENCE_MRNA_LEN
            fraction_lost = 1.0 - ratio
            fraction_lost = max(0.0, fraction_lost)
            impact = "Truncated protein"

        perc_lost = fraction_lost * 100

        st.markdown(f"**NMD?:** {nmd}")
        st.markdown(f"**Protein impact:** {impact}")
        st.markdown(f"**Approx. protein lost:** {fraction_lost:.2f} ({perc_lost:.1f}%)")
