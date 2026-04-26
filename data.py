# data.py - NMD Predictor Gene Library
# All gene-specific data lives here, requires last c. position for NMD, final coding nucleotide c. position, and final coding AA p. position
# Protein domains are a static pull from UniProt matched to the given transcript and approximate exon boundaries estimated from NCBI RefSeq for visual aid only

TRANSCRIPTS = {
    "ASXL1_NM_015338.5": (
        "NM_015338.5", 1669, 4623, 1541,
    ),
    "TET2_NM_001127208.2": (
        "NM_001127208.2", 4487, 6006, 2002,
    ),
    "DNMT3A_NM_175629.2": (
        "NM_175629.2", 2547, 4395, 912,
    ),
    "TP53_NM_000546.5": (
        "NM_000546.5", 1050, 2591, 393,
    ),
}

PROTEIN_DOMAINS = {
    "ASXL1_NM_015338.5": [
        {"name": "HARE-HTH", "start": 11,  "end": 86,   "color": "#4CAF50"},
        {"name": "DEUBAD",   "start": 255, "end": 364,  "color": "#2196F3"},
        {"name": "PHD",      "start": 1503,"end": 1540, "color": "#FF9800"},
    ],
    "TET2_NM_001127208.2": [
        {"name": "Cys-rich", "start": 1130, "end": 1300, "color": "#9C27B0"},
        {"name": "Tet_JBP",  "start": 1300, "end": 2002, "color": "#FF5722"},
    ],
    "DNMT3A_NM_175629.2": [
        {"name": "PWWP",  "start": 280, "end": 367, "color": "#4CAF50"},
        {"name": "ADD",   "start": 600, "end": 712, "color": "#2196F3"},
        {"name": "MTase", "start": 634, "end": 912, "color": "#FF9800"},
    ],
    "TP53_NM_000546.5": [
        {"name": "Transactivation", "start": 1,   "end": 61,  "color": "#4CAF50"},
        {"name": "DNA-binding",     "start": 94,  "end": 312, "color": "#2196F3"},
        {"name": "Tetramerisation", "start": 325, "end": 356, "color": "#FF9800"},
    ],
}

EXONS = {
    "ASXL1_NM_015338.5": [
        (1,19,"ex1"), (20,47,"ex2"), (48,48,"ex3"), (49,84,"ex4"),
        (85,124,"ex5"), (125,157,"ex6"), (158,188,"ex7"), (189,239,"ex8"),
        (240,294,"ex9"), (295,327,"ex10"), (328,362,"ex11"), (363,573,"ex12"),
        (574,1541,"ex13"),
    ],
    "TET2_NM_001127208.2": [
        (1,60,"ex1"), (61,117,"ex2"), (118,173,"ex3"), (174,240,"ex4"),
        (241,317,"ex5"), (318,393,"ex6"), (394,473,"ex7"), (474,550,"ex8"),
        (551,607,"ex9"), (608,1496,"ex10"), (1497,2002,"ex11"),
    ],
    "DNMT3A_NM_175629.2": [
        (1,120,"ex1"), (121,220,"ex2"), (221,300,"ex3"), (301,380,"ex4"),
        (381,460,"ex5"), (461,540,"ex6"), (541,620,"ex7"), (621,700,"ex8"),
        (701,780,"ex9"), (781,860,"ex10"), (861,912,"ex11-23"),
    ],
    "TP53_NM_000546.5": [
        (1,42,"ex1"), (43,82,"ex2"), (83,93,"ex3"), (94,140,"ex4"),
        (141,186,"ex5"), (187,222,"ex6"), (223,260,"ex7"), (261,292,"ex8"),
        (293,312,"ex9"), (313,356,"ex10"), (357,393,"ex11"),
    ],
}

EDUCATIONAL_FACTS = [
    "Nonsense-mediated decay (NMD) typically degrades transcripts with premature termination codons more than 50–55 nucleotides upstream of the last exon-exon junction.",
    "The 'last exon rule' means truncating variants in the final exon often escape NMD and produce truncated proteins.",
    "Tumor suppressor genes are particularly sensitive to loss-of-function variants, making accurate NMD prediction clinically important.",
    "Frameshift variants near the C-terminus are less likely to trigger NMD than those early in the protein.",
    "The position of a variant relative to key protein domains is often more biologically important than the exact percentage of protein lost.",
]
