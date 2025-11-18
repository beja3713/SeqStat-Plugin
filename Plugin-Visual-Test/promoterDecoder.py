"""
promoterDecoder.py
Streamlined, safe promoter scoring module (RegulonDB-aware, sigma-scanning).
Handles any input sequence and returns structured promoter hits.
"""

import pandas as pd
import csv

#Import "promotorSet.tsv"

def load_promoter_set(file="promotorSet.tsv"):
    rows = []

    with open(file, "r", encoding="utf-8") as f:
        lines = [line.rstrip("\n") for line in f]

    for line in lines:
        if not line.startswith("#"):
            header = line.split("\t")
            break

    for line in lines:
        if line.startswith("#") or line == "" or line == "\t".join(header):
            continue

        parts = line.split("\t")

        if len(parts) < len(header):
            parts += ["N/A"] * (len(header) - len(parts))

        row = {header[i]: parts[i] if parts[i] else "N/A"
               for i in range(len(header))}
        rows.append(row)

    return rows


#Sigma consensus motifs

SIGMA_CONSENSUS = {
    "sigma70": {"10": "TATAAT", "35": "TTGACA", "spacing": (16, 18)},
    "sigma38": {"10": "TATACT", "35": "CTTGAA", "spacing": (16, 18)},
    "sigma32": {"10": "CCCCAT", "35": "TTGAAA", "spacing": (12, 15)},
    "sigma24": {"10": "GGAACT", "35": "GAACTT", "spacing": (14, 19)},
    "sigma28": {"10": "TAAA", "35": "CTAAA", "spacing": (13, 15)},
    "sigma54": {"10": None, "35": None, "spacing": None}
}

#Utility functions

def similarity(seq1, seq2):
    """Returns similarity fraction between two sequences safely."""
    try:
        seq1 = str(seq1).upper()
        seq2 = str(seq2).upper()
        if len(seq1) != len(seq2) or not seq1:
            return 0
        return sum(c1 == c2 for c1, c2 in zip(seq1, seq2)) / len(seq1)
    except Exception:
        return 0

def spacing_score(pos35, pos10, ideal_range):
    """Scores spacing between -35 and -10 motifs."""
    try:
        if pos35 is None or pos10 is None or ideal_range is None:
            return 0
        observed = int(pos10) - int(pos35)
        low, high = ideal_range
        if low <= observed <= high:
            return 1
        return max(0, 1 - 0.1 * abs(observed - (low + high)/2))
    except Exception:
        return 0

def score_against_sigma(box10, box35, pos10, pos35, sigma_name):
    cons = SIGMA_CONSENSUS[sigma_name]
    if sigma_name == "sigma54":
        return {"sigma": sigma_name, "score": 0, "10_sim": 0, "35_sim": 0, "spacing": 0}
    s10 = similarity(box10, cons["10"])
    s35 = similarity(box35, cons["35"])
    spc = spacing_score(pos35, pos10, cons["spacing"])
    return {"sigma": sigma_name, "score": s10*2 + s35*2 + spc, "10_sim": s10, "35_sim": s35, "spacing": spc}

def score_all_sigmas(row):
    """Scores a single promoter row against all sigma factors."""
    box10 = row.get("boxMinus10seq") if pd.notna(row.get("boxMinus10seq")) else None
    box35 = row.get("boxMinus35seq") if pd.notna(row.get("boxMinus35seq")) else None
    pos10 = pos35 = None
    try:
        pos10 = int(str(row.get("boxMinus10pos")).split("-")[0])
    except Exception:
        pass
    try:
        pos35 = int(str(row.get("boxMinus35pos")).split("-")[0])
    except Exception:
        pass

    results = [score_against_sigma(box10, box35, pos10, pos35, sigma) for sigma in SIGMA_CONSENSUS]
    ranked = sorted(results, key=lambda x: x["score"], reverse=True)
    return {"ranked_sigma_matches": ranked, "best_sigma": ranked[0]["sigma"], "best_score": ranked[0]["score"]}

#Promoter dataset handling

def load_regulondb_table(promoters_file):
    columns = [
        "pmId", "pmName", "strand", "posTSS", "sigmaFactor", "pmSequence",
        "firstGeneName", "distToFirstGene", "boxMinus10pos", "boxMinus10seq",
        "boxMinus35pos", "boxMinus35seq", "pmEvidence", "addEvidence",
        "confidenceLevel", "pmids"
    ]
    df = pd.read_csv(promoters_file, sep="\t", comment="#", engine="python", dtype=str, on_bad_lines='skip')
    for col in columns:
        if col not in df.columns:
            df[col] = pd.NA
    df = df.iloc[:, :len(columns)]
    df.columns = columns
    df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    return df

# ---------------------- Promoter scoring interface ------------------------

def get_promoter_hits(promoter_df, top_n=5):
    """Return top N scored promoters with sigma factor and associated gene info."""
    hits = []
    for _, row in promoter_df.iterrows():
        sigma_scores = score_all_sigmas(row)
        top_sigma = sigma_scores["ranked_sigma_matches"][0]
        hits.append({
            "promoter_name": row.get("pmName", "unknown"),
            "best_sigma": top_sigma["sigma"],
            "best_score": sigma_scores["best_score"],
            "box10_sim": top_sigma["10_sim"],
            "box35_sim": top_sigma["35_sim"],
            "associated_gene": row.get("firstGeneName", "")
        })
    hits_sorted = sorted(hits, key=lambda r: r["best_score"], reverse=True)
    return hits_sorted[:top_n]

# ---------------------- DNA scanning (optional) ---------------------------

def scan_dna_for_promoters(dna_seq, promoter_df):
    """
    Scan a DNA sequence for promoter motifs. Returns top matches.
    """
    dna_seq = dna_seq.upper()
    hits = []
    for _, row in promoter_df.iterrows():
        box10 = row.get("boxMinus10seq")
        box35 = row.get("boxMinus35seq")
        if pd.isna(box10) or pd.isna(box35):
            continue
        box10, box35 = box10.upper(), box35.upper()
        window_size = len(box10) + len(box35) + 18
        best_score = -1
        best_pos = 0
        for i in range(len(dna_seq) - window_size + 1):
            window = dna_seq[i:i+window_size]
            s10 = similarity(window[:len(box10)], box10)
            s35 = similarity(window[-len(box35):], box35)
            score = s10 + s35
            if score > best_score:
                best_score = score
                best_pos = i
        hits.append({"promoter_name": row.get("pmName", "unknown"), "best_score": best_score, "position": best_pos})
    return sorted(hits, key=lambda r: r["best_score"], reverse=True)[:5]

def score_all_sigmas_on_sequence(dna_seq):
    results = []
    for sigma_name, sigma_motif in sigma_consensus_dict.items():
        box10_sim, box35_sim = score_sigma_motifs(dna_seq, sigma_motif)
        results.append({
            "sigma": sigma_name,
            "box10_sim": box10_sim,
            "box35_sim": box35_sim
        })
    return results

