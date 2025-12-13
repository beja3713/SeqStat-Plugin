import math
from collections import defaultdict, Counter
import numpy as np
import statistics

def compute_score_stats(score_list):
    if not score_list:
        return {"mean": None, "median": None, "mode": None, "stdev": None, "range": None}

    try:
        mean_val = round(statistics.mean(score_list), 6)
        median_val = round(statistics.median(score_list), 6)
        # If all values are unique, no mode exists
        try:
            mode_val = round(statistics.mode(score_list), 6)
        except statistics.StatisticsError:
            mode_val = "No unique mode"
        stdev_val = round(statistics.stdev(score_list), 6) if len(score_list) > 1 else 0
        range_val = round(max(score_list) - min(score_list), 6)
    except Exception as e:
        mean_val = median_val = mode_val = stdev_val = range_val = str(e)

    return {
        "mean": mean_val,
        "median": median_val,
        "mode": mode_val,
        "stdev": stdev_val,
        "range": range_val
    }

#Import "promotorSet.tsv"

def load_promoter_set(file="promoterSet.tsv"):
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

#Extract motiffs sequences

def extract_motif(rows, box_column):
    motifs = []
    for row in rows:
        motif_get = row.get(box_column)
        if motif_get and motif_get != "N/A":
            motifs.append(motif_get.upper())
    return motifs

#Build pwm

def build_pwm(motifs, pseudocount=1):
    if not motifs:
        return None

    motif_length = len(motifs[0])
    nucleotides = ['A', 'C', 'G', 'T']
    pwm = {n: np.zeros(motif_length) for n in nucleotides}

    # Count nucleotides at each position
    for motif_get in motifs:
        for i, nt in enumerate(motif_get):
            if nt in nucleotides:
                pwm[nt][i] += 1

    # Add pseudocounts and normalize
    for i in range(motif_length):
        total = sum(pwm[n][i] for n in nucleotides) + pseudocount * 4
        for n in nucleotides:
            pwm[n][i] = (pwm[n][i] + pseudocount) / total

    return pwm

#Score DNA sequence

def score_sequence_pwm(seq, pwm, bg=0.25):
    seq = seq.upper()
    score = 0
    for i, nt in enumerate(seq):
        if nt in pwm:
            score += math.log2(pwm[nt][i] / bg)
        else:
            score += math.log2(1e-6 / bg)
    return score


#Promoter Hits

def find_promoter_hits(sequence, pwm_dict, threshold, top_n=200):
    from Bio.Seq import Seq
    hits = []
    scores = []

    seqs_to_scan = [sequence.upper(), str(Seq(sequence).reverse_complement())]

    for seq in seqs_to_scan:
        for name, pwm in pwm_dict.items():
            motif_len = len(next(iter(pwm.values())))
            for i in range(len(seq) - motif_len + 1):
                subseq = seq[i:i+motif_len]
                score = score_sequence_pwm(subseq, pwm)
                if score >= threshold:
                    hits.append(f"{name} at {i}-{i+motif_len}")
                    scores.append(score)

    # Sort by score descending and pick top N
    sorted_hits = sorted(zip(hits, scores), key=lambda x: x[1], reverse=True)
    top_hits = sorted_hits[:top_n]

    # Separate hits and scores for template
    hits_out, scores_out = zip(*top_hits) if top_hits else ([], [])
    scores_out = [round(s, 3) for s in scores_out]

    return hits_out, scores_out

#Global strand score

def compute_global_promoter_strength(sequence, pwm_dict):
    sequence = sequence.upper()
    from Bio.Seq import Seq
    seqs_to_scan = [sequence, str(Seq(sequence).reverse_complement())]

    promoter_scores = {}

    for promoter_name, pwm in pwm_dict.items():
        motif_len = len(next(iter(pwm.values())))
        total_score = 0
        total_frames = 0

        for seq in seqs_to_scan:
            for i in range(len(seq) - motif_len + 1):
                subseq = seq[i:i+motif_len]
                score = score_sequence_pwm(subseq, pwm)
                total_score += score
                total_frames += 1

        global_score = total_score / total_frames if total_frames > 0 else 0
        promoter_scores[promoter_name] = global_score

    # Sort promoters by global score descending and take top 10
    top_global_scores = sorted(promoter_scores.items(), key=lambda x: x[1], reverse=True)[:100]

    return top_global_scores




