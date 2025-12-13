from flask import Flask, request, abort, render_template
import os
import sys
import traceback

import promoterDecoder as pdmod
from promoterML import predict_sequence

app = Flask(__name__)


@app.route("/status")
def status():
    return("The Visualisation Test Plugin Flask Server is up and running")


@app.route("/evaluate", methods=["POST"])
def evaluate():
    data = request.get_json(force=True)
    rdf_type = data['type']

    # ~~~~~~~~~~~~ REPLACE THIS SECTION WITH OWN RUN CODE ~~~~~~~~~~~~~~~~~~~
    # uses rdf types
    accepted_types = {'Activity', 'Agent', 'Association', 'Attachment',
                      'Collection', 'CombinatorialDerivation', 'Component',
                      'ComponentDefinition', 'Cut', 'Experiment',
                      'ExperimentalData', 'FunctionalComponent',
                      'GenericLocation', 'Implementation', 'Interaction',
                      'Location', 'MapsTo', 'Measure', 'Model', 'Module',
                      'ModuleDefinition', 'Participation', 'Plan', 'Range',
                      'Sequence', 'SequenceAnnotation', 'SequenceConstraint',
                      'Usage', 'VariableComponent'}

    acceptable = rdf_type in accepted_types

    # # to ensure it shows up on all pages
    # acceptable = True
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~ END SECTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if acceptable:
        return f'The type sent ({rdf_type}) is an accepted type', 200
    else:
        return f'The type sent ({rdf_type}) is NOT an accepted type', 415


@app.route("/run", methods=["POST"])
def run():
    import requests
    import re

    data = request.get_json(force=True)

    top_level_url = data['top_level']
    complete_sbol = data['complete_sbol']
    instance_url = data['instanceUrl']
    size = data['size']
    rdf_type = data['type']
    shallow_sbol = data['shallow_sbol']

    try:
        #Fetch the XML SBOL file
        response = requests.get(complete_sbol, timeout=10)
        if response.status_code != 200:
            return f"<h3>Failed to fetch SBOL file from: {complete_sbol}</h3>", 400

        sbol_text = response.text

        #Extract the sequence from <sbol:elements>...</sbol:elements>
        match = re.search(r"<sbol:elements>(.*?)</sbol:elements>", sbol_text, re.DOTALL)
        if match:
            dnaSeq = match.group(1).strip()
        else:
            dnaSeq = "No <sbol:elements> sequence found in SBOL file."

        #Translate DNA sequence
        dnaSeq = dnaSeq.upper()
        seqSize = (len(dnaSeq)%3 == 0)

        dnaMapColor = {"A":"red", "G":"green", "C":"yellow", "T":"blue"}
        dnaSeqColor = "".join(
            f'<span style="color:{dnaMapColor.get(base, "black")};">{base}</span>'
            for base in dnaSeq)

        rnaSeq = dnaSeq.replace("T", "U")
        aminoMap = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
                    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
                    "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
                    "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
                    
                    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
                    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
                    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
                    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",

                    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
                    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
                    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
                    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",

                    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
                    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
                    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
                    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
        aminoSeq = ""
        for i in range(0, len(rnaSeq) -2, 3):
            codon = rnaSeq[i:i+3]
            aminoSeq += aminoMap.get(codon, "?")

        aminoMapColor = {"F":"lightblue", "L":"lightgreen", "S":"orchid", "Y":"darkorange",
                         "*":"white", "C":"firebrick", "W":"purple", "P":"blue", "H":"red",
                         "Q":"darkorchid", "R":"darkgoldenrod", "I":"yellow", "M":"orange",
                         "T":"green", "N":"mediumblue", "K":"darkgreen", "V":"salmon",
                         "A":"gold", "D":"goldenrod", "E":"chocolate", "G":"darkgrey"}
        aminoSeqColor = "".join(
            f'<span style="color:{aminoMapColor.get(base, "black")};">{base}</span>'
            for base in aminoSeq)

        # Predict promoter using trained model
        ml_label, ml_prob, _ = predict_sequence(dnaSeq)

        
        # Use promoterDecoder.py
        promoter_rows = pdmod.load_promoter_set("promoterSet.tsv")
        pwm_dict = {}
        for row in promoter_rows:
            promoter_name = row.get("2)pmName", "Unknown")
            box_minus10 = row.get("10)boxMinus10seq", "").upper()
            box_minus35 = row.get("12)boxMinus35seq", "").upper()
            if box_minus10 and box_minus10 != "N/A":
                pwm_dict["-10 " + promoter_name] = pdmod.build_pwm([box_minus10])
            if box_minus35 and box_minus35 != "N/A":
                pwm_dict["-35 " + promoter_name] = pdmod.build_pwm([box_minus35])

        #Global promoter scores
        global_score_list = []
        if pwm_dict:
            global_scores = pdmod.compute_global_promoter_strength(dnaSeq, pwm_dict)
            global_score_dict = {}
            for name, score in global_scores:
                rounded = round(score, 2)
                global_score_dict.setdefault(rounded, []).append(name)

            top_scores = sorted(global_score_dict.keys(), reverse=True)[:5]
            global_score_list = [{"score": s, "promoters": global_score_dict[s]} for s in top_scores]

        #Promoter hits
        promoter_hits, consensus_scores = [], []
        if pwm_dict:
            promoter_hits, consensus_scores = pdmod.find_promoter_hits(dnaSeq, pwm_dict, threshold=-10.0)

        sorted_hits = sorted(zip(promoter_hits, consensus_scores), key=lambda x: x[1], reverse=True)
        score_dict = {}
        for hit_str, score in sorted_hits:
            rounded = round(score, 2)
            box_and_promoter, rest = hit_str.split(" at ")
            box, promoter_name = box_and_promoter.split(" ", 1)
            position = rest.split("-")[0]
            entry = {"Promoter": promoter_name, "Box": box, "Position": position}
            score_dict.setdefault(rounded, []).append(entry)

        top_scores = sorted(score_dict.keys(), reverse=True)[:100]
        score_columns = [{"score": s, "hits": score_dict[s]} for s in top_scores]
        max_hits_per_score = max((len(c["hits"]) for c in score_columns), default=0)

        # Flatten global scores (weight by number of promoters)
        all_global_scores = []
        for entry in global_score_list:
            all_global_scores.extend([entry["score"]] * len(entry["promoters"]))

        global_stats = pdmod.compute_score_stats(all_global_scores)

        # Flatten hit scores (weight by number of hits)
        all_hit_scores = []
        for col in score_columns:
            all_hit_scores.extend([col["score"]] * len(col["hits"]))

        hit_stats = pdmod.compute_score_stats(all_hit_scores)

        return render_template(
            "sequence.html",
            top_level_url=top_level_url,
            instance_url=instance_url,
            rdf_type=rdf_type,
            sequence_length=len(dnaSeq),
            dna_seq_color=dnaSeqColor,
            amino_seq_color=aminoSeqColor,
            ml_label=ml_label,
            ml_prob=round(ml_prob, 3),
            global_score_list=global_score_list,
            score_columns=score_columns,
            max_rows=max_hits_per_score,
            global_stats=global_stats,
            hit_stats=hit_stats
        )

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        lnum = exc_tb.tb_lineno
        abort(400, f'Exception: {e}, File: {fname}, line {lnum}')
