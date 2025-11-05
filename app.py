from flask import Flask, request, abort
import os
import sys
import traceback


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
            sequence = match.group(1).strip()
        else:
            sequence = "No <sbol:elements> sequence found in SBOL file."

        #Translate DNA sequence
        sequence = sequence.upper()
        seqStart = sequence.startswith("ATG")
        codon = 'Met'
        seqSize = (len(sequence)%3 == 0)

        colorMap = {"A": "red", "G": "green", "C": "yellow", "T": "blue"}
        colorSeq = "".join(
            f'<span style="color:{colorMap.get(base, "black")};">{base}</span>'
            for base in sequence)
        

        #Build an HTML page to display the sequence
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <style>
                body {{ font-family: monospace; margin: 2em; background: #f8f9fa; color: #333; }}
                h1 {{ font-family: Arial, sans-serif; }}
                pre {{
                    background: #fff;
                    border: 1px solid #ccc;
                    padding: 1em;
                    border-radius: 5px;
                    overflow-x: auto;
                    max-width: 100%;
                    white-space: pre-wrap;
                    word-wrap: break-word;
                }}
                .meta {{
                    font-family: Arial, sans-serif;
                    margin-bottom: 1em;
                }}
            </style>
        </head>
        <body>
            <h1>Gene Sequence Visualization</h1>
            <div class="meta">
                <p><b>Top Level URL:</b> <a href="{top_level_url}">{top_level_url}</a></p>
                <p><b>Instance:</b> <a href="{instance_url}">{instance_url}</a></p>
                <p><b>Type:</b> {rdf_type}</p>
                <p><b>Sequence Length:</b> {len(sequence)} bases</p>
            </div>
            <pre>{colorSeq}</pre>
        </body>
        </html>
        """

        return html

    except Exception as e:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        lnum = exc_tb.tb_lineno
        abort(400, f'Exception: {e}, File: {fname}, line {lnum}')
