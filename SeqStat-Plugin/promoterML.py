import numpy as np
import joblib
import os

# Load your pre-trained Random Forest model
model_file = "promoter_model.pkl"
if not os.path.exists(model_file):
    raise FileNotFoundError(f"Model file not found: {model_file}")
rf_model = joblib.load(model_file)

# One-hot encoding function
def one_hot_encode(seq):
    mapping = {'A':[1,0,0,0], 'T':[0,1,0,0], 'C':[0,0,1,0], 'G':[0,0,0,1]}
    return np.array([mapping.get(base.upper(), [0,0,0,0]) for base in seq])

def predict_sequence(dna_seq, window_size=81, step_size=1, threshold=0.5):
    dna_seq = dna_seq.upper()
    n = len(dna_seq)
    if n < window_size:
        # pad short sequences with 'N's to match window size
        dna_seq = dna_seq + "N" * (window_size - n)
        n = len(dna_seq)

    max_prob = -1
    max_label = "Non-promoter"
    hits = []

    for start in range(0, n - window_size + 1, step_size):
        window = dna_seq[start:start + window_size]
        X_window = one_hot_encode(window).flatten().reshape(1, -1)

        # Predict
        prob = rf_model.predict_proba(X_window)[0][1]  # probability of promoter
        label = "Promoter" if prob >= threshold else "Non-promoter"

        hits.append({"start": start, "end": start + window_size, "prob": prob, "label": label})

        if prob > max_prob:
            max_prob = prob
            max_label = label

    return max_label, max_prob, hits
