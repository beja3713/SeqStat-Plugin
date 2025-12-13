# promoter_training.py
import pandas as pd
import numpy as np
import random
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, roc_auc_score
import joblib

#Load the RegulonDB TSV
tsv_file = "promoterSet.tsv"  # Replace with your file path
promoter_df = pd.read_csv(tsv_file, sep="\t", comment="#")

# Inspect columns
print("Columns:", promoter_df.columns.tolist())
print("Total promoters:", len(promoter_df))

#Filter for strong/confirmed promoters
promoter_df = promoter_df[promoter_df['15)confidenceLevel'].isin(['S', 'C'])]
print("Filtered promoters:", len(promoter_df))

#Prepare positive and negative sequences
positive_seqs = promoter_df['6)pmSequence'].tolist()

# Generate negative sequences by shuffling positive sequences
def shuffle_sequence(seq):
    if not isinstance(seq, str):
        return ""
    import random
    seq = list(seq)
    random.shuffle(seq)
    return ''.join(seq)

negative_seqs = [shuffle_sequence(seq) for seq in positive_seqs if isinstance(seq, str)]

# Combine sequences and labels
X_seqs = positive_seqs + negative_seqs
y = [1]*len(positive_seqs) + [0]*len(negative_seqs)

print(f"Total sequences: {len(X_seqs)}, Positives: {len(positive_seqs)}, Negatives: {len(negative_seqs)}")

#One-hot encode sequences
def one_hot_encode(seq):
    mapping = {'A':[1,0,0,0], 'T':[0,1,0,0], 'C':[0,0,1,0], 'G':[0,0,0,1]}
    return np.array([mapping.get(base, [0,0,0,0]) for base in seq.upper()])

# Flatten sequences for scikit-learn
# Filter out bad sequences and corresponding labels
filtered_data = [(seq, label) for seq, label in zip(X_seqs, y) if isinstance(seq, str) and seq.strip() != ""]
X_seqs, y = zip(*filtered_data)  # unzip into sequences and labels
X = np.array([one_hot_encode(seq).flatten() for seq in X_seqs])
y = np.array(y)
print("Encoded shape:", X.shape)

#Split into train/test
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

#Train Random Forest
model = RandomForestClassifier(n_estimators=200, max_depth=10, random_state=42)
model.fit(X_train, y_train)

#Evaluate
y_pred = model.predict(X_test)
y_prob = model.predict_proba(X_test)[:,1]

print("\nClassification Report:\n", classification_report(y_test, y_pred))
print("ROC-AUC:", roc_auc_score(y_test, y_prob))

#Save model
model_file = "promoter_model.pkl"
joblib.dump(model, model_file)
print(f"Model saved to {model_file}")
