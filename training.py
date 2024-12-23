"""
Training script. Generates transition and emission matricies.
"""

from dotenv import load_dotenv
import os
import numpy as np
import os
from collections import defaultdict, deque

# Load environment variables
load_dotenv()
EXON_FILE = os.getenv('OUTPUT_EXON_FILEPATH')
DNA_FILE = os.getenv('FASTA_FILEPATH')
NUM_TRAINING_LABELS = int(os.getenv('NUM_TRAINING_LABELS'))

NUM_STATES = 8
NUM_EMISSIONS = 4

class Labels:
    INTERGENIC = 0
    INTRON_PHASE_0 = 1
    INTRON_PHASE_1 = 2
    INTRON_PHASE_2 = 3
    INITIAL_EXON = 4
    FINAL_EXON = 5
    INTERNAL_EXON = 6
    SINGLE_EXON = 7

class Neucleotides:
    A = 0
    C = 1
    G = 2
    T = 3

class HMM:
    def __init__(self, num_states, num_observations, order):
        self.num_states = num_states
        self.transition_matrix = np.zeros([num_states] * (order))
        self.emission_matrix = np.zeros((num_states, num_observations))
        self.num_observations = num_observations
        self.smoothing = 1e-6
        self.order = order - 1

    def train(self, sequences, labels):
        for seq, lbl in zip(sequences, labels):
            for i in range(len(lbl) - self.order):
                prev_states = tuple(lbl[i:i + self.order])
                next_state = lbl[i + self.order]
                self.transition_matrix[(*prev_states, next_state)] += 1

            for i, base in enumerate(seq):
                current_state = lbl[i]
                self.emission_matrix[current_state, base] += 1

        transition_sums = self.transition_matrix.sum(axis=-1, keepdims=True) + self.smoothing * self.num_states
        self.transition_matrix = (self.transition_matrix + self.smoothing) / transition_sums

        for state in range(self.num_states):
            total = self.emission_matrix[state].sum() + self.smoothing * self.num_observations
            self.emission_matrix[state] = (self.emission_matrix[state] + self.smoothing) / total


    def trainLabel(self, sequences, labels):
        for seq, lbl in zip(sequences, labels):
            unique_queue = deque(maxlen=self.order)

            for i in range(len(lbl)):
                current_label = lbl[i]

                for prev_state in unique_queue:
                    print([prev_state, current_label])
                    self.transition_matrix[prev_state, current_label] += 1
                    

                if not unique_queue or unique_queue[-1] != current_label:
                    unique_queue.append(current_label)

                base = seq[i]
                self.emission_matrix[current_label, base] += 1

        transition_sums = self.transition_matrix.sum(axis=-1, keepdims=True) + self.smoothing * self.num_states
        self.transition_matrix = (self.transition_matrix + self.smoothing) / transition_sums

        for state in range(self.num_states):
            total = self.emission_matrix[state].sum() + self.smoothing * self.num_observations
            self.emission_matrix[state] = (self.emission_matrix[state] + self.smoothing) / total


def parse_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as f:
        current_seq_id = None
        current_seq = []

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq_id is not None and current_seq_id.isdigit():
                    sequences[current_seq_id] = ''.join(current_seq)
                current_seq_id = line[1:3].replace(" ", "")
                current_seq = []
                
            else:
                current_seq.append(line)

    return sequences


def parse_exons(file_path):
    exons = defaultdict(list)
    currExons = []
    prev_seq_id = ''

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                exons[prev_seq_id].append(currExons)
                currExons = []
                continue

            seq_id, start, end = line.split()
            start, end = int(start), int(end)
            prev_seq_id = seq_id
            currExons.append((start, end))

    exons[prev_seq_id].append(currExons)
    return exons


def build_training_data(sequences, exons):
    region_data = {'sequences': [], 'labels': []}

    base_to_int = {"A": Neucleotides.A, "C": Neucleotides.C, "G": Neucleotides.G, "T": Neucleotides.T}

    for seq_id, seq in sequences.items():
        seq_length = len(seq)
        seq_labels = [0] * seq_length
        seq_integers = [base_to_int[base] for base in seq if base in base_to_int]

        for sublist in exons.get(seq_id, []):
            valid_exons = [(start, end) for start, end in sublist if start <= end]

            if len(valid_exons) == 1:
                start, end = valid_exons[Labels.INTERGENIC]
                for i in range(start - 1, end):
                    if 0 <= i < seq_length:
                        seq_labels[i] = Labels.SINGLE_EXON
            elif len(valid_exons) > 1:
                for i, (start, end) in enumerate(valid_exons):
                    for j in range(start - 1, end):
                        if 0 <= j < seq_length:
                            if i == 0:
                                seq_labels[j] = Labels.INITIAL_EXON
                            elif i == len(valid_exons) - 1:
                                seq_labels[j] = Labels.FINAL_EXON 
                            else:
                                seq_labels[j] = Labels.INTERNAL_EXON

            for i in range(len(valid_exons) - 1):
                _, end_of_current_exon = valid_exons[i]
                start_of_next_exon, _ = valid_exons[i + 1]
                for j in range(end_of_current_exon, start_of_next_exon - 1): 
                    if 0 <= j < seq_length:
                        intron_pos = j - end_of_current_exon
                        if intron_pos == 0:
                            seq_labels[j] = Labels.INTRON_PHASE_0
                        elif intron_pos == 1:
                            seq_labels[j] = Labels.INTRON_PHASE_1
                        elif intron_pos == 2:
                            seq_labels[j] = Labels.INTRON_PHASE_2

        region_data['sequences'].append(seq_integers)
        region_data['labels'].append(seq_labels)

    return region_data

if __name__ == "__main__":
    current_dir = os.getcwd()
    fasta_file = os.path.join(current_dir, DNA_FILE)
    exon_file = os.path.join(current_dir, EXON_FILE)

    # Parse & Structure data
    sequences = parse_fasta(fasta_file)
    exons = parse_exons(exon_file)
    region_data = build_training_data(sequences, exons)

    print('Starting training...')

    training_sequences = region_data['sequences'][:NUM_TRAINING_LABELS]
    training_labels = region_data['labels'][:NUM_TRAINING_LABELS]

    for i in range(7):
        model = HMM(NUM_STATES, NUM_EMISSIONS, i + 1)
        model.train(training_sequences, training_labels)
        np.save(f'matrices/transition_matrix_{i}', np.log(model.transition_matrix))
        print(f'Model {i} complete.')	
    
    np.save('matrices/emission_matrix',np.log(model.emission_matrix))
    print('Saved all .npy files.')

