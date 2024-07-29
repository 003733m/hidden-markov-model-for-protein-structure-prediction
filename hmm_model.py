import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def read_tsv(file_path):
    sequences = []
    labels = []

    with open(file_path, 'r') as file: #Reading the tsv file here.
        for line in file:
            if line.startswith("Entry"):
                continue
            data = line.strip().split('\t')
            if len(data) >= 4:
                sequence = data[2]
                label_info = data[3:]

                sequences.append(sequence)
                labels.append(label_info)

    return sequences, labels
def plot_hmm(transition_matrix, emission_matrix):
    states = ['Helix', 'Strand', 'Turn', 'Unknown']

    G = nx.DiGraph()

    # Add nodes
    for state in states:
        G.add_node(state)

    # Add edges
    for i, from_state in enumerate(states):
        for j, to_state in enumerate(states):
            weight = transition_matrix[i, j]
            G.add_edge(from_state, to_state, weight=weight)

    # Position nodes
    pos = {'Helix': (0, 1), 'Strand': (1, 1), 'Turn': (0, 0), 'Unknown': (1, 0)}

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=700, node_color='skyblue')

    # Draw edges
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=weights, edge_color='gray', arrowsize=20)

    # Label nodes
    nx.draw_networkx_labels(G, pos, font_size=10, font_color='black')

    plt.title('HMM Transition Diagram')
    plt.show()

class HMM:
    def __init__(self, num_states, num_symbols):
        self.num_states = num_states
        self.num_symbols = num_symbols

        # Initializing transition and emission matrices with random values
        self.transition_matrix = np.random.rand(num_states, num_states)
        self.transition_matrix /= np.sum(self.transition_matrix, axis=1, keepdims=True)

        self.emission_matrix = np.random.rand(num_states, num_symbols)
        self.emission_matrix /= np.sum(self.emission_matrix, axis=1, keepdims=True)

    def forward(self, observations):
        T = len(observations)
        alpha = np.zeros((T, self.num_states))

        # Initialization
        alpha[0] = self.emission_matrix[:, observations[0]]

        # Forward algo
        for t in range(1, T):
            alpha[t] = np.dot(alpha[t - 1], self.transition_matrix) * self.emission_matrix[:, observations[t]]

        return alpha

    def backward(self, observations):
        T = len(observations)
        beta = np.zeros((T, self.num_states))

        # Initialization
        beta[-1] = 1.0

        # Backward algo
        for t in range(T - 2, -1, -1):
            beta[t] = np.dot(self.transition_matrix, self.emission_matrix[:, observations[t + 1]] * beta[t + 1])

        return beta

    def update_params(self, observations, alpha, beta):
        T = len(observations)
        xi = np.zeros((T - 1, self.num_states, self.num_states))
        gamma = alpha * beta

        for t in range(T - 1):
            numerator = alpha[t][:, np.newaxis] * self.transition_matrix * self.emission_matrix[:,
                                                                           observations[t + 1]] * beta[t + 1]
            denominator = np.sum(np.sum(numerator, axis=0), axis=0)
            denominator=denominator+20

            if np.any(denominator == 0):
                # Handle the case where the denominator is zero
                self.transition_matrix = np.random.rand(self.num_states, self.num_states)
                self.transition_matrix /= np.sum(self.transition_matrix, axis=1, keepdims=True)
                return

            xi[t] = (numerator + 1) / (denominator + self.num_states)  # Add pseudo-count of 1 to numerator


    def train(self, observations, max_iters=100, tol=1e-6):
        norm_changes = []  # to store norm changes between iterations
        for iter in range(max_iters):
            if iter == 0:  # Print results only for the first iteration
                print(f'Iteration {iter + 1}, Transition Matrix:')
                print(self.transition_matrix)
                print('Emission Matrix:')
                print(self.emission_matrix)

            old_transition_matrix = np.copy(self.transition_matrix)
            old_emission_matrix = np.copy(self.emission_matrix)

            # E-step
            alpha = self.forward(observations)
            beta = self.backward(observations)

            # M-step
            self.update_params(observations, alpha, beta)

            # Check for convergence
            transition_norm = np.linalg.norm(self.transition_matrix - old_transition_matrix)
            emission_norm = np.linalg.norm(self.emission_matrix - old_emission_matrix)
            norm_changes.append((transition_norm, emission_norm))

            if iter == 0:  # Print results only for the first iteration
                print(f'Updated Transition Matrix:')
                print(self.transition_matrix)
                print('Updated Emission Matrix:')
                print(self.emission_matrix)

            if transition_norm < tol and emission_norm < tol:
                break

        return norm_changes



    def predict(self, test_sequence): # Predicting the next protein sequence
        numerical_test_sequence = np.array([ord(aa) - ord('A') for aa in test_sequence]) # It's processing with A value of ASCII

        alpha_test = self.forward(numerical_test_sequence)
        beta_test = self.backward(numerical_test_sequence)

        predicted_labels = np.argmax(alpha_test * beta_test, axis=1)
        predicted_states = ['Helix', 'Strand', 'Turn', 'Unknown']
        predicted_labels = [predicted_states[label] for label in predicted_labels]

        return predicted_labels

# Our path
file_path = "Training_Dataset.tsv"

sequences, labels = read_tsv(file_path)

num_states = 4  # Helix, Strand, Turn, Unknown
num_symbols = 26  # 26 amino acids (A-Z)
hmm_model = HMM(num_states, num_symbols)

flat_sequences = np.array([ord(aa) - ord('A') for sequence in sequences for aa in sequence])


norm_changes = hmm_model.train(flat_sequences[:1000])

#UBE2C protein sequence
test_sequence = "MASQNRDPAATSVAAARKGAEPSGGAARGPVGKRLQQELMTLMMSGDKGISAFPESDNLF" \
                 "KWVGTIHGAAGTVYEDLRYKLSLEFPSGYPYNAPTVKFLTPCYHPNVDTQGNICLDILKE" \
                 "KWSALYDVRTILLSIQSLLGEPNIDSPLNTHAAELWKNPTAFKKYLQETYSKQVTSQEP"


predicted_labels = hmm_model.predict(test_sequence)


print("\nPredicted Labels for UBE2C_Human Protein Sequence:")
print(predicted_labels)

plot_hmm(hmm_model.transition_matrix, hmm_model.emission_matrix)
