# Hidden-markov-model-for-protein-structure-prediction
Hidden Markov Model Based Secondary Structure Prediction of Proteins
# Background Information
Basically, a hidden Markov model (HMM) is one in which you observe a sequence of emissions, but do not know the sequence of states the model went through to generate the emissions. We analyze the hidden Markov models to recover the sequence of states from the observed data.
![image](https://github.com/user-attachments/assets/ee9e145a-a008-4b32-ac78-48f64c6d6017)
# Steps
Developing an HMM based predictor to predict the secondary structural regions of proteins as alpha helix (H), beta strand (E), turn (T), and unknown (U) and apply it on the amino acid sequence of the UBE2C protein. Development of an SS predictor includes:
i) Construction of a predictive model and training the model with labeled reference data,
ii) Calculating its prediction performance on labeled test data (i.e., UBE2C_Human protein)
iii) Comparing its performance with a baseline method (i.e., Chou-Fasman) to observe if your approach adds value to SSE prediction.
I created my own HMM algorithm here to predict emission and transition values. I used the training dataset here but it’s incredibly long and my editor and computer can not always enough to give healthy outputs so I use  one thousands of the data. Especially, I’ve just used a little bit smaller data to train here to process my algorithm clearly. Also I researched some libraries like hmmlearn in python and look at the function’s features to get to know and understand more.

Constructing the predictive model using an HMM with 4 states: (1) helix, (2) strand, (3) turn, (4) the unknown state (+ the start & end states). Calculate the transition and emission probabilities (add pseudo-counts of adding 1 to numerator and 20 to denominator for emission) using the known SS information in the given training dataset (i.e., "BBM411_Assignment2_ Q3_TrainingDataset.txt").  Using the necessary algorithm to analyze the input sequence and predict the most probable path that will emit that sequence (in terms of SSE states).
Format of the tab-delimited training dataset (columns; 1:unirprot id, 2:protein name, 3:sequence, 4:helix, 5: strand, 6:turn): O00244\tATOX1_HUMAN\tMPKHEFSVDMTCGGCAEAVSRVLNKLGGVKYDIDLPNKKVCIESEH\tHELIX 13..26; HELIX 48..56;\tSTRAND 3..8; STRAND 28..34; STRAND 39..46;STRAND 62..66;\tTURN 35..38; TURN 57..59;
In the training dataset file, there is one protein per row (including its sequence), and their residues/amino acids are assigned into three states (“Helix”, “Beta strand”, and “Turn”) with the notation “HELIX 213..229” which means amino acids from 213 to 229 (including the boundaries) belong to the helix class. There are multiple regional assignments for each class for most of the proteins (e.g., for protein X, 3 different regions
are assigned to Helix, 4 different regions are assigned to beta strand, etc.). Multiple assignments of the same class are separated from each other by the “;” character.
Residues that were left out of the 3 class assignments in the file should be considered to belong to the “Unknown” class. Use these SSE class assignments to calculate the probability values of the model.
