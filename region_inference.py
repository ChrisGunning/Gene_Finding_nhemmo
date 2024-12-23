"""
Inference script. Predicts genetic regions.
"""

from dotenv import load_dotenv
import os
import argparse
import numpy as np
import csv

# Load environment variables
load_dotenv()
MAX_MARKOV_ORDER = int(os.getenv('MAX_MARKOV_ORDER'))

INTERGENIC = 0
INTRON_0 = 1
INTRON_1 = 2
INTRON_2 = 3
INITIAL_EXON = 4
FINAL_EXON = 5
INTERNAL_EXON = 6
SINGLE_EXON = 7

A = 0
C = 1
G = 2
T = 3

EMISSIONS = {'A': A, 'C': C, 'G': G, 'T': T}

INTERGENIC_ORDER = [3,2,2,2,2,1]
INTRON_0_ORDER = [1,2,3,2,2,1]
INTRON_1_ORDER = [1,1,2,3,3,2]
INTRON_2_ORDER = [1,2,3,2,2,1]
INITIAL_EXON_ORDER = [1,1,2,3,4,3]
FINAL_EXON_ORDER = [1,1,2,3,3,2]
INTERNAL_EXON_ORDER = [1,1,2,3,2,1]
SINGLE_EXON_ORDER = [1,2,2,3,4,4]

if os.path.exists("matrices/transition_matrix_1.npy"):
    markov1 = "matrices/transition_matrix_1.npy"
if os.path.exists("matrices/transition_matrix_2.npy"):
    markov2 = "matrices/transition_matrix_2.npy"
if os.path.exists("matrices/transition_matrix_3.npy"):
    markov3 = "matrices/transition_matrix_3.npy"
if os.path.exists("matrices/transition_matrix_4.npy"):
    markov4 = "matrices/transition_matrix_4.npy"
if os.path.exists("matrices/transition_matrix_5.npy"):
    markov5 = "matrices/transition_matrix_5.npy"
if os.path.exists("matrices/transition_matrix_6.npy"):
    markov6 = "matrices/transition_matrix_6.npy"
if os.path.exists("matrices/transition_matrix_7.npy"):
    markov7 = "matrices/transition_matrix_7.npy"
if os.path.exists("matrices/transition_matrix_8.npy"):
    markov8 = "matrices/transition_matrix_8.npy"
if os.path.exists("matrices/emission_matrix_0.npy"):
    emission_file = "matrices/emission_matrix_0.npy"

ORIGINAL_VALS_DICT = {}

'''Reads the fasta file and outputs the sequence to analyze.
Arguments:
	filename: name of the fasta file
Returns:
	s: string with relevant sequence
'''
def read_fasta(filename):
    with open(filename, "r") as f:
        s = ""
        for l in f.readlines()[1:]:
            line = l.strip()
            if line[0] == ">":
                break
            s += line
    return s


'''
Reads the transition and emission probabilities from .npy files.
Returns a dictionary of (Matrix Name, 2D NumPy Array Names).
'''
def read_probabilities(transition_list, emission_filename):
    dictionary = {}
    for order, filename in enumerate(transition_list):
        dictionary["Transition Order " + str(order+1)] = np.load(filename)
    dictionary["Emissions"] = np.load(emission_filename)

    return dictionary

'''Recursive helper to execute calculations based on varying markov model order.'''
def viterbi_helper_recursive(t, n, weights, state,trans_probs, obs, emiss_probs, count, calc_so_far = None):
    if n == 0:
        #base case
        if None in calc_so_far:
            print("ERROR")

        last_calc = [calc_so_far[i] + trans_probs[i][state] for i in range(len(calc_so_far))]
        result = 0
        for i in range(len(calc_so_far)):
            result += weights[i]*np.exp(last_calc[i])
        return np.log(result), []

    else:
        # stack calculation so far for each markov order
        calc_so_far_intergenic = [None for i in range(len(calc_so_far))]
        calc_so_far_intron0 = [None for i in range(len(calc_so_far))]
        calc_so_far_intron1 = [None for i in range(len(calc_so_far))]
        calc_so_far_intron2 = [None for i in range(len(calc_so_far))]
        calc_so_far_initial_exon = [None for i in range(len(calc_so_far))]
        calc_so_far_final_exon = [None for i in range(len(calc_so_far))]
        calc_so_far_internal_exon =  [None for i in range(len(calc_so_far))]
        calc_so_far_single_exon = [None for i in range(len(calc_so_far))]

        done = False
        for i in range(len(calc_so_far)):
            if calc_so_far[i] is None:
                if done:
                    break
                else:
                    done = True
                
                calc_so_far_intergenic[i] = ORIGINAL_VALS_DICT[(t, INTERGENIC)] + emiss_probs[INTERGENIC][EMISSIONS[obs[t]]]
                calc_so_far_intron0[i] = ORIGINAL_VALS_DICT[(t, INTRON_0)]+ emiss_probs[INTRON_0][EMISSIONS[obs[t]]]
                calc_so_far_intron1[i] = ORIGINAL_VALS_DICT[(t, INTRON_1)]+ emiss_probs[INTRON_1][EMISSIONS[obs[t]]]
                calc_so_far_intron2[i] = ORIGINAL_VALS_DICT[(t, INTRON_2)]+ emiss_probs[INTRON_2][EMISSIONS[obs[t]]] 
                calc_so_far_initial_exon[i] = ORIGINAL_VALS_DICT[(t, INITIAL_EXON)]+ emiss_probs[INITIAL_EXON][EMISSIONS[obs[t]]]
                calc_so_far_final_exon[i] = ORIGINAL_VALS_DICT[(t, FINAL_EXON)]+ emiss_probs[FINAL_EXON][EMISSIONS[obs[t]]] 
                calc_so_far_internal_exon[i] =  ORIGINAL_VALS_DICT[(t, INTERNAL_EXON)]+ emiss_probs[INTERNAL_EXON][EMISSIONS[obs[t]]]
                calc_so_far_single_exon[i] = ORIGINAL_VALS_DICT[(t, SINGLE_EXON)]+ emiss_probs[SINGLE_EXON][EMISSIONS[obs[t]]]
            else:
                calc_so_far_intergenic[i] = calc_so_far[i] + emiss_probs[INTERGENIC][EMISSIONS[obs[t]]]
                calc_so_far_intron0[i] = calc_so_far[i] + emiss_probs[INTRON_0][EMISSIONS[obs[t]]]
                calc_so_far_intron1[i] = calc_so_far[i]+ emiss_probs[INTRON_1][EMISSIONS[obs[t]]]
                calc_so_far_intron2[i] = calc_so_far[i]+ emiss_probs[INTRON_2][EMISSIONS[obs[t]]] 
                calc_so_far_initial_exon[i] = calc_so_far[i]+ emiss_probs[INITIAL_EXON][EMISSIONS[obs[t]]]
                calc_so_far_final_exon[i] = calc_so_far[i]+ emiss_probs[FINAL_EXON][EMISSIONS[obs[t]]] 
                calc_so_far_internal_exon[i] =  calc_so_far[i]+ emiss_probs[INTERNAL_EXON][EMISSIONS[obs[t]]]
                calc_so_far_single_exon[i] = calc_so_far[i]+ emiss_probs[SINGLE_EXON][EMISSIONS[obs[t]]]

        #transition prob so far for each markov order
        trans_probs_intergenic = [trans_probs[i] for i in range(len(trans_probs))]
        trans_probs_intron0 = [trans_probs[i] for i in range(len(trans_probs))]
        trans_probs_intron1 = [trans_probs[i] for i in range(len(trans_probs))]
        trans_probs_intron2 = [trans_probs[i] for i in range(len(trans_probs))]
        trans_probs_initial_exon = [trans_probs[i] for i in range(len(trans_probs))]
        trans_probs_final_exon = [trans_probs[i] for i in range(len(trans_probs))]
        trans_probs_internal_exon =  [trans_probs[i] for i in range(len(trans_probs))]
        trans_probs_single_exon = [trans_probs[i] for i in range(len(trans_probs))]

        for i in range(count):
            trans_probs_intergenic[i] = trans_probs[i][INTERGENIC]
            trans_probs_intron0[i] = trans_probs[i][INTRON_0]
            trans_probs_intron1[i] = trans_probs[i][INTRON_1]
            trans_probs_intron2[i] = trans_probs[i][INTRON_2]
            trans_probs_initial_exon[i] = trans_probs[i][INITIAL_EXON]
            trans_probs_final_exon[i] = trans_probs[i][FINAL_EXON]
            trans_probs_internal_exon [i]=  trans_probs[i][INTERNAL_EXON]
            trans_probs_single_exon[i] = trans_probs[i][SINGLE_EXON]


        #recursive calls

        intergenic_score, intergenic_states = viterbi_helper_recursive(
                t+1, n-1, weights, state, trans_probs_intergenic, obs, emiss_probs, 
                count+1, calc_so_far=calc_so_far_intergenic)
        if intergenic_states is None:
            intergenic_states = []      
        intergenic_states.insert(0,INTERGENIC)

        intron0_score, intron0_states = viterbi_helper_recursive(
                t+1, n-1, weights, state, trans_probs_intron0, obs, emiss_probs, 
                count+1, calc_so_far=calc_so_far_intron0)
        if intron0_states is None:
            intron0_states = []      
        intron0_states.insert(0,INTRON_0)

        intron1_score, intron1_states = viterbi_helper_recursive(
                t+1, n-1, weights, state, trans_probs_intron1, obs, emiss_probs, 
                count+1, calc_so_far=calc_so_far_intron1)
        if intron1_states is None:
            intron1_states = []      
        intron1_states.insert(0,INTRON_1)

        intron2_score, intron2_states = viterbi_helper_recursive(
            t+1, n-1, weights, state, trans_probs_intron2, obs, emiss_probs, 
            count+1, calc_so_far=calc_so_far_intron2)
        if intron2_states is None:
            intron2_states = []  
        intron2_states.insert(0,INTRON_2)

        initial_exon_score, initial_exon_states = viterbi_helper_recursive(
            t+1, n-1, weights, state, trans_probs_initial_exon, obs, emiss_probs, 
            count+1, calc_so_far=calc_so_far_initial_exon)
        if initial_exon_states is None:
            initial_exon_states = []  
        initial_exon_states.insert(0,INITIAL_EXON)
        
        final_exon_score, final_exon_states = viterbi_helper_recursive(
            t+ 1, n-1, weights, state, trans_probs_final_exon, obs, emiss_probs, 
            count+1, calc_so_far=calc_so_far_final_exon)
        if final_exon_states is None:
            final_exon_states = []  
        final_exon_states.insert(0,FINAL_EXON)

        internal_exon_score, internal_exon_states = viterbi_helper_recursive(
            t+  1, n-1, weights, state, trans_probs_internal_exon, obs, emiss_probs, 
            count+1, calc_so_far=calc_so_far_internal_exon)
        if internal_exon_states is None:
            internal_exon_states = []  
        internal_exon_states.insert(0,INTERNAL_EXON)

        single_exon_score, single_exon_states = viterbi_helper_recursive(
                t+1, n-1, weights, state,trans_probs_single_exon, obs, emiss_probs, 
                count+1, calc_so_far=calc_so_far_single_exon)
        if single_exon_states is None:
            single_exon_states = [] 
        single_exon_states.insert(0,SINGLE_EXON) 

        compare_scores = [(intergenic_score, intergenic_states),
                        (intron0_score, intron0_states),
                        (intron1_score, intron1_states),
                        (intron2_score, intron2_states),
                        (initial_exon_score, initial_exon_states),
                        (final_exon_score, final_exon_states),
                        (internal_exon_score, internal_exon_states),
                        (single_exon_score, single_exon_states)]
        
        return max(compare_scores)


'''
Returns the max value and max state given an input state at time t.
'''
def viterbi_helper(t, state, emiss_probs, trans_probs_dict, weights, obs):
    markov_order = MAX_MARKOV_ORDER
    if markov_order > t:
        markov_order = t

    trans_probs = [trans_probs_dict["Transition Order " + str(i)] for i in range(markov_order,0,-1)]
    calc_so_far = [None for i in range(markov_order)]
    
    result = viterbi_helper_recursive(t-markov_order+1, markov_order, weights, 
                                      state, trans_probs, obs, emiss_probs,1, 
                                      calc_so_far)
    return result


''' Outputs the Viterbi decoding of a given observation. Using the interpolated
    Markov Models up to MAX_MARKOV_ORDER.
'''
def viterbi(obs, dictionary, init_probs_other, init_probs_start):
    #adjust markov order weights based on max markov order
    INTERGENIC_WEIGHTS = np.array(INTERGENIC_ORDER[:MAX_MARKOV_ORDER])/sum(INTERGENIC_ORDER[:MAX_MARKOV_ORDER])
    INTRON_0_WEIGHTS = np.array(INTRON_0_ORDER[:MAX_MARKOV_ORDER])/sum(INTRON_0_ORDER[:MAX_MARKOV_ORDER])
    INTRON_1_WEIGHTS = np.array(INTRON_1_ORDER[:MAX_MARKOV_ORDER])/sum(INTRON_1_ORDER[:MAX_MARKOV_ORDER])
    INTRON_2_WEIGHTS = np.array(INTRON_2_ORDER[:MAX_MARKOV_ORDER])/sum(INTRON_2_ORDER[:MAX_MARKOV_ORDER])
    INITIAL_EXON_WEIGHTS = np.array(INITIAL_EXON_ORDER[:MAX_MARKOV_ORDER])/sum(INITIAL_EXON_ORDER[:MAX_MARKOV_ORDER])
    FINAL_EXON_WEIGHTS = np.array(FINAL_EXON_ORDER[:MAX_MARKOV_ORDER])/sum(FINAL_EXON_ORDER[:MAX_MARKOV_ORDER])
    INTERNAL_EXON_WEIGHTS = np.array(INTERNAL_EXON_ORDER[:MAX_MARKOV_ORDER])/sum(INTERNAL_EXON_ORDER[:MAX_MARKOV_ORDER])
    SINGLE_EXON_WEIGHTS = np.array(SINGLE_EXON_ORDER[:MAX_MARKOV_ORDER])/sum(SINGLE_EXON_ORDER[:MAX_MARKOV_ORDER])

    emiss_probs = dictionary["Emissions"]
    trans_probs_dict = dictionary

    # Backtracing
    L = len(obs)
    intergenic_backtracing = [None for i in range(L-1)]
    intron_phase0_backtracing = [None for i in range(L-1)]
    intron_phase1_backtracing = [None for i in range(L-1)]
    intron_phase2_backtracing = [None for i in range(L-1)]
    exon_initial_backtracing = [None for i in range(L-1)]
    exon_final_backtracing = [None for i in range(L-1)]
    exon_inter_backtracing = [None for i in range(L-1)]
    exon_single_backtracing = [None for i in range(L-1)]
    backtracing_lst = [intergenic_backtracing, intron_phase0_backtracing, intron_phase1_backtracing, intron_phase2_backtracing,
                       exon_initial_backtracing, exon_final_backtracing, exon_inter_backtracing, exon_single_backtracing]

    # Initialization
    emission = obs[0]

    if len(obs) >= 3 and obs[0:3] == "ATG":
        init_probs = init_probs_start
    else:
        init_probs = init_probs_other

    v_intergenic_1 = init_probs[INTERGENIC] + emiss_probs[INTERGENIC][EMISSIONS[emission]]
    v_intron0_1 = init_probs[INTRON_0] + emiss_probs[INTRON_0][EMISSIONS[emission]]
    v_intron1_1 = init_probs[INTRON_1] + emiss_probs[INTRON_1][EMISSIONS[emission]]
    v_intron2_1 = init_probs[INTRON_2] + emiss_probs[INTRON_2][EMISSIONS[emission]]
    v_exon_initial_1 = init_probs[INITIAL_EXON] + emiss_probs[INITIAL_EXON][EMISSIONS[emission]]
    v_exon_final_1 = init_probs[FINAL_EXON] + emiss_probs[FINAL_EXON][EMISSIONS[emission]]
    v_exon_inter_1 = init_probs[INTERNAL_EXON] + emiss_probs[INTERNAL_EXON][EMISSIONS[emission]]
    v_exon_single_1 = init_probs[SINGLE_EXON] + emiss_probs[SINGLE_EXON][EMISSIONS[emission]]

    # Iteration
    v_intergenic_t = v_intergenic_1
    v_intron0_t = v_intron0_1
    v_intron1_t = v_intron1_1
    v_intron2_t = v_intron2_1
    v_exon_initial_t = v_exon_initial_1
    v_exon_final_t = v_exon_final_1
    v_exon_inter_t = v_exon_inter_1
    v_exon_single_t = v_exon_single_1

    for t in range(1, L):
        print(str(t) + " / " + str(L))
        ORIGINAL_VALS_DICT[(t, INTERGENIC)] = v_intergenic_t
        ORIGINAL_VALS_DICT[(t, INTRON_0)] = v_intron0_t
        ORIGINAL_VALS_DICT[(t, INTRON_1)] = v_intron1_t
        ORIGINAL_VALS_DICT[(t, INTRON_2)] = v_intron2_t
        ORIGINAL_VALS_DICT[(t, INITIAL_EXON)] = v_exon_initial_t
        ORIGINAL_VALS_DICT[(t, FINAL_EXON)] = v_exon_final_t
        ORIGINAL_VALS_DICT[(t, INTERNAL_EXON)] = v_exon_inter_t
        ORIGINAL_VALS_DICT[(t, SINGLE_EXON)] = v_exon_single_t

        # original_vals = [original_v_intergenic_t, original_v_intron0_t, original_v_intron1_t, original_v_intron2_t, original_v_exon_initial_t, original_v_exon_final_t, original_v_exon_inter_t, original_v_exon_single_t]

        emission = obs[t]

        max_intergenic_val, max_intergenic_state = viterbi_helper(t, INTERGENIC, emiss_probs, trans_probs_dict, INTERGENIC_WEIGHTS, obs)
        v_intergenic_t = max_intergenic_val
        backtracing_lst[INTERGENIC][t-1] = max_intergenic_state

        max_intron0_val, max_intron0_state = viterbi_helper(t, INTRON_0, emiss_probs, trans_probs_dict, INTRON_0_WEIGHTS, obs)
        v_intron0_t = max_intron0_val
        backtracing_lst[INTRON_0][t-1] = max_intron0_state

        max_intron1_val, max_intron1_state = viterbi_helper(t, INTRON_1, emiss_probs, trans_probs_dict, INTRON_1_WEIGHTS, obs)
        v_intron1_t = max_intron1_val
        backtracing_lst[INTRON_1][t-1] = max_intron1_state

        max_intron2_val, max_intron2_state = viterbi_helper(t, INTRON_2, emiss_probs, trans_probs_dict, INTRON_2_WEIGHTS, obs)
        v_intron2_t = max_intron2_val
        backtracing_lst[INTRON_2][t-1] = max_intron2_state

        max_exon_initial_val, max_exon_initial_state = viterbi_helper(t, INITIAL_EXON, emiss_probs, trans_probs_dict, INITIAL_EXON_WEIGHTS, obs)
        v_exon_initial_t = max_exon_initial_val
        backtracing_lst[INITIAL_EXON][t-1] = max_exon_initial_state

        max_exon_final_val, max_exon_final_state = viterbi_helper(t, FINAL_EXON, emiss_probs, trans_probs_dict, FINAL_EXON_WEIGHTS, obs)
        v_exon_final_t = max_exon_final_val
        backtracing_lst[FINAL_EXON][t-1] = max_exon_final_state

        max_exon_inter_val, max_exon_inter_state = viterbi_helper(t, INTERNAL_EXON, emiss_probs, trans_probs_dict, INTERNAL_EXON_WEIGHTS, obs)
        v_exon_inter_t = max_exon_inter_val
        backtracing_lst[INTERNAL_EXON][t-1] = max_exon_inter_state

        max_exon_single_val, max_exon_single_state = viterbi_helper(t, SINGLE_EXON, emiss_probs, trans_probs_dict, SINGLE_EXON_WEIGHTS, obs)
        v_exon_single_t = max_exon_single_val
        backtracing_lst[SINGLE_EXON][t-1] = max_exon_single_state


    # Termination
    l = []

    compare_scores = [(v_intergenic_t, INTERGENIC),
                          (v_intron0_t, INTRON_0),
                          (v_intron1_t, INTRON_1),
                          (v_intron2_t, INTRON_2),
                          (v_exon_initial_t, INITIAL_EXON),
                          (v_exon_final_t, FINAL_EXON),
                          (v_exon_inter_t, INTERNAL_EXON),
                          (v_exon_single_t, SINGLE_EXON)]
    
    
    max_val, max_state = max(compare_scores)
    p = max_val
    l.insert(0,max_state)
    
    i = L-2
    while i >= 0:
        insert_lst = backtracing_lst[l[0]][L-2-i]
        for j in range(len(insert_lst)-1,-1,-1):
            i -= 1
            l.insert(0,insert_lst[j])

    return l, p


def main():
    transition_list = [markov1, markov2, markov3, markov4, markov5, markov6]

    fasta_file = os.getenv('FASTA_FILEPATH')
    out_file = os.getenv('INFERENCE_FILEPATH')

    dictionary = read_probabilities(transition_list, emission_file)
    obs_sequence = read_fasta(fasta_file)
    print("Read probabilities and FASTA.")
    

    initial_probabilies = np.log(np.array([0.9,0.1/8,0.1/8,0.1/8,0.1/8,0.1/8,0.1/8,0.1/8]))
    initial_probabilies_start = np.log(np.array([0.75,0.05/5,0.05/5,0.05/5,0.1,0.05/5,0.05/5,0.1]))

    sequence, p = viterbi(obs_sequence, dictionary, initial_probabilies, initial_probabilies_start)

    with open(out_file, 'w') as file:
        writer = csv.writer(file)
        writer.writerow(sequence)  # Write the list as a single row

    print("Viterbi probability in log scale: {:.2f}".format(p))


if __name__ == "__main__":
    main()
    print("done")
