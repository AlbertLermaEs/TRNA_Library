#TRNA Library  
from itertools import combinations_with_replacement
from itertools import permutations
from Bio import pairwise2
import re
import pandas as pd
from nupack import *
# Module that generates combinations
a = ['A', 'U', 'G', 'C']
comb = list(combinations_with_replacement(a, 5))
# Module that generates permutations based on the previously generated combinations.
# Permutations are sent to a set to have exclusively unique elements.
perm = list(permutations(comb[0])) + list(permutations(comb[1])) + list(permutations(comb[2])) + list(
    permutations(comb[3])) + list(permutations(comb[4])) + list(permutations(comb[5])) + list(
    permutations(comb[6])) + list(permutations(comb[7])) + list(permutations(comb[8])) + list(
    permutations(comb[9])) + list(permutations(comb[10])) + list(permutations(comb[11])) + list(
    permutations(comb[12])) + list(permutations(comb[13])) + list(permutations(comb[14])) + list(
    permutations(comb[15])) + list(permutations(comb[16])) + list(permutations(comb[17])) + list(
    permutations(comb[18])) + list(permutations(comb[19])) + list(permutations(comb[20])) + list(
    permutations(comb[21])) + list(permutations(comb[22])) + list(permutations(comb[23])) + list(
    permutations(comb[24])) + list(permutations(comb[25])) + list(permutations(comb[26])) + list(
    permutations(comb[27])) + list(permutations(comb[28])) + list(permutations(comb[29])) + list(
    permutations(comb[30])) + list(permutations(comb[31])) + list(permutations(comb[32])) + list(
    permutations(comb[33])) + list(permutations(comb[34])) + list(permutations(comb[35])) + list(
    permutations(comb[36])) + list(permutations(comb[37])) + list(permutations(comb[38])) + list(
    permutations(comb[39])) + list(permutations(comb[40])) + list(permutations(comb[41])) + list(
    permutations(comb[42])) + list(permutations(comb[43])) + list(permutations(comb[44])) + list(
    permutations(comb[45])) + list(permutations(comb[46])) + list(permutations(comb[47])) + list(
    permutations(comb[48])) + list(permutations(comb[49])) + list(permutations(comb[50])) + list(
    permutations(comb[51])) + list(permutations(comb[52])) + list(permutations(comb[53])) + list(
    permutations(comb[54])) + list(permutations(comb[55]))
uniques = set(perm)
# The set is converted to a list
listofuniques = list(uniques)
# All the sequences of each tuple in the list are concatenated, resulting in a list of sequences.
seqlist = [''.join(tups) for tups in listofuniques]
# The previous modified code is used to obtain the reverses of each sequence.
revseq = [''.join(tups[::-1]) for tups in listofuniques]
# The Shine Dalgarno sequence, a loop of 10 A, and the anti SD sequence are defined.
SD = 'AGGAGG'
loop = 'AAAAAAAAAA'
antiSD = 'UCCU'
# Finally, the sequences from the newseq and revseq lists are concatenated, resulting in a list of 1048576 thermoregulators.
newseq = [x + antiSD + loop + SD for x in seqlist]
TRNAlibrary = []
for x in (newseq):
    for y in (revseq):
        TRNAlibrary.append(x + y)

# Define functions
# Function to calculate loop size
def calculate_loop_size(structure):
    # Find all occurrences of dots between parentheses using a regular expression
    matches = re.findall(r'\(\.+?\)', structure)

    # Get the loop size by summing the length of all found occurrences
    loop_size = sum(len(match) - 2 for match in matches)

    # Return the loop size
    return loop_size

# Function to calculate stem size
def calculate_stem_size(structure):
    stem_size = 0
    counter = 0

    for character in structure:
        if character == '(':
            counter += 1
        elif character == ')':
            if counter > 0:
                stem_size += counter
                counter = 0

    return stem_size

# Function GC percentage in stem
def calculate_GC_percentage_stem(sequence, structure):
    stem_sequence = ''.join([base for base, structure_base in zip(sequence, structure) if structure_base in ('(', ')')])

    total_bases = len(stem_sequence)
    gc_bases = stem_sequence.count('G') + stem_sequence.count('C')
    percentage_GC = (gc_bases / total_bases) * 100
    percentage_GC = round(percentage_GC, 2)  # Round to two decimal places

    return percentage_GC


# Internal loops function
def calculate_internal_loops(structure):
    # Find the index of the first "(" and the last "("
    first_open_parenthesis = structure.find("(")
    last_open_parenthesis = structure.rfind("(")

    # Find the index of the first ")" and the last ")"
    first_closed_parenthesis = structure.find(")")
    last_closed_parenthesis = structure.rfind(")")

    # Get the left stem and the right stem
    left_stem = structure[first_open_parenthesis:last_open_parenthesis + 1]
    right_stem = structure[first_closed_parenthesis:last_closed_parenthesis + 1]

    # Convert "(" and ")" to 1 and "." to 0 in the stems
    left_stem_converted = left_stem.replace("(", "1").replace(")", "1").replace(".", "0")
    right_stem_converted = right_stem.replace("(", "1").replace(")", "1").replace(".", "0")

    seq1 = left_stem_converted
    seq2 = right_stem_converted[::-1]

    # Perform global alignment using the Needleman-Wunsch algorithm
    alignments = pairwise2.align.globalxx(seq1, seq2)

    # Get the best alignment
    best_alignment = alignments[0]
    aligned_seq1 = best_alignment.seqA
    aligned_seq2 = best_alignment.seqB

    consensus_sequence = ""

    for i in range(len(aligned_seq1)):
        if aligned_seq1[i] == aligned_seq2[i]:
            consensus_sequence += aligned_seq1[i]
        else:
            consensus_sequence += "0"

    interruptions = 0
    cont_ones = 0

    for i in range(len(consensus_sequence)):
        if consensus_sequence[i] == "1":
            cont_ones += 1
        elif consensus_sequence[i] == "0" and cont_ones > 0:
            interruptions += 1
            cont_ones = 0

    return interruptions

#Dataframe generation
temperatures = [20, 25, 30, 35, 40]
sequences = TRNAlibrary 

results = []
global_sequence_counter = 1

for sequence in sequences:
    strands = sequence
    temp_results = []

    for temperature in temperatures:
        my_model = Model(material='RNA', celsius=temperature)
        mfe_structures = mfe(strands, model=my_model)
        rounded_energy = round(mfe_structures[0].energy, 2)

        # Get additional characteristics
        structure = str(mfe_structures[0].structure)
        loop_size = calculate_loop_size(structure)
        stem_size = calculate_stem_size(structure)
        GC_percentage_stem = calculate_GC_percentage_stem(sequence, structure)
        internal_loops = calculate_internal_loops(structure)

        # Generate the sequence name
        trna_name = f"TRNA{global_sequence_counter}"

        temp_results.append({
            "Name": trna_name,
            "Sequence": sequence,
            "Temperature": temperature,
            "Energy": rounded_energy,
            "Structure": structure,
            "Loop_Size": loop_size,
            "Stem_Size": stem_size,
            "GC_Stem": GC_percentage_stem,
            "Internal_Loops": internal_loops
        })

    # Increment the global unique sequence counter
    global_sequence_counter += 1

    # Add temp_results to the main results list
    results.extend(temp_results)

# Create a DataFrame from the results list
df = pd.DataFrame(results)

# Save the DataFrame to a CSV file
df.to_csv("TRNA_library_LermaEscalera.csv", index=False)