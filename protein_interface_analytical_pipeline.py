# Import libraries
import pandas as pd
from Bio import PDB
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB import PDBIO, Select


# Define functions 

# LOADDING STRUCTURES FROM PDB FILES
def load_structure(pdb_id, pdb_file):
    """
    Load a PDB structure using a PDB ID and file path.
    """
    parser = PDBParser()
    structure = parser.get_structure(pdb_id, pdb_file)
    return structure


# IDENTIFYING OVERLAPPING AAs BETWEEN TWO CHAINS/PROTEINS
def identify_interface_residues(structure, prot_a_chain_id, distance_threshold=3.3): # 4.0 accounts for vann der vaals, 3.2 only for hydrogen bonds
    """
    Identify interface residues between Prot_A and other proteins in the structure.
    Returns a set of tuples where each tuple contains three elements:
    (residue on Prot_A, residue on interacting chain, chain ID of interacting residue)
    """
    prot_a_chain = structure[0][prot_a_chain_id]
    atoms_prot_a = list(prot_a_chain.get_atoms())

    ns = NeighborSearch(atoms_prot_a)
    interface_residues = set()
    
    for chain in structure.get_chains():
        if chain.id == prot_a_chain_id:
            continue  # Skip Prot_A's chain

        for residue in chain:
            for atom in residue:
                close_residues = ns.search(atom.coord, distance_threshold, level='R')
                
                for prot_a_residue in close_residues:
                    interface_residues.add((prot_a_residue, residue, chain.id))

    return interface_residues


# KIND-TO-THE-EYE PRINTING OF INTERFACE RESIDUES
def print_interface_residues(interface_residues):
    for prot_a_residue, interacting_residue, interacting_chain_id in interface_residues:
        print(f"Prot_A Residue: {prot_a_residue.get_resname()} {prot_a_residue.get_id()[1]} (Chain {prot_a_residue.get_parent().id}), \
                Interacting Residue: {interacting_residue.get_resname()} {interacting_residue.get_id()[1]} (Chain {interacting_chain_id})")

    print("Unique residues in Prot_A chain:")

    # Store and count the unique residues in chain E and print them 
    unique_res_count = {}  # Dictionary to store and count unique residues
    for prot_a_residue, _, _ in interface_residues:
        residue_key = (prot_a_residue.get_resname(), prot_a_residue.get_id()[1], prot_a_residue.get_parent().id)
        if residue_key in unique_res_count:
            unique_res_count[residue_key] += 1
        else:
            unique_res_count[residue_key] = 1
    
    # Sort and print them with counts, based on residue number
    for res in sorted(unique_res_count.keys(), key=lambda x: x[1]):
        count = unique_res_count[res]
        print(f"Chain Prot_A ({res[2]}) residue: {res[0]} {res[1]}, Count: {count}")


# EXTRACT NUMERICAL POSITIONS OF THE INTERFACE RESIDUES, TO LATER FIND OVERLAPPING REGION IN A+C AND A+B
def chain_interface_residue_positions(interface_res, chain): 

    unique_res_pos = set()

    for prot_a_residue, interacting_residue, interacting_chain_id in interface_res:
        # print(prot_a_residue.get_parent().id)
        if prot_a_residue.get_parent().id == chain: # only store the ones from the chain we are interested in, spike is chain E in AC and chain R in AB
            unique_res_pos.add(prot_a_residue.get_id()[1])


    return unique_res_pos


# SUBSTRACT THE RESIDUES THAT ARE ALSO IN THE INTERFACE WITH A+C, TO NOT MESS WITH THAT BINDING 
def filter_out_residue_range(interface_residues, chain_id, start_position, end_position):
    """
    Filter out residues within a specified range from a set of interface residues.

    Parameters:
    interface_residues (set): A set of tuples, each containing two Residue objects and an additional piece of information.
    chain_id (str): The chain ID from which to exclude residues (e.g., 'R').
    start_position (int): The start position of the range to exclude.
    end_position (int): The end position of the range to exclude.

    Returns:
    filtered_residues (set): A set of tuples with the filtered residue pairs and the additional information.
    """
    filtered_residues = set()

    for residue_pair in interface_residues:
        residue1, residue2, additional_info = residue_pair
        filtered_pair = []

        for residue in [residue1, residue2]:
            residue_position = residue.get_id()[1]
            # Check if the residue is in the specified chain and outside the exclude range
            if residue.get_parent().id != chain_id or (residue_position < start_position or residue_position > end_position):
                filtered_pair.append(residue)

        # Only add pairs that still have both residues and the additional information
        if len(filtered_pair) == 2:
            filtered_residues.add((filtered_pair[0], filtered_pair[1], additional_info))

    return filtered_residues


# Extract sequence of the protein from the PDB structure
def get_sequence_from_pdb(structure, chain_id):
    # get the aa sequence from the structure 
     for model in structure:
        for chain in model:
            if chain.id == chain_id:
                return "".join([PDB.Polypeptide.three_to_one(residue.get_resname()) for residue in chain if PDB.Polypeptide.is_aa(residue, standard=True)])


# Perform point mutation PHE486>ALA486 (F486A) and provide the sequence
def mutate_sequence(sequence, mutations):
    """
    Mutate specific residues in a protein sequence.

    Parameters:
    sequence (str): Original amino acid sequence.
    mutations (dict): Dictionary where keys are residue positions (1-indexed) and values are the new amino acids.

    Returns:
    str: Mutated sequence.
    """
    sequence_list = list(sequence)
    for position, new_residue in mutations.items():
        sequence_list[position] = new_residue  # Converting 1-indexed position to 0-indexed not needed 
    return "".join(sequence_list)




# PIPELINE EXECUTION 



# LOAD STRUCTURES
strAC = load_structure("6M0J", "6M0J.pdb")  # spike (Prot_A) bound to receptor (C)
strAB = load_structure("7Z0X", "7Z0X.pdb")  # spike (Prot_A) bound to antibody (B)

# FIND AND PRINT INTERFACE RESIDUES 
interface_residues_strAC = identify_interface_residues(strAC, 'E')  # Replace 'A' with the chain ID of Prot_A in strAC
interface_residues_strAB = identify_interface_residues(strAB, 'R')  # Replace 'A' with the chain ID of Prot_A in strAB

# EXTRACT POSITION OF INTERFACE RESIDUES 
int_pos_AC = chain_interface_residue_positions(interface_residues_strAC, 'E')
int_pos_AB = chain_interface_residue_positions(interface_residues_strAB, 'R')

# Obtain the position of the interface amino acids common in the interaction in AB and AC
AB_AC_pos = set() # store in this list

for pos in sorted(int_pos_AB):
    if pos in int_pos_AC: # if the residue is in both lists, store it 
        print(pos)
        AB_AC_pos.add(pos) 


# Substract the region with amino acids that participate in AC binding from the R chain (in the AB complex)
# This will ensure that we can choose a mutation in a region that doesn't affect AC binding
f_interface_residues_AB = filter_out_residue_range(interface_residues_strAB, 'R', min(AB_AC_pos), max(AB_AC_pos))

print_interface_residues(f_interface_residues_AB)



# The main candidates for mutation are: 
# 4 H-bonds: SER 477, THR 478, VAL 483, GLU 484, TYR 489
# 5 H-bonds: GLY 485, HOH 622 (water)
# 11 H-bonds: PHE 486

# Strongest candidate: PHE486
# PHE486 is not an interface residue in Prot_A+Prot_C interaction
# PHE usually interacts with a hydrophobic pocket in the binding protein

# PHE486>ALA486 (also Phe486Ala or F486A) point mutation is a candidate: 
# can cut interactions with the bydrophobic pocket due to reduced size and still mantain non-polarity of that position



# Extract sequence of the protein from the PDB structure (specific chain)
sequence = get_sequence_from_pdb(strAB, 'R')  # Replace 'A' with the chain ID you're interested in
print(f'Original sequence chain R : \n', sequence)

# Perform point mutation PHE486>ALA486 (F486A) and provide the sequence

# Position annotation of chain R starts at position 333, so a conversion is necessary
def convert_pos(annotated_pos, starting_aa_num):
    pos_in_seq = annotated_pos - starting_aa_num # -1 becasue 333 is the first aa (e.g. 10th position will be 342, not 343)
    return pos_in_seq

cpos = convert_pos(486, 333)

mutations = {cpos: 'A'}  # mutate 486th residue to A, can also take more inputs
mutated_sequence = mutate_sequence(sequence, mutations)
print('Mutated sequence: \n', mutated_sequence)

print(sequence[cpos])
print(mutated_sequence[cpos])

# Now that the mutated sequence is obtained, the next steps are:
# Using AlphaFold for 3D protein prediction
# Molecular docking and calculation of change in binding affinity



# Extract sequence of the protein from the PDB structure (specific chain)
sequence = get_sequence_from_pdb(strAB, 'H')  # Replace 'H' with the chain ID you're interested in
print(f'Original sequence chain H : \n', sequence)

# Extract sequence of the protein from the PDB structure (specific chain)
sequence = get_sequence_from_pdb(strAB, 'L')  # Replace 'L' with the chain ID you're interested in
print(f'Original sequence chain L : \n', sequence)
