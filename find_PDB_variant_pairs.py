import sys
import os
import shutil
import re
import gzip
import pickle
import datetime
import time
import math
import contextlib
import multiprocessing
import argparse
# Updated biopython on Computerome:
# if os.uname()[1] == 'computerome02':
#     sys.path.insert(0, '/home/projects/cu_10020/apps/python3-site-packages/lib/python/')
sys.path.insert(0, '/home/projects/cu_10020/apps/python3-site-packages/lib/python/')
from Bio.PDB import *

# Build commandline parser:
parser = argparse.ArgumentParser(description="Mine PDB for variant chain pairs.")

# Arguments:
parser.add_argument(
    "-cache_dir",
    "--cache_directory",
    type=str,
    dest="cache_dir",
    metavar="FILE",
    help="Directory where various cached results will be saved/loaded as pickled dictionaries.",
    required=True)
parser.add_argument(
    "-ss_dis",
    "--ss_dis_file",
    type=str,
    dest="ss_dis_file",
    metavar="FILE",
    help="Sequence and disorder file from RCSB.org over all protein chains in the PDB.",
    required=True)
parser.add_argument(
    "-pdb",
    "--pdb_folder",
    type=str,
    dest="pdb_folder",
    metavar="DIR",
    help="The innermost folder that contains all the full PDB database. This folder is often called 'split' and has all the two character folders that devides the full PDB into an easily accessible structure.",
    required=True)
parser.add_argument(
    "-biolip",
    "--biolip_file",
    type=str,
    dest="biolip_fnam",
    metavar="FILE",
    help="BioLiP file containing all biological relevant ligands for all PDB chains. Found here: http://zhanglab.ccmb.med.umich.edu/BioLiP/",
    required=True)
parser.add_argument(
    "-scratch",
    "--scratch_dir",
    type=str,
    dest="scratch_dir",
    metavar="DIR",
    help="Directory for writing temporary files/folders.",
    required=True)
parser.add_argument(
    "-out",
    "--output_file",
    type=str,
    dest="result_file",
    metavar="FILE",
    help="Dataset of tab separated values with relevant features from the mined PDB chain pairs.",
    required=True)
parser.add_argument(
    "-np",
    "--n_processes",
    type=int,
    dest="np",
    metavar="int",
    help="Number of processes to use in a pool.")

# Set default arguments:
parser.set_defaults(
    cache_dir='/home/projects/cu_10020/kortemme_visit/pdb_mining',
    ss_dis_file='ss_dis.txt',
    pdb_folder='/home/projects/pr_46690/data/db/pdb/split',
    biolip_fnam='BioLiP_all.csv',
    scratch_dir='/home/projects/cu_10020/kortemme_visit/pdb_mining/scratch',
    result_file='phi_psi_vs_dist.tab',
    np=2)
args = parser.parse_args()

# System specific variables:
# cache_dir = '/Users/krdav/Dropbox/3Dmutation/PDB_single_mutation_effect'
# cache_dir = '/home/projects/cu_10020/kortemme_visit/pdb_mining'

# pdb_folder = '/Users/krdav/Dropbox/3Dmutation/PDB_single_mutation_effect/split'
# pdb_folder = '/home/projects/pr_46690/data/db/pdb/split'

# scratch_dir = '/Users/krdav/Dropbox/3Dmutation/PDB_single_mutation_effect/scratch'
# scratch_dir = '/home/projects/cu_10020/kortemme_visit/pdb_mining/scratch'

# Remove trialing slashes from directory names:
args.pdb_folder = args.pdb_folder.rstrip('/')
args.scratch_dir = args.scratch_dir.rstrip('/')

run_dir = os.getcwd()
if '/' not in args.result_file:
    args.result_file = run_dir + args.result_file

# Global variables:
residue_type_3to1_map = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "MSE": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
    "UNK": 'X',
}

AAletters = 'ACDEFGHIKLMNPQRSTVWY'

D_isomer_AA = ["DAL", "DCY", "DAP", "DGU", "DPH",
               "DGL", "DHI", "DIL", "DLY", "DLE",
               "DME", "DMS", "DAN", "DPR", "DGN",
               "DAR", "DSE", "DTH", "DVA", "DTR",
               "DTY"]

# Define HETATM of modified/non-canonical residues that are allowed and disallowed:
whitelist = ['MSE']

blacklist = D_isomer_AA + ['SEP', 'SLZ', 'TPO', '1MA', '2MG', '5MC',
                           '5MU', '7MG', 'H2U', 'M2G', 'OMC', 'OMG',
                           'PSU', 'HPH', 'CSO', 'CME', 'NLE', 'YCM',
                           'OIL', 'ABA', 'GOA', 'SMC', 'SNC', 'MHO',
                           'M3L', 'OCS', 'PCA', 'CSR', 'ACE', 'CAF',
                           'OCS', 'SEB', '4IN', 'CSS', '4FB', 'XX1',
                           'SNN', 'ACY', 'CYQ', 'NIY', 'BIF', 'CYQ',
                           'MES', 'CSD', 'ASB', 'KCX']


# For silencing BioPDB/DSSP/other output. Taken from:
# http://marc-abramowitz.com/archives/2013/07/19/python-context-manager-for-redirected-stdout-and-stderr/@contextlib.contextmanager
@contextlib.contextmanager
def stdchannel_redirected(stdchannel, dest_filename):
    """
    A context manager to temporarily redirect stdout or stderr

    e.g.:


    with stdchannel_redirected(sys.stderr, os.devnull):
        if compiler.has_function('clock_gettime', libraries=['rt']):
            libraries.append('rt')
    """

    try:
        oldstdchannel = os.dup(stdchannel.fileno())
        dest_file = open(dest_filename, 'w')
        os.dup2(dest_file.fileno(), stdchannel.fileno())

        yield
    finally:
        if oldstdchannel is not None:
            os.dup2(oldstdchannel, stdchannel.fileno())
        if dest_file is not None:
            dest_file.close()


def valid_PDB_dict(pdb_folder, res, cache_dir):
    '''
    Validate all PDB files according to a resolution criterium and that they are solved by X-RAY diffraction.
    Input pdb_folder is the folder that contains all the 2 letter folders splitting up all PDB entries according to
    the middle part of their name. A dictionary containing all the PDB ID's that have pass the quality control.
    '''
    # Go to the cache dir to fetch the possible cached dicts,
    # else create them a dump it here:
    os.chdir(cache_dir)
    # Make a timestamp:
    ts = time.time()
    # st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d')  # If things change from day to day
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m')
    # Use the timestamp in the filename to dump the results.
    # This should prevent to much recalculation:
    dict_name = 'QC_PDB_RES_' + str(res) + '_' + str(st) + '.p'

    # If the results are already calculated just load them and return:
    if os.path.isfile(dict_name):
        QC_PDB_dict = pickle.load(open(dict_name, "rb"))
        return(QC_PDB_dict)

    # Else validate all structures:
    QC_PDB_dict = dict()
    # Make an os.walk to go through all files in all folders:
    for dirname, dirnames, filenames in os.walk(pdb_folder):
        for filename in filenames:
            PDB_path = os.path.join(dirname, filename)
            PDB_id = PDB_path.split('/')[-1][3:7]
            xray = 0  # Swicth type variable
            resolution = 0   # Swicth type variable
            # Notice file getting opened is gzipped a file and therefore the input stream i binary:
            with gzip.open(PDB_path, 'r') as PDB_inf:
                lines = PDB_inf.readlines()
                for line in lines:
                    # Find the line that tells if the PBD file is from X-RAY diffraction,
                    # also translate the unicode string to binary before matching:
                    if line.startswith(b'EXPDTA    X-RAY DIFFRACTION'):
                        xray = 1
                    # Find the resolution line:
                    elif line.startswith(b'REMARK   2 RESOLUTION.'):
                        # Validate that the resolution is under the cutoff:
                        res = line[22:29]
                        res_str = str(res, 'utf-8')
                        try:
                            res_float = float(res_str)
                            if res_float < res:
                                resolution = 1
                        except:
                            pass

            # A PDB ID is kept if both xray and resolution flags checks out:
            keep = xray * resolution  # Both switches must be on to keep
            QC_PDB_dict[PDB_id] = keep
    # Lastly the results dictionary are dumped with Pickle:
    pickle.dump(QC_PDB_dict, open(dict_name, "wb"))
    return(QC_PDB_dict)


def parse_biolip(biolip_fnam, cache_dir):
    '''
    Parse the BioLiP database from a csv text file into a dictionary.
    Use the concatenation of the PDB ID and chain name as a dictionary key
    and a sorted list of all the associated ligands in a list as values.
    '''
    # Go to the cache dir to fetch the possible cached dicts,
    # else create them a dump it here:
    os.chdir(cache_dir)
    # Make a timestamp:
    ts = time.time()
    # st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d')  # If things change from day to day
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m')
    # Use the timestamp in the filename to dump the results.
    # This should prevent to much recalculation:
    biolip_dict_name = 'biolip_dict_' + str(st) + '.p'

    # If the results are already calculated just load them and return:
    if os.path.isfile(biolip_dict_name):
        biolip_dict = pickle.load(open(biolip_dict_name, "rb"))
        return(biolip_dict)

    # Do the parsing:
    biolip_dict = dict()
    with open(biolip_fnam, encoding='utf-8') as infile:
        infile_lines = infile.readlines()
        # Remove the header:
        infile_lines.pop(0)
        for line in infile_lines:
            line = line.rstrip('\n')
            cols = line.split(',')
            # Strip the annoying double quotes that comes with the csv file as default:
            cols = [el.strip('"') for el in cols]
            # Dictionary key is PDB ID + chain name:
            chainID = cols[1].lower() + cols[2].upper()
            ligand = cols[4]
            # Add the ligand to the list of ligands:
            if chainID in biolip_dict:
                biolip_dict[chainID].append(ligand)
            else:
                biolip_dict[chainID] = [ligand]

    # Sort the list of ligands for each chainID.
    # This makes the ligands lists comparable between PDB chains:
    for chainID in biolip_dict:
        biolip_dict[chainID].sort()

    # Lastly the resulting dictionary are dumped with Pickle and returned:
    pickle.dump(biolip_dict, open(biolip_dict_name, "wb"))
    return(biolip_dict)


def parse_ss_dis(ss_dis_file, cache_dir):
    '''
    Parse the ss_dis.txt file provided by RCSB which contains three things per PDB ID:
    1) the primary sequence
    2) the secondary sequence
    3) missing residues (disorder)
    Returns a nested dictionary which first key being the PDB ID
    and second key being the information type (1, 2 or 3).
    '''
    # Go to the cache dir to fetch the possible cached dicts,
    # else create them a dump it here:
    os.chdir(cache_dir)
    # Make a timestamp:
    ts = time.time()
    # st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d')  # If things change from day to day
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m')
    # Use the timestamp in the filename to dump the results.
    # This should prevent to much recalculation:
    ss_dis_dict_name = 'ss_dis_dict_' + str(st) + '.p'

    # If the results are already calculated just load them and return:
    if os.path.isfile(ss_dis_dict_name):
        ss_dis_dict = pickle.load(open(ss_dis_dict_name, "rb"))
        return(ss_dis_dict)

    # Do the parsing:
    ss_dis_dict = dict()
    with open(ss_dis_file, 'r') as infile:
        infile_lines = infile.readlines()
        for line in infile_lines:
            line = line.rstrip('\n')
            if line.startswith('>'):
                # Keep the PDB IDs in lower case letters with chain identifier as upper case:
                ID = line[1:5].lower() + line[6].upper()
                if ID not in ss_dis_dict:
                    ss_dis_dict[ID] = dict()
                # Find the sequence type 1, 2 or 3, and then make this key in the dictionary:
                stype = line[8:len(line)]
                ss_dis_dict[ID][stype] = ''
            else:
                ss_dis_dict[ID][stype] += line

    # Lastly the resulting dictionary are dumped with Pickle and returned:
    pickle.dump(ss_dis_dict, open(ss_dis_dict_name, "wb"))
    return(ss_dis_dict)


def seq_liglen_bin(min_seq_len, QC_PDB_dict, ss_dis_dict, biolip_dict, cache_dir):
    '''
    Bin PDB amino acid sequences based on their length.
    Also throw away PDBs that could not be validated
    because of not being high enough resolution, not long enough etc.
    '''
    # Go to the cache dir to fetch the possible cached dicts,
    # else create them a dump it here:
    os.chdir(cache_dir)
    # Make a timestamp:
    ts = time.time()
    # st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d')  # If things change from day to day
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m')
    # Use the timestamp in the filename to dump the results.
    # This should prevent to much recalculation:
    seqlen_dict_name = 'seq_liglen_dict_' + str(st) + '.p'

    # If the results are already calculated just load them and return:
    if os.path.isfile(seqlen_dict_name):
        seq_liglen_dict = pickle.load(open(seqlen_dict_name, "rb"))
        return(seq_liglen_dict)

    # Reverse (sequence: chainID) and dedup the ss_dis_dict:
    seq2chains = dedup_ss_dis(ss_dis_dict)
    seq_liglen_dict = dict()  # Nested dict: dict[ligkey][len][chainID] = sequence
    bad_list = list()         # For debugging only, keeps track of bad PDB chains
    not_in_QC = list()        # For debugging only, keeps track of bad PDB chains

    for seq, chainIDs in seq2chains.items():
        # Replace SEC (U) and PYL (O) with X:
        seq = seq.replace("O", "X")
        seq = seq.replace("U", "X")
        # Require a pretty long protein to avoid small unstructured proteins:
        if len(seq) < min_seq_len:
            bad_list.append(chainIDs)
            continue
        # Require canonical amino acids or X for hetero atoms.
        # PDB entries with Z, B or J are very few and bad entries so discard them:
        if re.search('[^ACDEFGHIKLMNPQRSTVWYX]', seq):
            bad_list.append(chainIDs)
            continue

        # Now find the chains to keep, as annotated in the QC dictionary:
        chainID_list = chainIDs.split('-')
        keep_list = list()
        for chainID in chainID_list:
            pdb = chainID[0:-1]
            if pdb not in QC_PDB_dict:
                not_in_QC.append(pdb)
            elif QC_PDB_dict[pdb]:
                keep_list.append(chainID)

        # If no validated PDB IDs were appended to the list put the sequence into the "bad_list":
        if len(keep_list) == 0:
            bad_list.append(chainIDs)
            continue

        # Build the seq_liglen_dict, first pass will be with multiple identical sequences
        # having the same lig and len keys. This will be deduped in anther pass:
        seq_len = len(seq)
        for chainID in keep_list:
            # Define the ligang key:
            lig_key = ''
            if chainID in biolip_dict:
                lig_key = '-'.join(biolip_dict[chainID])
            else:
                lig_key = 'no_lig'

            if lig_key in seq_liglen_dict:
                if seq_len in seq_liglen_dict[lig_key]:
                    seq_liglen_dict[lig_key][seq_len][chainID] = seq
                else:
                    seq_liglen_dict[lig_key][seq_len] = dict()
                    seq_liglen_dict[lig_key][seq_len][chainID] = seq
            else:
                seq_liglen_dict[lig_key] = dict()
                if seq_len in seq_liglen_dict[lig_key]:
                    seq_liglen_dict[lig_key][seq_len][chainID] = seq
                else:
                    seq_liglen_dict[lig_key][seq_len] = dict()
                    seq_liglen_dict[lig_key][seq_len][chainID] = seq

    # Now do the deduping in a second pass:
    for lig_key, seq_len_dict in seq_liglen_dict.items():
        for len_key, seq_chain_dict in seq_len_dict.items():
            # Keep all the chainID's with the same sequence in a dictionary:
            seq2chains = dict()
            for chainID, seq in seq_chain_dict.items():
                if seq in seq2chains:
                    # Store all PDB IDs with same sequences as dash separated values:
                    seq2chains[seq] = '{}-{}'.format(seq2chains[seq], chainID)
                else:
                    seq2chains[seq] = chainID

            # Inverse the seq2chains dictionary:
            inv_seq2chains = {v: k for k, v in seq2chains.items()}
            seq_liglen_dict[lig_key][len_key] = dict()
            seq_liglen_dict[lig_key][len_key] = inv_seq2chains

    # Lastly the results dictionary are dumped with Pickle:
    pickle.dump(seq_liglen_dict, open(seqlen_dict_name, "wb"))
    return(seq_liglen_dict)


def find_chain_pairs(pair_dist, seq_len_dict, cache_dir):
    '''
    Calculate the hamming distance between sequences of the same length.
    Then return a list of the pairs that have the requested distance.
    The requested distance will, in most cases, be 1
    and the pairs returned will have a single point variant between them.
    Deprecated, but keep for using if a seq_len_dict is extracted from a seq_liglen_dict.
    Like: seq_liglen_dict['no_lig'] = seq_len_dict
    '''
    # Go to the cache dir to fetch the possible cached dicts,
    # else create them a dump it here:
    os.chdir(cache_dir)
    # Make a timestamp:
    ts = time.time()
    # st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d')
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m')  # If things change from day to day
    # Use the timestamp in the filename to dump the results.
    # This should prevent to much recalculation:
    pair_dict_name = 'mut_pairs' + str(pair_dist) + '_' + str(st) + '.p'
    pair_dist_name = 'mut_pairs_dist' + str(pair_dist) + '_' + str(st) + '.p'

    # If the results are already calculated just load them and return:
    if os.path.isfile(pair_dict_name) and os.path.isfile(pair_dist_name):
        pairs = pickle.load(open(pair_dict_name, "rb"))
        distance_distribution = pickle.load(open(pair_dist_name, "rb"))
        return(pairs, distance_distribution)

    pairs = list()  # List of tuples with the pairs and mutation positions
    distance_distribution = dict()  # Only for debugging and/or interesting stats
    # Pairwise seqeunce alignment to find the single substitution chain pairs:
    for len_key, len_dict in sorted(seq_len_dict.items(), key=lambda x: x[0], reverse=False):
        # Make two identical lists of the PDB_IDs/seq keys:
        PDB_ID_list = list(len_dict.keys())
        com_list = PDB_ID_list[:]
        for PDB_ID in PDB_ID_list:
            # Pop the top element of one of the lists to avoid double calculation
            # (under and over the diagonal of a all against all comparison matrix):
            com_list.pop(0)
            seq1 = seq_len_dict[len_key][PDB_ID]
            # Make the comparison:
            for PDB_ID_com in com_list:
                seq2 = seq_len_dict[len_key][PDB_ID_com]
                # Calculate the hamming distance and mismatch positions, disregarding X's as mismatches:
                pos_list = mutation_positions_noX(seq1, seq2)
                ham_dist = len(pos_list)
                # Add this information to the distribution of distances:
                if ham_dist in distance_distribution:
                    distance_distribution[ham_dist] += 1
                else:
                    distance_distribution[ham_dist] = 1
                # If this is the requested distance add the pair to the list:
                if ham_dist == pair_dist:
                    pairs.append((PDB_ID, PDB_ID_com, pos_list))

    # Lastly the results dictionary are dumped with Pickle:
    pickle.dump(pairs, open(pair_dict_name, "wb"))
    pickle.dump(distance_distribution, open(pair_dist_name, "wb"))
    return(pairs, distance_distribution)


def find_chain_pairs2(pair_dist, seq_liglen_dict, cache_dir):
    '''
    Calculate the hamming distance between sequences of the same length.
    Then return a list of the pairs that have the requested distance.
    The requested distance will, in most cases, be 1
    and the pairs returned will have a single point variant between them.
    '''
    # Go to the cache dir to fetch the possible cached dicts,
    # else create them a dump it here:
    os.chdir(cache_dir)

    # Make a timestamp:
    ts = time.time()
    # st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d')  # If things change from day to day
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m')
    # Use the timestamp in the filename to dump the results.
    # This should prevent to much recalculation:
    pair_dict_name = 'mut_pairs_v2_' + str(pair_dist) + '_' + str(st) + '.p'
    pair_dist_name = 'mut_pairs_dist_v2_' + str(pair_dist) + '_' + str(st) + '.p'

    # If the results are already calculated just load them and return:
    if os.path.isfile(pair_dict_name) and os.path.isfile(pair_dist_name):
        pairs = pickle.load(open(pair_dict_name, "rb"))
        distance_distribution = pickle.load(open(pair_dist_name, "rb"))
        return(pairs, distance_distribution)

    pairs = list()  # List of tuples with the pairs and mutation positions
    distance_distribution = dict()  # Only for debugging and/or interesting stats
    for lig_key, seq_len_dict in seq_liglen_dict.items():
        # Pairwise seqeunce alignment to find the single substitution chain pairs:
        for len_key, len_dict in sorted(seq_len_dict.items(), key=lambda x: x[0], reverse=False):
            # Make two identical lists of the PDB_IDs/seq keys:
            PDB_ID_list = list(len_dict.keys())
            com_list = PDB_ID_list[:]
            for PDB_ID in PDB_ID_list:
                # Pop the top element of one of the lists to avoid double calculation
                # (under and over the diagonal of a all against all comparison matrix):
                com_list.pop(0)
                seq1 = seq_len_dict[len_key][PDB_ID]
                # Make the comparison:
                for PDB_ID_com in com_list:
                    seq2 = seq_len_dict[len_key][PDB_ID_com]
                    pos_list = mutation_positions_noX(seq1, seq2)
                    ham_dist = len(pos_list)
                    # Calculate the hamming distance and mismatch positions, disregarding X's as mismatches:
                    if ham_dist in distance_distribution:
                        distance_distribution[ham_dist] += 1
                    else:
                        distance_distribution[ham_dist] = 1
                    # If this is the requested distance add the pair to the list:
                    if ham_dist == pair_dist:
                        pairs.append((PDB_ID, PDB_ID_com, pos_list))

    # Lastly the results dictionary are dumped with Pickle:
    pickle.dump(pairs, open(pair_dict_name, "wb"))
    pickle.dump(distance_distribution, open(pair_dist_name, "wb"))
    return(pairs, distance_distribution)


def dedup_ss_dis(ss_dis_dict):
    '''
    Reverse the ss_dis dictionary to get the primary sequence as keys
    and all the PDB IDs having this sequences as a single string value.
    '''
    seq2chains = dict()
    for chainID, subdict in ss_dis_dict.items():
        seq = subdict['sequence']
        if seq in seq2chains:
            # Store all PDB IDs with same sequences as dash separated values:
            seq2chains[seq] = '{}-{}'.format(seq2chains[seq], chainID)
        else:
            seq2chains[seq] = chainID
    return(seq2chains)


# To be deleted:
def dedup_seq_liglen_dict(seq_liglen_dict):
    '''
    Currently deprecated since this function is moved into the seq_liglen_bin function.
    '''
    for lig_key, seq_len_dict in seq_liglen_dict.items():
        for len_key, seq_chain_dict in seq_len_dict.items():
            seq2chains = dict()
            for chainID, seq in seq_chain_dict.items():
                if seq in seq2chains:
                    # Store all PDB IDs with same sequences as dash separated values:
                    seq2chains[seq] = '{}-{}'.format(seq2chains[seq], chainID)
                else:
                    seq2chains[seq] = chainID

            # Inverse the seq2chains dictionary:
            inv_seq2chains = {v: k for k, v in seq2chains.items()}
            seq_liglen_dict[lig_key][len_key] = dict()
            seq_liglen_dict[lig_key][len_key] = inv_seq2chains

    return(seq_liglen_dict)


# To be deleted:
def dedup_seq_len_dict(seq_len_dict):
    '''
    Currently deprecated since this function is moved into the seq_liglen_bin function.
    '''
    for len_key, seq_chain_dict in seq_len_dict.items():
        seq2chains = dict()
        for chainID, seq in seq_chain_dict.items():
            if seq in seq2chains:
                # Store all PDB IDs with same sequences as dash separated values:
                seq2chains[seq] = '{}-{}'.format(seq2chains[seq], chainID)
            else:
                seq2chains[seq] = chainID

        # Inverse the seq2chains dictionary:
        inv_seq2chains = {v: k for k, v in seq2chains.items()}
        seq_len_dict[len_key] = dict()
        seq_len_dict[len_key] = inv_seq2chains

    return(seq_len_dict)


def fasta2dict(fasta_file):
    '''
    Read fasta file into dictionary, header as key and sequence as value.
    '''
    fasta_dict = dict()
    with open(fasta_file, 'r') as infile:
        infile_lines = infile.readlines()
        for line in infile_lines:
            line = line.strip()
            if line.startswith('>'):
                head = line
                fasta_dict[head] = ''
            else:
                fasta_dict[head] += line
    return(fasta_dict)


def mutation_positions_noX(seq1, seq2):
    '''
    Find and return the positions of difference between two strings.
    This is used for finding the position of the mutation on sequence level.
    For finding it on crystal level potential missing residues must be taken into account
    and therefore X residues are considered wildcards.
    '''
    assert(len(seq1) == len(seq2))
    # Converted to a more simple list comprehension:
    pos_list = [pos for pos in range(0, len(seq1)) if seq1[pos] != seq2[pos]]

# To be deleted:
#    pos_list = list()
#    for pos in range(0, len(seq1)):
#        if seq1[pos] == 'X' or seq2[pos] == 'X':
#            pass
#        elif seq1[pos] != seq2[pos]:
#            pos_list.append(pos)
#        else:
#            pass
    return(pos_list)


def create_pair_folder(scratch_dir, pair_number, pair, pdb_folder):
    '''
    Creates a folder for each pair and moves the necessary files
    to do all calculations in that local folder and hence enable parallelization.
    '''
    # Name the folder according to the pair number:
    pair_folder = '{}/{:06d}'.format(scratch_dir, pair_number)
    # Remove existing folder if any:
    if os.path.exists(pair_folder):
        shutil.rmtree(pair_folder)
    os.makedirs(pair_folder)

    # Make a list of all the PDBs in the pair tuple:
    pdbs = pair[0].split('-') + pair[1].split('-')
    # Then loop through all these PDBs and move them to the folder:
    for pdb_file in pdbs:
        chain_id = pdb_file[-1]
        pdb_file = pdb_file[0:-1]   # Cut away the chain ID
        folder_key = pdb_file[1:3]  # Notice this is PDB convention and for custom files this should be a adapted
        pdb_path = '{}/{}/pdb{}.ent.gz'.format(pdb_folder, folder_key, pdb_file)
        fh_out = open('{}/{}'.format(pair_folder, pdb_file), 'w')
        # Notice this is opened directly as gzipped and therefore comes in binary:
        try:  # If the file exists
            with gzip.open(pdb_path, 'r') as fh_in:
                lines = fh_in.readlines()
                for line in lines:
                    # Decode from binary to unicode:
                    line = str(line, 'utf-8')
                    # Then only keep 'ATOM' cards or HETATM on the whitelist:
                    if line[0:6] == 'ATOM  ' or (line[0:6] == 'HETATM' and line[17:20] in whitelist):
                        print(line, file=fh_out, end='')
                    # If residue in blacklist skip the pair:
                    elif line[0:6] == 'HETATM' and line[21] == chain_id and line[17:20] in blacklist:
                        print('Residue found in blacklist', pdb_file)
                        return(False)
        except Exception as e:
            print(e)
            print(pdb_path, 'probably missing')
            return(False)

        fh_out.close()
    return(pair_folder)


def correct_pdb_obj_pairs(res_list1, res_list2, mut_pos_seq, ss_dis_dict, id1, id2):
    '''Correct missing residues in two pdb objects simultaneously by trimming the sequence in both
    if a residue is missing. Return the two correted lists of residues
    and an updated mutation position relative to the residue lists.
    '''
    # Get the primary sequence of each pdb object as it is given by the ss_dis dictionary:
    seq1 = ss_dis_dict[id1]['sequence']
    seq2 = ss_dis_dict[id2]['sequence']
    # Also find the amino acid on the mutation position:
    mAA_seq_1 = seq1[mut_pos_seq]
    mAA_seq_2 = seq2[mut_pos_seq]
    assert(len(seq1) == len(seq2))
    # And the disorder between sequence and crystal structure:
    dis1 = ss_dis_dict[id1]['disorder']
    dis2 = ss_dis_dict[id2]['disorder']
    assert(len(dis1) == len(dis2))
    assert(len(dis1) == len(seq1))

    # Assert any oddities in the disorder string and/or crystal structure:
    res_count1 = dis1.count('-')
    res_count2 = dis2.count('-')
    X1_count = seq1.count('X')
    X2_count = seq2.count('X')
    assert((len(res_list1) + X1_count) == res_count1)
    assert((len(res_list2) + X2_count) == res_count2)

# To delete:
##### Remove this out-commented part if things are running as usual:
    # Now trim the two residue obejcts:
#    res_idx = 0
#    for i in range(0, len(dis1)):
#        if dis1[i] == 'X' and dis2[i] == 'X':    # Both are missing
#            continue
#            res_idx -= 1
#        elif dis1[i] == 'X' and seq2[i] != 'X':  # Missing residue in structure 1, but found in structure 2
#            del res_list2[res_idx]
#            res_idx -= 1
#            #if   ## something with "mut_pos"
#        elif dis2[i] == 'X' and seq1[i] != 'X':  # Missing residue in structure 2, but found in structure 1
#            del res_list1[res_idx]
#            res_idx -= 1
#        elif seq1[i] == 'X' and dis2[i] != 'X':  # Missing residue in structure 1, but found in structure 2
#            del res_list2[res_idx]
#            res_idx -= 1
#        elif seq2[i] == 'X' and dis1[i] != 'X':  # Missing residue in structure 2, but found in structure 1
#            del res_list1[res_idx]
#            res_idx -= 1

#        res_idx += 1

    # Now trim the two residue obejcts:
    res_idx = 0
    for i in range(0, len(dis1)):
        # Both are missing
        if dis1[i] == 'X' and dis2[i] == 'X':
            continue
            res_idx -= 1
        # Missing residue in structure 1, but found in structure 2:
        elif (dis1[i] == 'X' and seq2[i] != 'X') or (seq1[i] == 'X' and dis2[i] != 'X'):
            del res_list2[res_idx]
            res_idx -= 1
        # Missing residue in structure 2, but found in structure 1:
        elif (dis2[i] == 'X' and seq1[i] != 'X') or (seq2[i] == 'X' and dis1[i] != 'X'):
            del res_list1[res_idx]
            res_idx -= 1

        res_idx += 1

    # Get the trimmed sequence as a string:
    res_seq_list1 = [residue_type_3to1_map[res_obj.get_resname()] for res_obj in res_list1]
    res_list_seq1 = ''.join(res_seq_list1)
    # Then join to a string:
    res_seq_list2 = [residue_type_3to1_map[res_obj.get_resname()] for res_obj in res_list2]
    res_list_seq2 = ''.join(res_seq_list2)

    # Check that the mutation is still there and not trimmed away:
    try:
        mut_list_crystal = mutation_positions_noX(res_list_seq1, res_list_seq2)
    except:
        if len(mut_list_crystal) == 0:
            return('no_mut', 'no_mut', 'no_mut')

    # Use the amino acid identity on the mutation position as and anchor point,
    # to check whether everything went as it should in the trimming:
    mut_pos_crystal = mut_list_crystal[0]
    mAA_crystal_1 = res_list_seq1[mut_pos_crystal]
    mAA_crystal_2 = res_list_seq2[mut_pos_crystal]
    assert(mAA_seq_1 == mAA_crystal_1)
    assert(mAA_seq_2 == mAA_crystal_2)

    return(res_list1, res_list2, mut_pos_crystal)


def Calpha_from_residue_list(residue_list):
    '''
    Returns a list of the carbon alpha objects from a biopython PDB object.
    '''
    Ca_list = list()
    for residue in residue_list:
        if 'CA' in residue:
            Ca_list.append(residue['CA'])
    return(Ca_list)


def Calpha_from_pair_residue_list(res_list1, res_list2):
    '''
    Returns two lists of the carbon shared alpha objects between two biopython PDB object.
    '''
    assert(len(res_list1) == len(res_list2))
    Ca_list1 = list()
    Ca_list2 = list()
    for idx in range(0, len(res_list1)):
        if 'CA' in res_list1[idx] and 'CA' in res_list2[idx]:
            Ca_list1.append(res_list1['CA'])
            Ca_list2.append(res_list2['CA'])
    return(Ca_list1, Ca_list2)


def BBatom_from_residue_list(residue_list):
    '''
    Returns a list of the backbone atom objects from a biopython PDB object.
    '''
    Ba_list = list()
    for residue in residue_list:
        if 'N' in residue:
            Ba_list.append(residue['N'])
        if 'CA' in residue:
            Ba_list.append(residue['CA'])
        if 'C' in residue:
            Ba_list.append(residue['C'])
        if 'O' in residue:
            Ba_list.append(residue['O'])
    return(Ba_list)


def BBatom_from_pair_residue_list(res_list1, res_list2):
    '''
    Returns two lists of the backbone atom shared between two Biopython PDB object.
    '''
    Ba_list1 = list()
    Ba_list2 = list()
    for idx in range(0, len(res_list1)):
        if 'N' in res_list1[idx] and 'N' in res_list2[idx]:
            Ba_list1.append(res_list1[idx]['N'])
            Ba_list2.append(res_list2[idx]['N'])
        if 'C' in res_list1[idx] and 'C' in res_list2[idx]:
            Ba_list1.append(res_list1[idx]['C'])
            Ba_list2.append(res_list2[idx]['C'])
        if 'CA' in res_list1[idx] and 'CA' in res_list2[idx]:
            Ba_list1.append(res_list1[idx]['CA'])
            Ba_list2.append(res_list2[idx]['CA'])
        if 'O' in res_list1[idx] and 'O' in res_list2[idx]:
            Ba_list1.append(res_list1[idx]['O'])
            Ba_list2.append(res_list2[idx]['O'])
    return(Ba_list1, Ba_list2)


def rot_residue_list(residue_list, rotran):
    '''
    Rotate all the atoms of a given PDB structure loaded into a biopython PDB obejct.
    '''
    # Here the "rotran" is a rotational matrix from the superposition protocol:
    rot, tran = rotran
    rot = rot.astype('f')
    tran = tran.astype('f')
    for residue in residue_list:
        for atom in residue:
            atom.transform(rot, tran)


def calc_phi(residue_list, residue_numb):
    '''
    Calculate and return the phi dihedral angel of a residue.
    '''
    if 'C' in residue_list[residue_numb - 1] and 'N' in residue_list[residue_numb] and 'CA' in residue_list[residue_numb] and 'C' in residue_list[residue_numb]:
        vector1 = residue_list[residue_numb - 1]['C'].get_vector()
        vector2 = residue_list[residue_numb]['N'].get_vector()
        vector3 = residue_list[residue_numb]['CA'].get_vector()
        vector4 = residue_list[residue_numb]['C'].get_vector()
        # calc_dihedral is part of Bio.PDB:
        phi_torsion = calc_dihedral(vector1, vector2, vector3, vector4)
        return(phi_torsion)
    else:
        # If any missing bacbone atoms skip this residue:
        return('NA')


def calc_psi(residue_list, residue_numb):
    '''
    Calculate and return the psi dihedral angel of a residue.
    '''
    if 'N' in residue_list[residue_numb] and 'CA' in residue_list[residue_numb] and 'C' in residue_list[residue_numb] and 'N' in residue_list[residue_numb + 1]:
        vector1 = residue_list[residue_numb]['N'].get_vector()
        vector2 = residue_list[residue_numb]['CA'].get_vector()
        vector3 = residue_list[residue_numb]['C'].get_vector()
        vector4 = residue_list[residue_numb + 1]['N'].get_vector()
        # calc_dihedral is part of Bio.PDB:
        psi_torsion = calc_dihedral(vector1, vector2, vector3, vector4)
        return(psi_torsion)
    else:
        # If any missing bacbone atoms skip this residue:
        return('NA')


def calc_model_torsions(residue_list):
    '''
    Calculate the phi/psi torsions for all residues
    except the first and last residue in a Biopython PDB object.
    '''
    torsion_list = list()
    last_residue_idx = len(residue_list) - 1
    for residue_idx, residue in enumerate(residue_list):
        if residue_idx > 0 and residue_idx < last_residue_idx:
            phi = calc_phi(residue_list, residue_idx)
            psi = calc_psi(residue_list, residue_idx)
            torsion_list.append((phi, psi))
        else:  # For first and last residue just append "blank" values
            torsion_list.append(('NA', 'NA'))
# To be deleted:
# So far just skip the terminal residues:
#        elif residue_idx == 0:
#            psi = calc_psi(chain_obj, residue_idx)
#        elif residue_idx == last_residue_idx:
#            phi = calc_phi(chain_obj, residue_idx)

    return(torsion_list)


def calc_diff_torsions(torsion1, torsion2):
    '''
    Calculate the difference in two lists of torsion angles
    as returned from the "calc_model_torsions" function.
    '''
    assert(len(torsion1) == len(torsion2))
    torsion_diff = list()
    for idx in range(0, len(torsion1)):
        phi1 = torsion1[idx][0]
        phi2 = torsion2[idx][0]
        psi1 = torsion1[idx][1]
        psi2 = torsion2[idx][1]
        # Skip those all torsions for a residue missing either one of them:
        if 'NA' in [phi1, phi2, psi1, psi2]:
            torsion_diff.append(('NA', 'NA'))
            continue

        # Calculate the radial distance between phi1 and phi2:
        if phi1 < 0 and phi2 > 0:
            phi1 = math.pi + -1 * phi1
        elif phi1 > 0 and phi2 < 0:
            phi2 = math.pi + -1 * phi2
        d_phi = abs(phi1 - phi2) * 360 / (2 * math.pi)
        # The difference cannot be greater than a half circle (180 degrees):
        if d_phi > 180:
            d_phi = 360 - d_phi

        # Calculate the radial distance between psi1 and psi2:
        if psi1 < 0 and psi2 > 0:
            psi1 = math.pi + -1 * psi1
        elif psi1 > 0 and psi2 < 0:
            psi2 = math.pi + -1 * psi2
        d_psi = abs(psi1 - psi2) * 360 / (2 * math.pi)
        # The difference cannot be greater than a half circle (180 degrees):
        if d_psi > 180:
            d_psi = 360 - d_psi

        torsion_diff.append((d_phi, d_psi))
    return(torsion_diff)


def align_and_find_torsions(res_list1, res_list2):
    '''
    Superimpose two biopython PDB objects and then calculate the difference in torsions for each residue.
    '''
    # Superimposer is part of Bio.PDB:
    sup = Superimposer()
    # First object is fixed the second is rotated:
    fixed, moving = BBatom_from_pair_residue_list(res_list1, res_list2)

    # Do the superposition and get the rotation matrix:
    sup.set_atoms(fixed, moving)
    rotran = sup.rotran
    # Rotate pdb object two:
    rot_residue_list(res_list2, rotran)

    # Calculate the phi/psi torsions:
    torsion1 = calc_model_torsions(res_list1)
    torsion2 = calc_model_torsions(res_list2)

    # Calculate the difference in torsions:
    torsion_diff = calc_diff_torsions(torsion1, torsion2)
    return(torsion_diff)


def dist_from_mut(res_list, mut_pos):
    '''
    Calculate the distance from all residues to the mutated position in the wild type protein.
    Distance measure is based on C-alpha distance in 3D space and positions in 1D space.
    '''
    dist_list = list()
    Calpha = Calpha_from_residue_list(res_list)
    for i in range(0, len(Calpha)):
        dist_3D = Calpha[i] - Calpha[mut_pos]
        dist_lin = abs(mut_pos - i)
        dist_list.append([dist_3D, dist_lin])
    return(dist_list)


def create_atom_line(atom, hetfield, segid, atom_number, resname, resseq, icode, chain_id, charge="  "):
    '''
    Create the ATOM line in a PDB file from the atom object and its parent pdb object information.
    This function is taken directly from the biopython source code and slightly modified.
    '''
    from Bio.Data.IUPACData import atom_weights  # Allowed Elements
    _ATOM_FORMAT_STRING = "%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%s%6.2f      %4s%2s%2s\n"

    if hetfield != " ":
        record_type = "HETATM"
    else:
        record_type = "ATOM  "
    if atom.element:
        element = atom.element.strip().upper()
        if element.capitalize() not in atom_weights:
            raise ValueError("Unrecognised element %r" % atom.element)
        element = element.rjust(2)
    else:
        element = "  "
    name = atom.get_fullname()
    altloc = atom.get_altloc()
    x, y, z = atom.get_coord()
    bfactor = atom.get_bfactor()
    occupancy = atom.get_occupancy()
    try:
        occupancy_str = "%6.2f" % occupancy
    except TypeError:
        if occupancy is None:
            occupancy_str = " " * 6
            import warnings
            from Bio import BiopythonWarning
            warnings.warn("Missing occupancy in atom %s written as blank" %
                          repr(atom.get_full_id()), BiopythonWarning)
        else:
            raise TypeError("Invalid occupancy %r in atom %r"
                            % (occupancy, atom.get_full_id()))

    args = (record_type, atom_number, name, altloc, resname, chain_id,
            resseq, icode, x, y, z, occupancy_str, bfactor, segid,
            element, charge)
    return(_ATOM_FORMAT_STRING % args)


def write_reslist_as_pdbfile(res_list, chain_id, outname):
    '''
    Dump a list of residue objects as a PDB file.
    '''
    atom_number = 1
    with open(outname, 'w') as fh:
        for residue in res_list:
            hetfield, resseq, icode = residue.get_id()
            resname = residue.get_resname()
            segid = residue.get_segid()
            for atom in residue.get_unpacked_list():
                s = create_atom_line(atom, hetfield, segid, atom_number, resname,
                                     resseq, icode, chain_id)
                fh.write(s)
                atom_number += 1


def add_dssp_to_reslist(pair_folder, single_pair, res_list1, res_list2, parser):
    '''
    Run DSSP on a list of residue objects by first dumping them as ATOM cards to a PDB file
    and then running DSSP as a subprocess.
    Finally add this DSSP information to each residue object and return the updated residue lists.
    '''
    # First dump the list of residues as a PDB file:
    dest1 = pair_folder + '/' + single_pair[0] + '_trimmed'
    dest2 = pair_folder + '/' + single_pair[1] + '_trimmed'
    write_reslist_as_pdbfile(res_list1, 'A', dest1)
    write_reslist_as_pdbfile(res_list2, 'A', dest2)

    # Then parse them as pdb objects:
    s1 = parser.get_structure(single_pair[0], dest1)
    m1 = s1[0]
    s2 = parser.get_structure(single_pair[1], dest2)
    m2 = s2[0]

    # And run DSSP:
    # Notice that the two returned dictionaries are passed to throwaway variables.
    # DSSP is part of Bio.PDB:
    _dssp1 = DSSP(m1, dest1, 'dssp', 'Wilke', 'PDB')
    _dssp2 = DSSP(m2, dest2, 'dssp', 'Wilke', 'PDB')

    # Extract the list of residue objects:
    res_list1_dssp = list(list(s1[0].get_chains())[0].get_residues())
    res_list2_dssp = list(list(s2[0].get_chains())[0].get_residues())
    assert(len(res_list1_dssp) == len(res_list1))
    assert(len(res_list2_dssp) == len(res_list2))

    # Add the xtra data from the old residue lists:
    for idx in range(len(res_list1_dssp)):
        # Store the monomeric DSSP values by another name:
        # res_list1_dssp[idx].xtra['SS_DSSP_trim'] = res_list1_dssp[idx].xtra['SS_DSSP']  # The real secondary structure annotation comes from the full crystal
        # res_list2_dssp[idx].xtra['SS_DSSP_trim'] = res_list2_dssp[idx].xtra['SS_DSSP']  # The real secondary structure annotation comes from the full crystal
        if 'EXP_DSSP_ASA' not in res_list1_dssp[idx].xtra:
            print('res_list1, idx:', idx)
        else:
            res_list1_dssp[idx].xtra['EXP_DSSP_ASA_trim'] = res_list1_dssp[idx].xtra['EXP_DSSP_ASA']
        if 'EXP_DSSP_ASA' not in res_list2_dssp[idx].xtra:
            print('res_list1, idx:', idx)
        else:
            res_list2_dssp[idx].xtra['EXP_DSSP_ASA_trim'] = res_list2_dssp[idx].xtra['EXP_DSSP_ASA']
        res_list1_dssp[idx].xtra['EXP_DSSP_RASA_trim'] = res_list1_dssp[idx].xtra['EXP_DSSP_RASA']
        res_list2_dssp[idx].xtra['EXP_DSSP_RASA_trim'] = res_list2_dssp[idx].xtra['EXP_DSSP_RASA']
        # Then update the DSSP values with the full PDB entry values:
        res_list1_dssp[idx].xtra.update(res_list1[idx].xtra)
        res_list2_dssp[idx].xtra.update(res_list2[idx].xtra)

    return(res_list1_dssp, res_list2_dssp)


def add_dssp_to_pdb_obj(pair_folder, single_pair, pdb_obj1, pdb_obj2):
    '''
    Run DSSP on a pdb objects created from a PDB file.
    Add this DSSP information to each residue object and return the updated pdb object.
    '''
    # Get filename to run DSSP on:
    dest1 = pair_folder + '/' + single_pair[0][0:-1]
    dest2 = pair_folder + '/' + single_pair[1][0:-1]

    # Run DSSP:
    # Notice that the two returned dictionaries are parsed to throwaway variables.
    # DSSP is part of Bio.PDB:
    _dssp1 = DSSP(pdb_obj1[0], dest1, 'dssp', 'Wilke', 'PDB')
    _dssp2 = DSSP(pdb_obj2[0], dest2, 'dssp', 'Wilke', 'PDB')

    return(pdb_obj1, pdb_obj2)


def bfactor_sidechain(res_obj):
    '''
    Calculate the average b-factor of the side chain atoms in a biopython residue obejct.
    '''
    bfactor = 0
    count = 0
    for atom in res_obj:
        # Continue if backbone atom:
        if atom.get_name() in ['CA', 'C', 'N', 'O']:
            continue
        bfactor += atom.get_bfactor()
        count += 1

    if count > 0:
        bfactor /= count
    else:
        bfactor = 'NA'
    return(bfactor)


def bfactor_backbone(res_obj):
    '''
    Calculate the average b-factor of the backbone chain atoms in a biopython residue obejct.
    '''
    bfactor = 0
    count = 0
    for atom in res_obj:
        # Continue if sidechain atom:
        if atom.get_name() not in ['CA', 'C', 'N', 'O']:
            continue
        bfactor += atom.get_bfactor()
        count += 1

    if count > 0:
        bfactor /= count
    else:
        bfactor = 'NA'
    return(bfactor)


def print_res(single_pair, mut_dist, torsion_diff, res_list1, res_list2, mut_pos_crystal):
    '''
    Builds up a long string of all the formated results for a single pair.
    Then returns this string to let other processes print it.
    '''
    print_str = ''
    # Define all the static information, specific for the mutation,
    # but common for all the residues in the pair:
    pair_name = '-'.join(single_pair)
    prot_len = len(mut_dist)
    AA_wt = res_list1[mut_pos_crystal].get_resname()
    AA_mut = res_list2[mut_pos_crystal].get_resname()
    ss_wt = res_list1[mut_pos_crystal].xtra['SS_DSSP']
    ss_mut = res_list2[mut_pos_crystal].xtra['SS_DSSP']
    ASA_wt = res_list1[mut_pos_crystal].xtra['EXP_DSSP_ASA']
    ASA_mut = res_list2[mut_pos_crystal].xtra['EXP_DSSP_ASA']
    RASA_wt = res_list1[mut_pos_crystal].xtra['EXP_DSSP_RASA']
    RASA_mut = res_list2[mut_pos_crystal].xtra['EXP_DSSP_RASA']
    ASA_trim_wt = res_list1[mut_pos_crystal].xtra['EXP_DSSP_ASA_trim']
    ASA_trim_mut = res_list2[mut_pos_crystal].xtra['EXP_DSSP_ASA_trim']
    RASA_trim_wt = res_list1[mut_pos_crystal].xtra['EXP_DSSP_RASA_trim']
    RASA_trim_mut = res_list2[mut_pos_crystal].xtra['EXP_DSSP_RASA_trim']
    ichain_numb_wt = res_list1[mut_pos_crystal].xtra['ichain_numb']
    ichain_numb_mut = res_list2[mut_pos_crystal].xtra['ichain_numb']
    idist_wt = res_list1[mut_pos_crystal].xtra['idist']
    idist_mut = res_list2[mut_pos_crystal].xtra['idist']
    ires_wt = res_list1[mut_pos_crystal].xtra['ires']
    ires_mut = res_list2[mut_pos_crystal].xtra['ires']
    iatom_self_wt = res_list1[mut_pos_crystal].xtra['iatom_self']
    iatom_self_mut = res_list2[mut_pos_crystal].xtra['iatom_self']
    iatom_wt = res_list1[mut_pos_crystal].xtra['iatom']
    iatom_mut = res_list2[mut_pos_crystal].xtra['iatom']

    # Then add the residue specific information.
    # Skip first and last residue since these do not have both phi and psi:
    for i in range(1, prot_len - 1):
        # Define all the residue specific information:
        _3D_dist = mut_dist[i][0]
        _1D_dist = mut_dist[i][1]
        N_dist = i
        C_dist = prot_len - (i + 1)
        dphi = torsion_diff[i][0]
        dpsi = torsion_diff[i][1]
        AA_res = res_list1[i].get_resname()
        ss1_res = res_list1[i].xtra['SS_DSSP']
        ss2_res = res_list2[i].xtra['SS_DSSP']
        ASA1_res = res_list1[i].xtra['EXP_DSSP_ASA']
        ASA2_res = res_list2[i].xtra['EXP_DSSP_ASA']
        RASA1_res = res_list1[i].xtra['EXP_DSSP_RASA']
        RASA2_res = res_list2[i].xtra['EXP_DSSP_RASA']
        ASA1_trim_res = res_list1[i].xtra['EXP_DSSP_ASA_trim']
        ASA2_trim_res = res_list2[i].xtra['EXP_DSSP_ASA_trim']
        RASA1_trim_res = res_list1[i].xtra['EXP_DSSP_RASA_trim']
        RASA2_trim_res = res_list2[i].xtra['EXP_DSSP_RASA_trim']
        bfactor_back1_res = bfactor_backbone(res_list1[i])
        bfactor_back2_res = bfactor_backbone(res_list2[i])
        bfactor_side1_res = bfactor_sidechain(res_list1[i])
        bfactor_side2_res = bfactor_sidechain(res_list2[i])
        next2cut = res_list1[i].xtra['next2cut']  # Same for both res_list 1 and 2
        idist1_res = res_list1[i].xtra['idist']
        idist2_res = res_list2[i].xtra['idist']
        ires1_res = res_list1[i].xtra['ires']
        ires2_res = res_list2[i].xtra['ires']
        iatom_self1_res = res_list1[i].xtra['iatom_self']
        iatom_self2_res = res_list2[i].xtra['iatom_self']
        iatom1_res = res_list1[i].xtra['iatom']
        iatom2_res = res_list2[i].xtra['iatom']

        # Gather all information and add it to the return string:
        print_list = [pair_name,
                      prot_len,
                      AA_wt,
                      AA_mut,
                      ss_wt,
                      ss_mut,
                      ASA_wt,
                      ASA_mut,
                      RASA_wt,
                      RASA_mut,
                      ASA_trim_wt,
                      ASA_trim_mut,
                      RASA_trim_wt,
                      RASA_trim_mut,
                      _3D_dist,
                      _1D_dist,
                      N_dist,
                      C_dist,
                      dphi,
                      dpsi,
                      ichain_numb_wt,
                      ichain_numb_mut,
                      idist_wt,
                      idist_mut,
                      ires_wt,
                      ires_mut,
                      iatom_self_wt,
                      iatom_self_mut,
                      iatom_wt,
                      iatom_mut,
                      AA_res,
                      ss1_res,
                      ss2_res,
                      ASA1_res,
                      ASA2_res,
                      RASA1_res,
                      RASA2_res,
                      ASA1_trim_res,
                      ASA2_trim_res,
                      RASA1_trim_res,
                      RASA2_trim_res,
                      bfactor_back1_res,
                      bfactor_back2_res,
                      bfactor_side1_res,
                      bfactor_side2_res,
                      next2cut,
                      idist1_res,
                      idist2_res,
                      ires1_res,
                      ires2_res,
                      iatom_self1_res,
                      iatom_self2_res,
                      iatom1_res,
                      iatom2_res]
        print_list = [str(el) for el in print_list]
        # print(*print_list, sep='\t', end='\n', file=fh_out)
        print_str += '\t'.join(print_list) + '\n'
    return(print_str)


def cut_bad_crystal(res_list1, res_list2, mut_pos_crystal, single_pair):
    '''
    Find missing residues in the trimmed crystal structure
    and discard them if these missing residues are close to the mutation
    or opens a big gap in the structure. Also trim away short stretches of terminal residues
    when a missing residue is found close to the end of the protein sequence.
    '''
    # Define the distance cutoffs:
    max_dist = 12   # Equivalent to 3 stretched out residues (max 3 residue gaps allowed not considering proline)
    cut_dist = 4    # This distance indicates a cut in the sequence (missing residue in the crystal structure)
    min_dist = 2.5  # Most will not go below 3.5, but Proline probably will and therefore as low as 2.5 is allowed
    min_dist_from_defect = 10  # Minimun allow distance distance to a crystal defect
    min_frag_len = 20          # Mininum length between termina and cut before terminal fragment is removed

    pair_name = pair_name = '-'.join(single_pair)
    Calpha = Calpha_from_residue_list(res_list1)  # Use C-alpha to measure distance
    n_res = len(Calpha)
    cut_offset = 0  # Offset to take into account the truncation of the residue list when cutting the N-terminal off
    for i in range(2, n_res):
        i -= cut_offset
        res_list1[i].xtra['next2cut'] = 0  # Start by assuming no cut
        res_list2[i].xtra['next2cut'] = 0  # Start by assuming no cut
        dist_3D = Calpha[i] - Calpha[i - 1]
        if dist_3D > cut_dist and dist_3D < max_dist:
            res_list1[i - 1].xtra['next2cut'] = 1  # Indicate a cut
            res_list1[i].xtra['next2cut'] = 1      # Indicate a cut
            res_list2[i - 1].xtra['next2cut'] = 1  # Indicate a cut
            res_list2[i].xtra['next2cut'] = 1      # Indicate a cut
            print('Missing residue detected between residue', i - 1, 'and', i, 'in pair', pair_name)
            mut_dist1 = Calpha[i - 1] - Calpha[mut_pos_crystal]
            mut_dist2 = Calpha[i] - Calpha[mut_pos_crystal]
            if mut_dist1 < min_dist_from_defect or mut_dist2 < min_dist_from_defect:
                print('The missing residue is too close to the mutating residue. The pair is discarded:', pair_name)
                return('bad_crystal', 'bad_crystal')
            # Now check if the break creates a very short fragment, if so delete this fragment:
            if i < min_frag_len:
                # Trim away the small N-terminal fragment:
                res_list1 = res_list1[i:]
                res_list2 = res_list2[i:]
                # Check that the mutation is not inside the discarded fragment:
                if mut_pos_crystal < i:
                    print('The mutation is observed inside a short fragment broken by a missing residue. The pair is discarded:', pair_name)
                    return('bad_crystal', 'bad_crystal')
                cut_offset = i  # Add an offset to account for the truncation of the N-terminal
            elif (n_res - (i - 1)) < min_frag_len:
                # Trim away the small C-terminal fragment:
                res_list1 = res_list1[0:i]
                res_list2 = res_list2[0:i]
                # Check that the mutation is not inside the discarded fragment:
                if mut_pos_crystal > i:
                    print('The mutation is observed inside a short fragment broken by a missing residue. The pair is discarded:', pair_name)
                    return('bad_crystal', 'bad_crystal')
                break  # No reason to ontinue loop because the remainder C-terminal have been cut away
        elif dist_3D < min_dist:
            print('The distance between adjacent C-alpha atoms is suspiciously low (<{:.1f}A). The pair is discarded: {}, with a distance of {:.1f} between residue {} and {}'.format(min_dist, pair_name, dist_3D, i - 1, i))
            return('bad_crystal', 'bad_crystal')
        elif dist_3D > max_dist:
            print('The distance between adjacent C-alpha atoms is suspiciously high (>{:.1f}A). This probably means that there is a large +3 residue gap in the crystal structure. The pair is discarded: {}, with a distance of {:.1f} between residue {} and {}'.format(max_dist, pair_name, dist_3D, i - 1, i))
            return('bad_crystal', 'bad_crystal')

    return(res_list1, res_list2)


# To delete:
def dist_between_Calpha(res_list):
    '''
    Find the distance between C-alpha and return it as a list.
    This function is unused and deprecated.
    '''
    # Add linear (sequence based) distance
    dist_list = list()
    Calpha = Calpha_from_residue_list(res_list)
    for i in range(1, len(Calpha)):
        dist_3D = Calpha[i] - Calpha[i - 1]
        dist_list.append(dist_3D)
    return(dist_list)


def remove_homodimers_in_pairs(pairs):
    '''
    Takes the list of pair tuples and removed mulitple chains from the same PDB file.
    '''
    new_pairs = list()  # New list to the pair tuples with max on PDB entry
    for pair in pairs:
        p1l = pair[0].split('-')  # Split to a list of chains
        p2l = pair[1].split('-')  # Split to a list of chains
        p1s = pair[0]  # For search also keep the string of chains
        p2s = pair[1]  # For search also keep the string of chains

        # First step remove those homo dimeric chains with the same sequence
        # appearing only in one of the pair partners.
        # Remove duplicate PDB entries in pair partner 1:
        p1d = {pp[0:-1]: pp[-1] for pp in p1l}
        # Remove duplicate PDB entries in pair partner 2:
        p2d = {pp[0:-1]: pp[-1] for pp in p2l}

        # Use list mutation instead of assigning a new list:
        p1l[:] = [p + c for p, c in p1d.items()]
        p2l[:] = [p + c for p, c in p2d.items()]

        # Now remove those PDB entries that are the same for both
        # pair partner 1 and pair partner 2:
        p1l_iter = p1l[:]  # Copy for iterations
        for p1 in p1l_iter:
            if p1[0:-1] in p1s and p1[0:-1] in p2s:
                # Delete from the longest pair partner list:
                if len(p1l) > len(p2l):
                    p1l.remove(p1)  # Because p1 is used for iteration this can be found easily
                else:
                    # Notice here that it is a little more difficult to find the pair partner
                    # do delete because the chain name is different:
                    p2l = [p2 for p2 in p2l if p1[0:-1] != p2[0:-1]]

        # Check that we still have a pair, if yes then update it:
        if len(p1l) > 0 and len(p2l) > 0:
            new_pair = ('-'.join(p1l), '-'.join(p2l), pair[2])
            new_pairs.append(new_pair)
        # If no pair no nothing
        else:
            pass

    return(new_pairs)


def pair_stat(pairs):
    '''
    Small stats function printing to STDOUT.
    '''
    high = 0
    high_comp = 0
    total_comp = 0
    total_chains = 0
    numb_rep = 0
    for pair in pairs:
        p1l = pair[0].split('-')
        p2l = pair[1].split('-')
        if len(p1l) > high:
            high = len(p1l)
        if len(p2l) > high:
            high = len(p2l)

        total_chains += len(p1l) + len(p2l)
        comp = len(p1l) * len(p2l)
        total_comp += comp
        if comp > high_comp:
            high_comp = comp

        if len(p1l) > 1:
            numb_rep += len(p1l) - 1
        if len(p2l) > 1:
            numb_rep += len(p2l) - 1

    print('Highest number of pair partners i.e. same sequence chains from different PDB files with a mutant partner:', high)
    print('Highest number of pairwise comparisons necessary for a given pair in the pair list:', high_comp)
    print('Total number of chains summed over all the pairs:', total_chains)
    print('Total number of comparisons between chains necessary for all pairs:', total_comp)
    print('Number of replicate pair partners i.e. chains that can be used for determining a background difference:', numb_rep)


def get_interactions(res_list, pdb_obj, sele_chain):
    '''
    Look for possible protein-protein interactions between any residue in a res_list
    and other residues on other chains, but in the same PDB file as the input res_list originates from.
    Interactions are only found based on distance cutoff no chemical interactions
    are taken into consideration when searching.
    '''
    # res.xtra['ichain_numb'] = int
    # res.xtra['idist'] = float / NA
    # res.xtra['ichain'] = chain identifier / NA
    # res.xtra['ires'] = residue identifier / NA
    # res.xtra['iatom_self'] = atom type / NA
    # res.xtra['iatom'] = atom type of partner atom / NA
    n_chains = 0
    for chain in pdb_obj[0]:
        # Skip the self chain:
        if chain.get_id() == sele_chain:
            continue
        else:
            n_chains += 1

        dist_cut = 20  # Interactions cannot be further away than 20A backbone atom to backbone atom
        # R -> E = 7.1A + 4.9A = 12A + bond max 15A to give a bit of freedom go to 20A
        min_idist = 999  # Set this high to pick up first idist as the min_idist
        for res in res_list:
            # Find a backbone atom that can be used to find the distance.
            # CA is prefered:
            if 'CA' in res:
                res_comp = 'CA'
            elif 'N' in res:
                res_comp = 'N'
            elif 'C' in res:
                res_comp = 'C'
            # If no backbone atoms are found no interactions are annotated:
            else:
                continue

            for dres in chain:  # All residues on the chain for comparison
                if res_comp in dres:  # All residues in the chosen chain
                    idist = res[res_comp] - dres[res_comp]
                    # Only when the backbone distance is smaller than the cutoff,
                    # continue and find the distance between all atoms:
                    if idist < dist_cut:
                        for atom in res:  # All atoms for the residue in the chain for comparison
                            for datom in dres:  # All atoms for the residue in the chosen chain
                                idist = atom - datom
                                # Update the residue object with the smallest distance:
                                if idist <= min_idist:
                                    res.xtra['idist'] = idist
                                    res.xtra['ichain'] = chain.get_id()
                                    res.xtra['ires'] = dres.get_resname()
                                    res.xtra['iatom_self'] = atom.get_name()
                                    res.xtra['iatom'] = datom.get_name()
                                    min_idist = idist

                # If the no backbone atom for comparison skip that residue:
                else:
                    continue

    # Add the number of different chains observed in the crystal lastly
    # and fill out missing fields:
    for res in res_list:
        res.xtra['ichain_numb'] = n_chains
        # Fill out the xtra fields that have not been assigned anything with NA.
        # Notice that this cannot be moved to the loop above because single chain PDB files would fail:
        if 'idist' not in res.xtra:
            res.xtra['idist'] = 'NA'
        if 'ichain' not in res.xtra:
            res.xtra['ichain'] = 'NA'
        if 'ires' not in res.xtra:
            res.xtra['ires'] = 'NA'
        if 'iatom_self' not in res.xtra:
            res.xtra['iatom_self'] = 'NA'
        if 'iatom' not in res.xtra:
            res.xtra['iatom'] = 'NA'

    return(res_list)


def print_header(result_file):
    '''
    Creates a results file with a header as the first line, ready for appending results.
    '''
    with open(result_file, 'w') as fh_out:
        # Print a header:
        print_list = ['pair_name',
                      'prot_len',
                      'AA_wt',
                      'AA_mut',
                      'ss_wt',
                      'ss_mut',
                      'ASA_wt',
                      'ASA_mut',
                      'RASA_wt',
                      'RASA_mut',
                      'ASA_trim_wt',
                      'ASA_trim_mut',
                      'RASA_trim_wt',
                      'RASA_trim_mut',
                      '_3D_dist',
                      '_1D_dist',
                      'N_dist',
                      'C_dist',
                      'dphi',
                      'dpsi',
                      'ichain_numb_wt',
                      'ichain_numb_mut',
                      'idist_wt',
                      'idist_mut',
                      'ires_wt',
                      'ires_mut',
                      'iatom_self_wt',
                      'iatom_self_mut',
                      'iatom_wt',
                      'iatom_mut',
                      'AA_res',
                      'ss1_res',
                      'ss2_res',
                      'ASA1_res',
                      'ASA2_res',
                      'RASA1_res',
                      'RASA2_res',
                      'ASA1_trim_res',
                      'ASA2_trim_res',
                      'RASA1_trim_res',
                      'RASA2_trim_res',
                      'bfactor_back1_res',
                      'bfactor_back2_res',
                      'bfactor_side1_res',
                      'bfactor_side2_res',
                      'next2cut',
                      'idist1_res',
                      'idist2_res',
                      'ires1_res',
                      'ires2_res',
                      'iatom_self1_res',
                      'iatom_self2_res',
                      'iatom1_res',
                      'iatom2_res']
        print(*print_list, sep='\t', end='\n', file=fh_out)


def mp_worker(pair_info):
    '''
    Worker function that wraps all the steps in creating the outputted result.
    The function returns the total output of a pair tuple as a string.
    The worker can be called by a handler and work in parallel on multiple pairs
    on multiple cores.
    NOTICE that if a pair tuple contains many pair partners and the number
    of pairwise comparisons are high, then this can be a significant memory
    problem.
    '''
    print_str = ''  # String to store all the output from one pair tuple
    # Get all information needed by unpacking the pair_info:
    pair_number, pair_tuple, ss_dis_dict, scratch_dir, pdb_folder = pair_info
    pdbs_tuple = (tuple(pair_tuple[0].split('-')), tuple(pair_tuple[1].split('-')))
    mut_pos_seq = pair_tuple[2][0]
    # Create a tmp folder for this pair and move the PDB files needed:
    pair_folder = create_pair_folder(scratch_dir, pair_number, pair_tuple, pdb_folder)
    # If something went wrong, e.g. the file does not exists:
    if not pair_folder:
        print('No pair folder, something went wrong', pdbs_tuple)
        return()
    os.chdir(pair_folder)

    # Loop through the all pairs:
    for i in range(len(pdbs_tuple[0])):
        for j in range(len(pdbs_tuple[1])):
            # Define the single pair to look at:
            single_pair = (pdbs_tuple[0][i], pdbs_tuple[1][j])
            # PDBParser is part of Bio.PDB:
            parser = PDBParser()
            # Silence the error messages:
            with stdchannel_redirected(sys.stderr, os.devnull):
                # Parse the pdbs into BioPDB objects:
                pdb_obj1 = parser.get_structure(single_pair[0], single_pair[0][0:-1])
                pdb_obj2 = parser.get_structure(single_pair[1], single_pair[1][0:-1])
                # Add DSSP information on the full PDB file:
                pdb_obj1, pdb_obj2 = add_dssp_to_pdb_obj(pair_folder, single_pair, pdb_obj1, pdb_obj2)
                # Extract the list of residue objects:
                id1 = pdb_obj1.get_id()
                id2 = pdb_obj2.get_id()
                chain1 = id1[-1]
                chain2 = id2[-1]
                res_list1 = list(pdb_obj1[0][chain1].get_residues())
                res_list2 = list(pdb_obj2[0][chain2].get_residues())

            # Start getting possible protein-protein interactions
            # on the untouched pdb object and residue list:
            try:
                res_list1 = get_interactions(res_list1, pdb_obj1, chain1)
                res_list2 = get_interactions(res_list2, pdb_obj2, chain2)
            except:
                print(single_pair, 'get_interactions')
                print(e)
                continue

            # Correct for missing residues:
            try:
                res_list1, res_list2, mut_pos_crystal = correct_pdb_obj_pairs(res_list1, res_list2, mut_pos_seq, ss_dis_dict, id1, id2)
                # If the mutation was in the missing residues:
                if 'no_mut' in [res_list1, res_list2, mut_pos_crystal]:
                    print('Mutation was in the missing residues', single_pair, mut_pos_seq)
                    continue
            except Exception as e:
                print(single_pair, 'correct_pdb_obj_pairs')
                print(e)
                continue

            # Cut out the crystal defects:
            res_list1, res_list2 = cut_bad_crystal(res_list1, res_list2, mut_pos_crystal, single_pair)
            if 'bad_crystal' in [res_list1, res_list2, mut_pos_crystal]:
                print('Bad crystal', single_pair, mut_pos_seq)
                continue

            # Now add DSSP information to the pruned residue list:
            try:
                res_list1, res_list2 = add_dssp_to_reslist(pair_folder, single_pair, res_list1, res_list2, parser)
            except Exception as e:
                print(single_pair, 'add_dssp_to_reslist')
                print(e)
                continue

            # Find difference in torsions between the pair chains:
            try:
                torsion_diff = align_and_find_torsions(res_list1, res_list2)
            except Exception as e:
                print(single_pair, 'align_and_find_torsions')
                print(e)
                continue

            # Generate information about the distance to the variant:
            try:
                mut_dist = dist_from_mut(res_list1, mut_pos_crystal)
            except Exception as e:
                print(single_pair, 'dist_from_mut')
                print(e)
                continue

            # Add all the results up and concatenate them with the previous results:
            try:
                # Before adding anything check if the res_list is still long enough:
                if len(res_list1) > 60 or len(res_list2) > 60:
                    result = print_res(single_pair, mut_dist, torsion_diff, res_list1, res_list2, mut_pos_crystal)
                    if result:
                        print_str += result
                else:
                    print('Pair deleted because it was too short just before printing:', single_pair)
            except Exception as e:
                print(single_pair, 'print_res')
                print(e)
                continue

    # Clean-up:
    shutil.rmtree(pair_folder)
    return(print_str)


def mp_handler(pairs, ss_dis_dict, scratch_dir, pdb_folder, result_file, np):
    '''
    Multi process handler. Start a pool of processes and queues them
    on a number of cores. Runs the mp_worker and prints the output sequentially.
    '''
    # Generate a list of info that is needed for each pair
    # to run indepedently by the worker process:
    pair_info_list = list()
    for pair_number, pair_tuple in enumerate(pairs):
        pair_info_list.append([pair_number, pair_tuple, ss_dis_dict, scratch_dir, pdb_folder])

    # Start the pool with X cores:
    pool = multiprocessing.Pool(np)
    # Print output header:
    print_header(result_file)
    # Continue printing the results from each process in the pool:
    with open(result_file, 'a') as fh:
        for result in pool.imap_unordered(mp_worker, pair_info_list):
            if result:
                fh.write(result)
    pool.close()
    pool.join()


if __name__ == "__main__":
    # Define input variables:
    res = 2.5           # Use minimum 2.5A resolution
    min_seq_len = 60    # Minimum sequence length to avoid small unstructured proteins
    pair_distance = 1   # The hamming distance of the pairs to mine for

    # Get the validated PDB IDs:
    QC_PDB_dict = valid_PDB_dict(args.pdb_folder, res, args.cache_dir)
    # print(QC_PDB_dict['1ubi'])
    # print(QC_PDB_dict['5a1a'])

    # Read the missing residue and secondary structure PDB information:
    ss_dis_dict = parse_ss_dis(args.ss_dis_file, args.cache_dir)
    # print(len(ss_dis_dict['4bjsA']['disorder']))
    # print(len(ss_dis_dict['4bjsA']['sequence']))
    # print(len(ss_dis_dict['4bjsA']['secstr']))

    # Read the BioLIP database:
    biolip_dict = parse_biolip(args.biolip_fnam, args.cache_dir)

    # Split the PDB chains into bins according to the length of their amino acid sequence:
    seq_liglen_dict = seq_liglen_bin(min_seq_len, QC_PDB_dict, ss_dis_dict, biolip_dict, args.cache_dir)

    ### Example code to extract pairs that does not have any ligands:
    # seq_len_dict = seq_liglen_dict['no_lig']
    # seq_len_dict = dedup_seq_len_dict(seq_len_dict)  # Currently deprecated since this function is moved into the seq_liglen_bin function
    # pairs1, distance_distribution = find_chain_pairs(pair_distance, seq_liglen_dict, args.cache_dir)
    ##########

    # Find the pairs of PDB IDs that deviate by a single amino acid:
    # seq_liglen_dict = dedup_seq_liglen_dict(seq_liglen_dict)  # Currently deprecated since this function is moved into the seq_liglen_bin function
    pairs1, distance_distribution = find_chain_pairs2(pair_distance, seq_liglen_dict, args.cache_dir)
    pairs1 = remove_homodimers_in_pairs(pairs1)

    # Do all the calculation in a parallel pool and print the results:
    mp_handler(pairs1, ss_dis_dict, args.scratch_dir, args.pdb_folder, args.result_file, args.np)
