
import pandas as pd
from Bio.Seq import Seq
from Bio import Align
import os
from tqdm import tqdm
import argparse

# # Load AIRR Data and Resources


def get_script_dir():
    """Returns the directory where the script is located."""
    return os.path.dirname(os.path.abspath(__file__))

def chain_dict_maker(txt_file):
    """Parses a text file into a dictionary of chain sequences."""
    list_of_chains = txt_file.split('>')
    chains_dict = {}
    
    for chain in range(1, len(list_of_chains)):
        chain_data = list_of_chains[chain].replace(' ', '').split('|')
        chain_data[-1] = chain_data[-1].replace('\n', '')
        AAseq = chain_data[-1]
        chain_name = chain_data[1]
        chains_dict[chain_name] = AAseq

    return chains_dict

def file_loader(directory):
    """Loads text files from a directory and processes them into a dictionary."""
    if not os.path.exists(directory):
        print(f"Error: The directory '{directory}' does not exist. Please provide the correct path.")
        return {}

    files_content = {}
    master_dict = {}

    # Iterate through all files in the directory
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)

        if filename.endswith('.txt'):
            with open(file_path, 'r', encoding='utf-8') as file:
                content = file.read()
            files_content[filename] = content  # Store file content with filename as key

    for key, txt_file in files_content.items():
        new_key = key[:-4]  # Remove '.txt' extension
        chain_dict = chain_dict_maker(txt_file)
        master_dict[new_key] = chain_dict

    return master_dict

# Argument parser setup
parser = argparse.ArgumentParser(description="Process AIRR file and generate full-length, endogenous codon TCR nucleotide sequences.")
parser.add_argument("--airr_path", "-a", type=str, help="Path to the AIRR file (TSV)")
parser.add_argument("--libraries_path", "-l", type=str, default=None, help="Path to the gene libraries directory (optional)")
parser.add_argument("--output_path", "-o", type=str, default=os.getcwd(), help="Path to save the output CSV file (default: current directory)")

args = parser.parse_args()

# Load the AIRR file
airr_path = args.airr_path

try:
    airr = pd.read_csv(airr_path, sep='\t')
except FileNotFoundError:
    print(f"Error: The file at {airr_path} does not exist.")
    exit(1)
except Exception as e:
    print(f"An error occurred while reading the AIRR file: {e}")
    exit(1)

# Determine the library path
if args.libraries_path:
    libraries_path = args.libraries_path
else:
    # Default to the 'data' directory next to the script
    libraries_path = os.path.join(get_script_dir(), "v_genes")

# Load the gene libraries into master_dict
master_dict = file_loader(libraries_path)


# # DEFINE FUNCTIONS


def get_seqs():

    #CDR3;
    seq1 = 'KYIHWYRQLPSQGPEYVIHGLTSNVNNRMASLAIAEDRKSSTLILHRATLRDAAVYYCILVNNNAGNMLTFGGGTRLMVKPHIQNPDPAVYQLRDSKSSDKSVCLFTDFD'
    cdr = seq1
    #seq_1 = 'gTATAGGGACAGGAAAGAAGATCACTCTGGAATGTTCTCAAACCATGGGCCATGACAAAATGTACTGGTATCAACAAGATACAGGAATGGAACTACAGCTCATCCACTATTCCTATGGAGTTAATTGCACAGAGAAGGGAGATCTTTCCTCTGAGTCAACAGTCTCCAGAATAAGGACGGAGCATTTTCCCCTGACCCTGGAGTCTGCCAGGCCCTCACATACCTCTCAGTACCTCTGTGCCAGCAGGTTAGCGGGAGCAGTGGATGAGCAGTTCTTCGGGCCAGGGACACGGCTCACCGTGCTAGAGGACCTGAAAAACGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATC'
    #seq1 = seq_1.lower()

    #TRAV;
    seq2 = 'MKLVTSITVLLSLGIMGDAKTTQPNSMESNEEEPVHLPCNHSTISGTDYIHWYRQLPSQGPEYVIHGLTSNVNNRMASLAIAEDRKSSTLILHRATLRDAAVYYCILRD'
    trv = seq2
    #TRAJ;
    seq3 = 'NNNAGNMLTFGGGTRLMVKP'
    trj = seq3
    #TRAC;
    seq4 = 'XIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESSCDVKLVEKSFETDTNLNFQNLSVIGFRILLLKVAGFNLLMTLRLWSS'
    trc = seq4
    #seq4:
    #seq_2 = 'atgactatcaggctcctctgctacatgggcttttattttctgggggcaggcctcatggaagctgacatctaccagaccccaagataccttgttatagggacaggaaagaagatcactctggaatgttctcaaaccatgggccatgacaaaatgtactggtatcaacaagatccaggaatggaactacacctcatccactattcctatggagttaattccacagagaagggagatctttcctctgagtcaacagtctccagaataaggacggagcattttcccctgaccctggagtctgccaggccctcacatacctctcagtacctctgtgccagcagtgaata'
    #seq2 = seq_2.lower()


    cdr = 'GAAGGGTACAAAGTCTCTCGAAAAGAGAAGAGGAATTTCCCCCTGATCCTGGAGTCGCCCAGCCCCAACCAGACCTCTCTGTACTTCTGTGCCAGCAGTTTATCGACAAGCTCGACCGGGGCCAACGTCCTGACTTTCGGGGCCGGCAGCAGGCTGACCGTGCTGGAGGACCTGAAAAACGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATC'
    trv = 'atgggcccccagctccttggctatgtggtcctttgccttctaggagcaggccccctggaagcccaagtgacccagaacccaagatacctcatcacagtgactggaaagaagttaacagtgacttgttctcagaatatgaaccatgagtatatgtcctggtatcgacaagacccagggctgggcttaaggcagatctactattcaatgaatgttgaggtgactgataagggagatgttcctgaagggtacaaagtctctcgaaaagagaagaggaatttccccctgatcctggagtcgcccagccccaaccagacctctctgtacttctgtgccagcagtttatc'
    trc = 'GAGGACCTGAAAAACGTGTTCC'



    return(cdr.upper(), trv.upper(), trc.upper())

def get_aligner():
    aligner = Align.PairwiseAligner()

    aligner.mode = 'global'
    # Set the scoring scheme
    aligner.open_gap_score = -10  # Penalty for opening a gap
    aligner.extend_gap_score = -1  # Penalty for extending a gap

    # Set the match and mismatch scores
    aligner.match_score = 1
    aligner.mismatch_score = -0.5

    # Prevent end gaps by forcing end gap extensions to be zero
    aligner.target_end_gap_score = 0
    aligner.query_end_gap_score = 0

    return(aligner)


# PANDAS-based consenesusizer (numpy was not faster)
# x needs to be a dataframe with one the column on the left containing consensus counts and the column on the right containing sequences
# note that consensus counts and operations performed on them are not useful for the current build and need to be cleaned out
# method is 'minimum' or 'padded'
    # minimum will trim characters from the left side of all sequences to the length of whichever is the shortest
    # padded will add 'N' characters to the left side of all sequences to make them the same length as the longest
# pad_char = 'N' or '-'
def consensusizer(x, method = 'padded', confidence = 0.8):
    pad_char = '-'
    # Pass in a df where the left column has consensus counts and the right column has sequences
    x.columns = ['consensus_count','sequence']
    sequences_df = x[['consensus_count','sequence']].copy()
    sequences_df = sequences_df[sequences_df['sequence'].str.contains('Error') == False].copy() # Remove sequences that contain string: "Error"
    sequences_df['SRP_seq'] = sequences_df['sequence']
    sequences_df['SRP_len'] = sequences_df['sequence'].apply(lambda x: len(x))
    minimum_len = sequences_df['SRP_len'].min()
    maximum_len = sequences_df['SRP_len'].max()
    sequences_df['SRP_min_seq'] = sequences_df['SRP_seq'].apply(lambda x: x[-minimum_len:])

    # Add '-' characters to pad length on the 5' so all sequences are the same length as the longest one.
    sequences_df['SRP_seq_padded'] = sequences_df['SRP_seq'].apply(lambda x: (maximum_len - len(x))*pad_char+x)

    consensus_string = ''
    
    if method == 'padded':
        working_column = 'SRP_seq_padded'
        working_length = maximum_len
    else:
        working_column = 'SRP_min_seq'
        working_length = minimum_len

    #debugger:
    #minimum_len = 2
    for i in range(0, working_length):
        # Create a "position" column (pos) and update with the nucleotides in each sequence at that postion
        sequences_df['pos'] = sequences_df[working_column].apply(lambda x: x[i]) 
        sequences_df['pos_consensus'] = sequences_df['consensus_count']
        # Group by nucleotides and sum the consensus count
        ###best_nuc = sequences_df.groupby('pos')['consensus_count'].sum().reset_index()
        best_nuc = sequences_df.groupby('pos').agg(
            consensus_count = ('consensus_count', 'sum'),
            sequence_count = ('pos', 'count')
        ).reset_index()
        best_nuc['prevalence_percent'] = best_nuc['sequence_count'].apply(lambda x: x/len(sequences_df))

        if confidence != 0:
            best_nuc_prevalence = best_nuc.loc[best_nuc['prevalence_percent'].idxmax(),'pos']
            best_nuc_consensus = best_nuc.loc[best_nuc['consensus_count'].idxmax(),'pos']
            if best_nuc_consensus != best_nuc_prevalence:
                # These are often not equal. Imagine in the padded case, the 5' prevalence max will be '-' even if a high consensus seq has a nuc
                pass
                #return('Error, somehow the consensus and prevalence methods are in disagreement')
            elif best_nuc.loc[best_nuc['prevalence_percent'].idxmax(),'prevalence_percent'] >= confidence:
                consensus_string += best_nuc_prevalence
            else:
                consensus_string += pad_char
        else:
            # Take nucleotide with most consensus
            best_nuc = best_nuc.loc[best_nuc['consensus_count'].idxmax(),'pos']
            consensus_string += best_nuc
    
    return(consensus_string)

# aligns every v_gene against the consensus sequence to decide which is best
# Not useful under current build, future plan is to have it be the default method when v_call is empty
    # Another fix that needs to be made: Modify alignment formatting in dictionary ouput so it is readable with unispace font (make it not a dictionary)
def v_gene_finder_concise(sequence):
    

    genes = master_dict['ALL_Vregion_nuc']

    genes = master_dict['ALL_nuc']
    
    aligner = get_aligner()

    gene_dict = {}
    alignment_dict = {}
    top_score = 0.0
    for gene_name in genes.keys():
        gene_sequence = genes[gene_name]
        alignment = aligner.align(sequence.upper(), gene_sequence.upper())[0]
        score = alignment.score
        gene_dict[gene_name] = score
        alignment_dict[gene_name] = str(alignment)
        # These lines used to return only the best scoring allele
        if score > top_score:
            top_score = score
            top_gene = gene_name
            best_align = alignment

    sorted_gene_dict = dict(sorted(gene_dict.items(), key=lambda item: item[1], reverse=True))
    top_5_alignments = {}
    count = 0
    for gene_name in sorted_gene_dict.keys():
        if count <= 5:
            top_5_alignments[gene_name] = alignment_dict[gene_name]
        else:
            break
        count += 1


    # Return a sort dict from best to worst alignment score:
    return(sorted_gene_dict,top_5_alignments)
    #Return only the top gene:
    #return(f"{top_gene}: {top_score}")


aligner = get_aligner()
alignment = aligner.align('atg','gta')[0]
#print(str(alignment))


# aligner constructor creates a long string for each given sequence with '-' padding characters where needed so that each string is the same length
# equal length strings are then stored in a df where each row is a sequence and each character is a column
# Columns are iterated over to identify a nucleotide and NOT a '-' character.
    # This strategy has a cdr3 bias. If a nucleotide is found in the cdr3 string, it will be used to construct the sequence
    # Given loss of sequencing accuracy for long sequences, there is some concern that incorrect nucleotides are chosen.
    # This is why a consensus step needs to be implemented either before or after sequence construction
    # Note that changing enumerate(cdr) to enumerate(trv) results in amino acid substitution in the cdr3. DO NOT DO THIS

def aligner_constructor_pairwise(seq0, seq1):
    #cdr, trv, trc = get_seqs()

    aligner = get_aligner()

    alignment = aligner.align(seq0,seq1)[0]
    seq0 = alignment[0]
    seq1 = alignment[1]

    aligned_df = pd.DataFrame({'first':list(seq0),
                               'second':list(seq1)
                               }).T # Put all the chains in a df and transpose so they can be indexed more cleanly

    # Iterate through every column in aligned df and find the row index that has a letter in it when the cdr string contains '-'
    # seq0 is the basis for string construction. Nucleotides present in seq0 will always be present in combined string
    combined_string = ''
    for i, letter in enumerate(seq0):
        if letter == '-' or letter == 'N': # If the cdr string says '-', look for another letter in the column
            options = aligned_df.loc[:,i] # Extract column
            try:
                new_letter = options[(options != '-') & (options != 'N')].iloc[0] # Find where in the column there's a nucleotide
            except:
                new_letter = '!'
        else:
            new_letter = letter
        #print(f'{i},{new_letter}')
        combined_string += new_letter

    return(alignment, seq0, seq1, combined_string)


def clean_airr(airr):
    airr_filt = airr[airr['sequence_alignment'].isna() == False].copy()
    airr_filt = airr_filt[airr_filt['junction'].isna() == False].copy()
    airr_filt = airr_filt[airr_filt['productive'] == True].copy()

    # This step is important. Clonotypes are defined by the "short" version of junction_aa
    # Justification is demonstrated in "temp_check" df below.
    airr_filt['junction_aa_short'] = airr_filt['junction_aa'].apply(lambda x: x[1:-1])

    airr_filt_size = airr_filt.groupby('junction_aa_short').agg(
        junction_aa_short_count = ('junction_aa_short', 'size')
    ).reset_index()
    airr_filt = airr_filt.merge(airr_filt_size)

    # These sequences taken from IMGT
    TRAC_nuc = 'natatccagaaccctgaccctgccgtgtaccagctgagagactctaaatccagtgacaag'
    TRAC_nuc = TRAC_nuc.upper()
    TRBC1_nuc = 'gaggacctgaacaaggtgttcccacccgaggtcgctgtgtttgagccatcagaagcagag'
    TRBC1_nuc = TRBC1_nuc.upper()
    TRBC2_nuc = 'gaggacctgaaaaacgtgttcccacccgaggtcgctgtgtttgagccatcagaagcagag'
    TRBC2_nuc = TRBC2_nuc.upper()

    # These are the exact sequences used to prime for C-iso PCR
    beta_c1_nuc = 'gaggacctgaacaaggtgttcc'
    beta_c1_nuc = beta_c1_nuc.upper()
    beta_c2_nuc = 'GAGGACCTGAAAAACGTGTTCC'
    alpha_c_nuc = 'ATCCAGAACCCTGACCCTGC'

    #IMGT choose
    nucleotides = [TRAC_nuc, TRBC1_nuc, TRBC2_nuc]
    #SRP choose
    nucleotides_SRP = [alpha_c_nuc, beta_c1_nuc, beta_c2_nuc]

    c_call = ['TRAC','TRBC1','TRBC2']
    c_genes = pd.DataFrame(data = {'c_call':c_call,'c_call_nucleotides':nucleotides,'c_call_nuc_SRP':nucleotides_SRP})

    airr_filt = airr_filt.merge(c_genes, on='c_call', how = 'left')
    
    airr_filt.fillna('')

    return(airr_filt)


# # Clean data


airr_filt = clean_airr(airr)


# # Begin construction


def consensus_wrapper(df):

    grouped = df.groupby('junction_aa_short')

    #df_out = pd.DataFrame(data={'junction_aa_short':df['junction_aa_short'].unique()})
    df_out = df.groupby('junction_aa_short')['junction_aa_short_count'].count().reset_index()


    df_out['sequence_consensus'] = ''

    print('CONSENSUSIZING')
    for i, (group_name, group_data) in tqdm(enumerate(grouped), total=len(grouped)):
        data = group_data[['consensus_count','sequence']]
        cons_seq = consensusizer(data, method='padded', confidence = 0.6)
        df_out.loc[df_out['junction_aa_short'] == group_name,'sequence_consensus'] = cons_seq

    return(df_out)


def v_search_wrapper(df_out):
    print('Searching for optimal V-gene')
    tqdm.pandas()  # Enable progress bar for apply
    df_out[['v_search_all','top_alignments']] = df_out['sequence_consensus'].progress_apply(v_gene_finder_concise).apply(pd.Series)
        
    df_out['v_search_top'] = df_out['v_search_all'].apply(lambda x: f"{list(x.keys())[0]}: {list(x.values())[0]}")

    return(df_out)


def v_call_consolidation(df_out, airr_filt):

    df_v_calls = airr_filt.groupby('junction_aa_short').agg( #['v_call'].apply(lambda x: x.mode()[0]).reset_index()
        v_call = ('v_call', lambda x: x.mode()[0]),
        v_call_list = ('v_call', lambda x: pd.unique(x)),
        c_call_nucleotides = ('c_call_nucleotides', lambda x: x.mode()[0])
    ).reset_index()
    step = df_out[['junction_aa_short']].merge(df_v_calls)
    df_out['v_call'] = step['v_call']
    df_out['c_call_nucleotides'] = step['c_call_nucleotides']
    #df_out['v_call_list'] = step['v_call_list']

    #df_out['v_call'] = df.groupby('junction_aa_short')['v_call'].apply(lambda x: x.value_counts().to_dict()).reset_index()['v_call']


    return(df_out)


# Performs construction. Works hard to figure out what v-gene nucleotides to use
# Most notable function inside is aligner_constructor_pairwise()
def construction_wrapper(df_out):
    
    pass_fail = 'FAIL'

    #print(f'ConstucTCR-ing: {len(df_out)} iterations will be performed')

    #df_out['warnings'] = ''
    #df_out['warned_sequences'] = ''

    new_cols = ['v_call_nuc','warnings','warned_sequences','FULLSEQ']
    for col in new_cols:
        df_out[col] = ''

    sequence = df_out['sequence_consensus'].values[0]

    v_gene = df_out['v_call'].values[0]

    try:    
        v_gene_nuc = master_dict['ALL_nuc'][v_gene]
        v_source = 'long'
    except:
        try:
            v_gene_nuc = master_dict['ALL_Vregion_nuc'][v_gene]
            v_source = 'short'
            df_out['warnings'] = df_out['warnings'] + 'Incompletely sequenced gene found upstream, construction will continue; '

        except:
            v_gene = list(v_gene_finder_concise(sequence)[0].keys())[0]
            v_gene_nuc = master_dict['ALL_nuc'][v_gene]
            # default v_source is fully sequenced genes in the concise function. Can add Vregion compatibility but would prefer not to
            v_source = 'long'

    # Store the nucleotides for whichever v_gene was decided on in a new column
    df_out['v_call_nuc'] = f"{v_gene_nuc.upper()}"

    try: 
        cdr = sequence
        trv = v_gene_nuc
        trc = df_out['c_call_nucleotides'].values[0]
        # Add trv nucleotides onto cdr string
        alignment1, cdr, trv, combined_string = aligner_constructor_pairwise(cdr.upper(),trv.upper())
        # Add trc nucleotides onto the combined string from the last construction
        alignment2, combined_string, trc, final_string = aligner_constructor_pairwise(combined_string.upper(),trc.upper())

        # The "ALL_Vregion_nuc" library from IMGT has genes that are not fully sequenced and therefore are missing nucleotides at the 5'
        # Constructing TCRs with these genes (that were identified by the program that generates the AIRR file) thus leaves sequences incomplete
        # The "if" block below finds the fully-sequences gene with the best alignment to whatever was made with the incomplete gene and adds nucleotides from the fully sequenced gene to the 5'
        if v_source == 'short':
            df_out['warned_sequences'] = df_out['warned_sequences'] + final_string
            # Find the best alignment
            v_gene = list(v_gene_finder_concise(sequence)[0].keys())[0]
            # retrieve nucleotides
            v_gene_nuc = master_dict['ALL_nuc'][v_gene]
            # Construct final string again by passing final string and the best v_gene_nuc match
            alignment3, incomplete_final_string, v_gene_nuc, final_string = aligner_constructor_pairwise(final_string.upper(),v_gene_nuc.upper())
        else:
            pass

    except ValueError:
        # If something fails, you'll get drowned in prints. They will have "-" characters generated by aligner_constructor so you can see how the alignment was generated
        print(sequence)
        print(cdr)
        print(trv)
        print(trc)
        print('')

    # exclamation points are what get returned by aligner constructor when both characters in a postion are in {'-','N'}
    # In this case, the cdr, trv alignment will be returned for inspection. Example of why this happens: 
    # Two cells have the same cdr3, one of them has a SNV. During consensusizer, neither will be over the 0.5 position occurance threshold and a "-" will be returned. 
    # The trv won't code for that position therefore it becomes a "!" during construction
    if '!' in final_string:
        df_out['warnings'] = df_out['warnings'] + 'Some positions inconclusive; '
        df_out['warned_sequences'] = df_out['warned_sequences'] + final_string
        df_out['FULLSEQ'] = str(alignment1)
        # Re-perform consensusizer using consensus count as deciding factor for '!' positions
    else:
        df_out['FULLSEQ'] = final_string
        pass_fail = 'pass'

    #v_and_cdr = region_finder_pro_x()
    
    return(pass_fail, df_out)


def df_out_gen(df):

    df_out = df.groupby('junction_aa_short')['junction_aa_short_count'].count().reset_index()
    df_v_calls = airr_filt.groupby('junction_aa_short').agg(
        v_call = ('v_call', lambda x: x.mode()[0]),
        v_call_list = ('v_call', lambda x: pd.unique(x)),
        c_call_nucleotides = ('c_call_nucleotides', lambda x: x.mode()[0])
    ).reset_index()
    step = df_out[['junction_aa_short']].merge(df_v_calls)
    df_out['v_call'] = step['v_call']
    df_out['c_call_nucleotides'] = step['c_call_nucleotides']

    return (df_out)


# This function operates on each junction_aa group.
    # ouputs are saved in df_slice, to be concatonated later
    # consensus_confidence_threshold is passed to consensusizer. 
        # This value determines the minimum frequency a nucleotide must be found at a position for it to be called the consensus nucleotide
def group_pipe(df_slice, group_data, cons_conf_threshold):

    # CONSENSUSIZE
    data = group_data[['consensus_count','sequence']]
    cons_seq = consensusizer(data, method='padded', confidence = cons_conf_threshold)
    df_slice['sequence_consensus'] = cons_seq

    # V SEARCH (Not imporant)
    sorted_gene_dict,top_5_alignments = v_gene_finder_concise(cons_seq)
    df_slice['v_search_all'] = [sorted_gene_dict]
    df_slice['top_alignments'] = [top_5_alignments]

    df_slice['v_search_top'] = df_slice['v_search_all'].apply(lambda x: f"{list(x.keys())[0]}: {list(x.values())[0]}")

    # CONSTRUCT
    pass_fail, df_slice = construction_wrapper(df_slice)
    
    df_slice['consensus_confidence_threshold'] = cons_conf_threshold
    df_slice['pass_fail'] = pass_fail

    return(df_slice, pass_fail)


# df is used to make grouped data
df = airr_filt
grouped = df.groupby('junction_aa_short')

# df_out has one junction per row. Should probably be named df_unique
df_out = df_out_gen(df)

df_concatr = pd.DataFrame()

threshold_options = [0.9,0.8,0.6,0.5]

print('Iterating over groups')
for i, (group_name, group_data) in tqdm(enumerate(grouped), total=len(grouped)):
    
    # Select the single row/junction to be worked on
    df_slice = df_out[df_out['junction_aa_short'] == group_name].copy()

    # Iterate over threshold values and run group_pipe(). 
        # If construction works, the loop is broken. If not, a lower threshold is tried.
    for consensus_confidence_threshold in threshold_options:
        df_slice, pass_fail = group_pipe(df_slice, group_data, consensus_confidence_threshold)
        if pass_fail == 'pass':
            break
        else:
            #print(f'Construction failed on {group_name} junction using {consensus_confidence_threshold}, trying again')
            pass
    
    # Concat the information generated in df_slice to the ouput dataframe df_concatr
    df_concatr = pd.concat([df_concatr, df_slice], ignore_index=True)




df_out = df_concatr
df_out_filt = df_out[~df_out['FULLSEQ'].str.contains('target', na=False)].copy()
df_out_errors = df_out[df_out['FULLSEQ'].str.contains('target', na=False)].copy()

df_out_filt['FULLSEQ_AA'] = df_out_filt['FULLSEQ'].apply(lambda x: str(Seq(x).translate()))
df_out_filt = df_out_filt.sort_values(by = 'junction_aa_short_count', ascending = False).reset_index(drop=True)


# Save the output file
output_file = os.path.join(args.output_path, "construct-tcr_OUTPUT.xlsx")
df_out_filt.to_excel(output_file, index=False)
print(f"Output saved to {output_file}")


