from typing import List, Dict
from path import Path
from io import StringIO

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from epitope_mutations import Epitope


DEFAULT_COLUMN_NAME_DICT_CD8 = {
        "epitope_id_col_name" : 'Epitope ID', 
        "sequence_col_name" : 'Sequence',
        "protein_col_name": 'Protein',
        "start_col_name" : 'Mapped Start Position',
        "end_col_name" : 'Mapped End Position',
        "HLA_restrictions_col_name" : 'MHC Restriction',
        "length_col_name" : 'length'}
DEFAULT_OUTPUT_COLUMN_ORDER = [
      'epitope_id',
      'protein',
      'start',
      'end',
      'length', 
      'original_sequence',
      'modified not final sequence',
      'mutated_sequence',
      'contains deletion', 
      'contains insertion',
      'mutations in this epitope',
      'HLA restrictions']

# TODO: The list of necessary keys needs some thinking. At the moment, when generating from .csv, it will just fill a column with None.
# And that is fine with the function below. That is likely not the way it should work.
def epis_from_dicts(dict_list:List[dict]) -> List[Epitope]:
    epilist = []
    for d in dict_list:
        keys = d.keys()
        assert all([
            'epitope_id' in keys,
            'protein' in keys,
            'sequence' in keys,
            'start' in keys,
            'end' in keys
            ]), 'Dictionary to create Epitope does not contain all necessary keywords.'
        epilist.append(Epitope(d))
    return epilist


def epilist_from_csv(
        filepath:Path, column_name_dict:dict=DEFAULT_COLUMN_NAME_DICT_CD8, delim:str=';') -> List[Epitope]:
    
    cnd = column_name_dict
    df = pd.read_csv(filepath, sep=delim)
    # print (df.columns)
    for column_name in column_name_dict.values():
        try:
            df[column_name]
        except:
            print(f'Data input error: assumed column name {column_name} does not match a column in input file.')
            df[column_name] = None
            
    # Need to get rid of lines, which contain no information or not sufficient information.
    df.dropna(inplace=True, subset=[cnd['sequence_col_name'],cnd['start_col_name']])
    #This is just a fix for the moment. Those should not be in the input data to begin with, but it's easy to fix here.
    df = df[df[cnd['start_col_name']] != 'mut']
    try:
        # TODO: Seems this is not actually working, but doesn't seem to matter for the moment.
        df.astype({cnd['start_col_name']:'int32', cnd['length_col_name']:'int32'}, copy=None)
    except ValueError:
        print('Start or length column cannot be converted into integer.')
    if not cnd['end_col_name'] in df.columns:
        df[cnd['end_col_name']] = df[cnd['start_col_name']]+df[cnd['length_col_name']] -1
    #return df
    df.rename(inplace=True, columns={
        cnd["epitope_id_col_name"]:'epitope_id',
        cnd["sequence_col_name"]:'sequence',
        cnd["protein_col_name"]: 'protein',
        cnd["start_col_name"]: 'start',
        cnd["end_col_name"]: 'end',
        cnd["length_col_name"]: 'length',
        cnd["HLA_restrictions_col_name"]: 'HLA_restrictions'
    })
    # print (df.dtypes)
    list_of_dicts = df.to_dict(orient='records')
    epilist = epis_from_dicts(list_of_dicts)
    return epilist


def mutationlist_from_csv(
        filepath:Path,
        delimiter:str=';',
        protein_col_name:str='protein',
        position_col_name:str='Position',
        mutation_col_name:str='mutation'
        )->List[dict]:
    df = pd.read_csv(filepath,sep=delimiter)
    for column_name in [protein_col_name,position_col_name, mutation_col_name]:
        try:
            df[column_name]
        except KeyError:
            print(f'Data input error: assumed column name {column_name} does not match a column in input file.')
            #df[column_name] = None
    df.rename(inplace=True, columns={
        protein_col_name:'protein',
        position_col_name:'position',
        mutation_col_name:'new'
    })
    df.replace('del','-', inplace=True)
    return df.to_dict(orient='records')

def read_sequences_from_fasta(fastafilepath:Path)->Dict[str,SeqRecord]:
    # Actually, the fasta format requires the sequence identifier in the line starting with '>'
    # to be without spaces. The SeqIO parser does cut off after the first whitespace.
    # Users might not adhere to this policy, so here is trying to prevent that error 
    # by replacing whitespaces (except trailing ones) with '_'.
    with open(fastafilepath,'r') as filehandle:
        all_lines = []
        for line in filehandle:
            if line.startswith('>'):
                mod_line = line.rstrip().replace(' ','_')
                all_lines.append(mod_line)
            else:
                all_lines.append(line)
    # Now generate a virtual file handle from that list of lines with line breaks and feed to biopython.
    filehandle = StringIO('\n'.join(all_lines))
    seq_gen = SeqIO.parse(filehandle,'fasta')
    seq_dict = SeqIO.to_dict(seq_gen)
    return seq_dict


def generate_mutated_sequences(original_sequence_dict:Dict[str,SeqRecord], mutations:list) -> Dict[str,str]:
    # Since methods to apply mutations are anyways already in the Epitope class, using this.
    pseudo_epis = []
    for key, seqcord in original_sequence_dict.items():
        pseudo_epis.append(Epitope({
            'epitope_id' : key,
            'protein' : key,
            'sequence' : str(seqcord.seq),
            'start' : 1,
            'end' : len(str(seqcord.seq)) +1,
            'HLA_restrictions' : ''

        }))
    output_dict = {}
    for pseudoepi in pseudo_epis:
        pseudoepi.apply_mutations(mutations)
        output_dict[pseudoepi.protein] = pseudoepi.mod_sequence.replace('-','')
    return output_dict

def reorder_dataframe_columns(df:pd.DataFrame,column_order:List[str]=DEFAULT_OUTPUT_COLUMN_ORDER)->pd.DataFrame:
    # Some columns might not exist like 'protein'. So I can't just apply the list as new order, but have to keep just those, which are actually in the dataset.
    new_columns = []
    for column in column_order:
        if column in df.columns:
            new_columns.append(column)
    df_out = df[new_columns]
    return df_out

def generate_stats_epi_mutations(df:pd.DataFrame) -> pd.DataFrame:
    no_unique_epis = df.shape[0]
    no_mutations = df[df['mutations in this epitope']>0].shape[0]
    no_deletions = df[df['contains deletion'] == True].shape[0]
    no_insertions = df[df['contains insertion'] == True].shape[0]
    stat_dict = {
        'Type':['All epitopes','With any mutations','With deletions','With insertions'],
        'Number': [no_unique_epis,no_mutations, no_deletions, no_insertions],
        'Of total': [no_unique_epis/no_unique_epis, no_mutations/no_unique_epis, no_deletions/no_unique_epis, no_insertions/no_unique_epis]

        }
    stats = pd.DataFrame(stat_dict)
    return stats

def generate_unique_epitope_df(non_unique_df:pd.DataFrame)->pd.DataFrame:
    # Since epitope ID can be just NaN for the whole column, this needs to be dropped, bc otherwise grouping doesn't work.
    non_unique_df_non_na = non_unique_df.dropna(axis=1)
    # Then let's find all columns existing except the HLA restrictions, which I want to aggregate.
    identical_columns = [col for col in non_unique_df_non_na.columns if not col == 'HLA restrictions']
    unique_df = non_unique_df_non_na.groupby(identical_columns, as_index=False).agg({'HLA restrictions':','.join})
    return unique_df