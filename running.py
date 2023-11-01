from pathlib import Path
from typing import Dict
import logging

import pandas as pd 

from code.epitope_mutations import Epitope
from code.input_output import mutationlist_from_csv, read_sequences_from_fasta, epilist_from_csv, generate_mutated_sequences,reorder_dataframe_columns, generate_stats_epi_mutations, generate_unique_epitope_df

COLUMN_NAME_DICT_CD4 = {
        "epitope_id_col_name" : 'Epitope ID', 
        "sequence_col_name" : 'Name',
        "protein_col_name": 'Molecule Parent',
        "start_col_name" : 'starting',
        "end_col_name" : 'Mapped End Position',
        "HLA_restrictions_col_name" : 'MHC restriction',
        "length_col_name" : 'length'}

def setup_logging(level:int,output_path, name_stem) -> logging.Logger :
  l = logging.getLogger('epitope_mutations_run')
  l.setLevel(level)
  fh = logging.FileHandler(Path(output_path).joinpath(name_stem +'.log'), mode='w')
  formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
  fh.setFormatter(formatter)
  l.addHandler(fh)
  return l

def run(epitope_path:str,mutationspath:str,original_sequences_path:str, output_path:str, name_stem:str) -> Dict[str,pd.DataFrame]:
  epitope_path = Path(epitope_path).expanduser()
  mutationspath = Path(mutationspath).expanduser()
  original_sequences_path = Path(original_sequences_path).expanduser()
  output_path = Path(output_path).expanduser()
  l = setup_logging(logging.DEBUG, output_path, name_stem)
  # Now using as raw_epis here and then later filter out those, where the sequence is actually not contained in the ancestral wuhan sequence.
  raw_epis = epilist_from_csv(epitope_path, COLUMN_NAME_DICT_CD4)
  l.info(f'Loaded list of epis from file {epitope_path}. Found {len(raw_epis)} epis before cross checking with full sequence.')
  l.info(f'Attempting to load mutations from file {mutationspath} assuming column names protein, Position, Residue.')
  mutations = mutationlist_from_csv(mutationspath, protein_col_name="protein",position_col_name="Position",mutation_col_name="Residue")
  l.info(f'Loading successful. Found {len(mutations)} mutations.')
  l.info(f'Attempting to load protein sequences from file {original_sequences_path}')
  sequence_dict =  read_sequences_from_fasta(original_sequences_path)
  l.info(f'Found {len(sequence_dict)} sequences in that file.')
  # First translate the orfab into the short nsp proteins, then cross check for existance of sequence in respective original sequence.
  l.info(f'For now {len(set([epi.protein for epi in raw_epis]))} distinct proteins are addressed by the epis. Attempting to translate some ORFS to proteins.')
  for epi in raw_epis:
    epi.translate_ORF_to_protein()
  l.info(f'Now {len(set([epi.protein for epi in raw_epis]))} distinct proteins are addressed by the epis.')
  # Now filtering out those epis, where the sequence is not contained in corresponding original wuhan sequence.
  # Check quickly if we attempt to work on epitopes of proteins, where we don't even have the original sequence.
  proteins = list(set([epi.protein for epi in raw_epis]))
  for protein in proteins:
    if protein not in sequence_dict.keys():
      l.critical(f"Error: At least one of the epitopes from protein {protein}, for which no sequence was found in sequence file. Aborting run.")
      quit()
  # Otherwise let's make sure each epitope is actually present in it's respective ancestral protein sequence.
  l.info(f'Checking for each of the {len(raw_epis)} epis if it is present in the ancestral sequence of its suggested protein.')
  epis = [epi for epi in raw_epis if epi.sequence in str(sequence_dict.get(epi.protein).seq)]
  l.info(f'Done. Now epilist has {len(epis)} entries.')
  l.info(f'Now creating the mutated sequences from the ancestral sequences and the mutations.')
  mutated_sequence_dict = generate_mutated_sequences(sequence_dict,mutations)
  l.info(f'And finally first applying the mutations to all epis and then modifying those, where indels happended by aligning against the full mutated protein sequences and padding.')
  for epi in epis: 
    epi.apply_mutations(mutations)
    epi.modify_indel_epitopes(mutated_sequence_dict[epi.protein])
  # From here mostly output and stat generation.
  l.info(f'Turning the epi objects into a dataframe with a reasonable column order.')
  all_epis = [epi.to_dict() for epi in epis ]
  all_epis_df = pd.DataFrame.from_records(all_epis)
  all_epis_df = reorder_dataframe_columns(all_epis_df)
  l.info(f'So far a single epitope with multiple HLA restrictions was treated as multiple epitopes. That way table has {all_epis_df.shape[0]} lines. Now generating also a dataframe with unique epitopes. ')
  # Creating also a dataframe with unique epitopes (duplications bc single epitope with multiple HLA restrictions are here treated as multiple epitopes)
  all_epis_df_unique = generate_unique_epitope_df(all_epis_df)
  l.info(f'Table with unique epitopes created resulting in {all_epis_df_unique.shape[0]} lines.')
  l.info(f'Now generating some statistics.') 
  # TODO: Could easily generate those per each protein if wanted. Just a question on how to arrange in output.
  # Maybe could use a multiindex with outer index being the protein??
  stats = generate_stats_epi_mutations(all_epis_df_unique)
  hist_data = all_epis_df_unique.groupby(by='mutations in this epitope').count().loc[:,'start']
  l.info(f'And finally creating outputfiles. Attempting to write into {output_path} with name stem {name_stem}. \n Writing a file for all epis, all unique epis standard statistics and histogramm data. All ending in .csv and being semicolon separated.')
  outputpath = Path(output_path)
  all_epis_df.to_csv(outputpath.joinpath(name_stem+'_all_non_unique_epitopes.csv'),sep=';',index=False)
  all_epis_df_unique.to_csv(outputpath.joinpath(name_stem+'_all_unique_epitopes.csv'),sep=';',index=False)
  stats.to_csv(outputpath.joinpath(name_stem+'_stats.csv'),sep=';',index=False)
  hist_data.to_csv(outputpath.joinpath(name_stem+'_histogramm_data.csv'),sep=';')

if __name__ == "__main__":
  dfs_CD4 = run(name_stem='CD4',epitope_path='../alba-project/run_data/CD4_epitopes_all_proteins_v2.csv',mutationspath='../alba-project/run_data/mutationlist.csv',original_sequences_path='../alba-project/run_data/ancestral_sequences_mod.txt',output_path='tmp/CD4_run')
  # dfs_CD8 = run(name_stem= 'CD8',epitope_path='run_data/CD8_epitopes_all_proteins_v2.csv',mutationspath='run_data/mutationlist.csv',original_sequences_path='run_data/ancestral_sequences_mod.txt',output_path='run_data/CD8_run')