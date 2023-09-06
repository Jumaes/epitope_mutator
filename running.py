from pathlib import Path
from typing import Dict

import pandas as pd 

from epitope_mutations import Epitope
from input_output import mutationlist_from_csv, read_sequences_from_fasta, epilist_from_csv, generate_mutated_sequences

COLUMN_NAME_DICT_CD4 = {
        "epitope_id_col_name" : 'Epitope ID', 
        "sequence_col_name" : 'Name',
        "protein_col_name": 'Molecule Parent',
        "start_col_name" : 'starting',
        "end_col_name" : 'Mapped End Position',
        "HLA_restrictions_col_name" : 'MHC restriction',
        "length_col_name" : 'length'}

def run(epitope_path:Path,mutationspath:Path,original_sequences_path:Path, output_path:str, name_stem:str) -> Dict[str,pd.DataFrame]:
  # Now using as raw_epis here and then later filter out those, where the sequence is actually not contained in the ancestral wuhan sequence.
  raw_epis = epilist_from_csv(epitope_path, COLUMN_NAME_DICT_CD4)
  mutations = mutationlist_from_csv(mutationspath, protein_col_name="protein",position_col_name="Position",mutation_col_name="Residue")
  sequence_dict =  read_sequences_from_fasta(original_sequences_path)
  # First translate the orfab into the short nsp proteins, then cross check for existance of sequence in respective original sequence.
  for epi in raw_epis:
    epi.translate_ORF_to_protein()
  # Now filtering out those epis, where the sequence is not contained in corresponding original wuhan sequence.
  # Check quickly if we attempt to work on epitopes of proteins, where we don't even have the original sequence.
  proteins = list(set([epi.protein for epi in raw_epis]))
  for protein in proteins:
    if protein not in sequence_dict.keys():
      print(f"Error: At least one of the epitopes from protein {protein}, for which no sequence was found in sequence file.")
      quit()
  # Otherwise let's make sure each epitope is actually present in it's respective ancestral protein sequence.
  epis = [epi for epi in raw_epis if epi.sequence in str(sequence_dict.get(epi.protein).seq)]
  #   Now creating the mutated sequences from the ancestral sequences and the mutations.
  mutated_sequence_dict = generate_mutated_sequences(sequence_dict,mutations)
  # And finally first applying the mutations to all epis and then modifying those, where indels happended by aligning against the full mutated protein sequences and padding.
  for epi in epis: 
    epi.apply_mutations(mutations)
    epi.modify_indel_epitopes(mutated_sequence_dict[epi.protein])
  # The rest is just outout stuff, conversion to dataframes and generating stats.
  # TODO: A lot of this should actually become functions in input_output
  outputpath = Path(output_path)
  all_epis = [epi.to_dict() for epi in epis ]
  ''' no_mutations = [epi.to_dict() for epi in epis if epi.mutation_counter == 0]
  with_mutations = [epi.to_dict() for epi in epis if epi.mutation_counter != 0]
  with_deletions = [epi.to_dict() for epi in epis if epi.has_del]
  with_insertions = [epi.to_dict() for epi in epis if epi.has_ins]
  print(f"Of {str(len(all_epis))} epis, {str(len(no_mutations))} have not mutations.")
  print(f"{str(len(with_mutations))} have any kind of mutation.")
  print(f"{str(len(with_deletions))} have at least one deletion.")
  print(f"{str(len(with_insertions))} have at least one insetion mutation.")
  print("Returning a dict of dataframes with first those without mutations, then any mutation, then deletions, then insertions.")
  no_mutations_df = pd.DataFrame.from_records(no_mutations)
  with_mutations_df = pd.DataFrame.from_records(with_mutations)
  with_deletions_df = pd.DataFrame.from_records(with_deletions)
  with_insertions_df = pd.DataFrame.from_records(with_insertions)'''
  all_epis_df = pd.DataFrame.from_records(all_epis)
  # df_dict = {'no_mutations':no_mutations_df, 'with_mutations':with_mutations_df,'with_deletions': with_deletions_df, 'with_insertions':with_insertions_df, 'all_epis':all_epis_df}  
  column_order = [
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
      'HLA restrictions'
      ]
  #for name,dataframe in df_dict.items():
    # Some columns might not exist like 'protein'. So I can't just apply the list as new order, but have to keep just those, which are actually in the dataset.
  new_columns = []
  for column in column_order:
    if column in all_epis_df.columns:
      new_columns.append(column)
  all_epis_df = all_epis_df[new_columns]
  all_epis_df.to_csv(outputpath.joinpath(name_stem+'_all_non_unique_epitopes.csv'),sep=';',index=False)
  all_epis_df_unique = all_epis_df.drop_duplicates(subset='original_sequence')
  all_epis_df_unique.to_csv(outputpath.joinpath(name_stem+'_all_unique_epitopes.csv'),sep=';',index=False)
  no_unique_epis = all_epis_df_unique.shape[0]
  no_mutations = all_epis_df_unique[all_epis_df_unique['mutations in this epitope']>0].shape[0]
  no_deletions = all_epis_df_unique[all_epis_df_unique['contains deletion'] == True].shape[0]
  no_insertions = all_epis_df_unique[all_epis_df_unique['contains insertion'] == True].shape[0]
  stat_dict = {
      'Type':['All epitopes','With any mutations','With deletions','With insertions'],
      'Number': [no_unique_epis,no_mutations, no_deletions, no_insertions],
      'Of total': [no_unique_epis/no_unique_epis, no_mutations/no_unique_epis, no_deletions/no_unique_epis, no_insertions/no_unique_epis]

    }
  stats = pd.DataFrame(stat_dict)
  stats.to_csv(outputpath.joinpath(name_stem+'_stats.csv'),sep=';',index=False)
  hist_data = all_epis_df_unique.groupby(by='mutations in this epitope').count().loc[:,'start']
  hist_data.to_csv(outputpath.joinpath(name_stem+'_histogramm_data.csv'),sep=';')

if __name__ == "__main__":
  dfs_CD4 = run(name_stem='CD4',epitope_path='run_data/CD4_epitopes_all_proteins_v2.csv',mutationspath='run_data/mutationlist.csv',original_sequences_path='run_data/ancestral_sequences_mod.txt',output_path='run_data/CD4_run')
  dfs_CD8 = run(name_stem= 'CD8',epitope_path='run_data/CD8_epitopes_all_proteins_v2.csv',mutationspath='run_data/mutationlist.csv',original_sequences_path='run_data/ancestral_sequences_mod.txt',output_path='run_data/CD8_run')
