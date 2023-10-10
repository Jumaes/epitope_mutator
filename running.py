from pathlib import Path
from typing import Dict

import pandas as pd 

from epitope_mutations import Epitope
from input_output import mutationlist_from_csv, read_sequences_from_fasta, epilist_from_csv, generate_mutated_sequences,reorder_dataframe_columns, generate_stats_epi_mutations, generate_unique_epitope_df

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
  # From here mostly output and stat generation.
  # Turning the epi objects into a dataframe with a reasonable column order.
  all_epis = [epi.to_dict() for epi in epis ]
  all_epis_df = pd.DataFrame.from_records(all_epis)
  all_epis_df = reorder_dataframe_columns(all_epis_df)
  # Creating also a dataframe with unique epitopes (duplications bc single epitope with multiple HLA restrictions are here treated as multiple epitopes)
  all_epis_df_unique = generate_unique_epitope_df(all_epis_df)
  # And now generating some stats 
  # TODO: Could easily generate those per each protein if wanted. Just a question on how to arrange in output.
  # Maybe could use a multiindex with outer index being the protein??
  stats = generate_stats_epi_mutations(all_epis_df_unique)
  hist_data = all_epis_df_unique.groupby(by='mutations in this epitope').count().loc[:,'start']

  outputpath = Path(output_path)
  all_epis_df.to_csv(outputpath.joinpath(name_stem+'_all_non_unique_epitopes.csv'),sep=';',index=False)
  all_epis_df_unique.to_csv(outputpath.joinpath(name_stem+'_all_unique_epitopes.csv'),sep=';',index=False)
  stats.to_csv(outputpath.joinpath(name_stem+'_stats.csv'),sep=';',index=False)
  hist_data.to_csv(outputpath.joinpath(name_stem+'_histogramm_data.csv'),sep=';')

if __name__ == "__main__":
  dfs_CD4 = run(name_stem='CD4',epitope_path='run_data/CD4_epitopes_all_proteins_v2.csv',mutationspath='run_data/mutationlist.csv',original_sequences_path='run_data/ancestral_sequences_mod.txt',output_path='run_data/CD4_run')
  dfs_CD8 = run(name_stem= 'CD8',epitope_path='run_data/CD8_epitopes_all_proteins_v2.csv',mutationspath='run_data/mutationlist.csv',original_sequences_path='run_data/ancestral_sequences_mod.txt',output_path='run_data/CD8_run')
