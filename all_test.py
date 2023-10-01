#! /home/julian/anaconda3/bin/python
from epitope_mutations import Epitope, DEFAULT_TRANSLATION_DICT
from input_output import epis_from_dicts, epilist_from_csv, mutationlist_from_csv, read_sequence_from_fasta

test_epitope_1 = {'epitope_id': 1,'protein':'S', 'start':10, 'end':19, 'length':10, 'sequence':'ALCTTSFWFH', 'HLA_restrictions': 'HLA class II'}
test_epitope_2 = {'epitope_id': 1,'protein':'ORF1ab', 'start':1, 'end':19, 'length':10, 'sequence':'ALCTTSFWFH', 'HLA_restrictions': 'HLA class II'}

mutation_start = [{'protein':'S', 'position': 10, 'new': 'T'}]
mutation_middle = [{'protein':'S', 'position': 15, 'new': 'T'}]
mutation_end = [{'protein':'S', 'position': 19, 'new': 'T'}]
mutation_end_wrong_protein = [{'protein':'A', 'position': 19, 'new': 'T'}
]
mutation_del = [{'protein':'S', 'position': 12, 'new': '-'}]
mutation_del_middle = [{'protein':'S', 'position': 12, 'new': '-'},
                       {'protein':'S', 'position': 15, 'new': 'T'}]

mutation_insert = [{'protein':'S', 'position': 14, 'new': 'CRC'}]
mutation_insert_middle = [{'protein':'S', 'position': 14, 'new': 'CRC'},
                          {'protein':'S', 'position': 15, 'new': 'W'}]

mutation_del_pos1 = [{'protein':'S', 'position': 10, 'new': '-'}]
mutation_del_pos2 = [{'protein':'S', 'position': 11, 'new': '-'}]

mutation_del_pos19 = [{'protein':'S', 'position': 19, 'new': '-'}]
mutation_del_pos18 = [{'protein':'S', 'position': 18, 'new': '-'}]

mutations_middle_front = [{'protein':'S', 'position': 12, 'new': '-'},
                          {'protein':'S', 'position': 13, 'new': '-'}
                          ]
mutations_middle_back = [{'protein':'S', 'position': 16, 'new': '-'},
                          {'protein':'S', 'position': 17, 'new': '-'}
                          ]

mutation_insert_start = [{'protein':'S', 'position': 10, 'new': 'CRC'}]
mutation_insert_end = [{'protein':'S', 'position': 16, 'new': 'CRC'}]
mutation_insert_dels_after = [{'protein':'S', 'position': 14, 'new': 'CRC'},
                          {'protein':'S', 'position': 15, 'new': '-'},
                          {'protein':'S', 'position': 16, 'new': '-'},
                          ]

mutation_insert_dels_before = [{'protein':'S', 'position': 14, 'new': 'CRC'},
                          {'protein':'S', 'position': 11, 'new': '-'},
                          {'protein':'S', 'position': 12, 'new': '-'},
                          ]


#                    0        1         2         3         4 
#                    1234567890123456789012345678901234567890
original_sequence = 'GGGGGGGGGALCTTSFWFHEEEEEEEEEEEEEEEEEEEEE'
seq_mut_start =     'GGGGGGGGGTLCTTSFWFHEEEEEEEEEEEEEEEEEEEEE'
seq_del_start =     'GGGGGGGGGLCTTSFWFHEEEEEEEEEEEEEEEEEEEEE'
seq_dels_start =    'GGGGGGGGGCTTSFWFHEEEEEEEEEEEEEEEEEEEEE'
seq_del_end =       'GGGGGGGGGALCTTSFWFEEEEEEEEEEEEEEEEEEEEE'
seq_dels_end =      'GGGGGGGGGALCTTSFWEEEEEEEEEEEEEEEEEEEEE'
seq_dels_middle_front =      'GGGGGGGGGALTSFWFHEEEEEEEEEEEEEEEEEEEEE'
seq_dels_middle_end = 'GGGGGGGGGALCTTSFHEEEEEEEEEEEEEEEEEEEEE'
#                    0        1         2         3         4 
#                    1234567890123456789012345678901234567890
seq_ins_middle =    'GGGGGGGGGALCTCRCTSFWFHEEEEEEEEEEEEEEEEEEEEE'
seq_ins_start =     'GGGGGGGGGCRCALCTTSFWFHEEEEEEEEEEEEEEEEEEEEE'
seq_ins_end =       'GGGGGGGGGALCTTSCRCFWFHEEEEEEEEEEEEEEEEEEEEE'
seq_ins_del_after = 'GGGGGGGGGALCTCRCTWFHEEEEEEEEEEEEEEEEEEEEE'
seq_ins_del_before ='GGGGGGGGGATCRCTSFWFHEEEEEEEEEEEEEEEEEEEEE'

#                    0        1         2         3         4 
#                    1234567890123456789012345678901234567890
mutated_seq_set2  = 'GGGGGGGTG-LCLTSCWTTGHRWLIGHKCCCCC-RRRRRRRR'
test_epi1_mod =   '-LCLTSCWTTGH'
test_epi1_final = 'GLCLTSCWTT'

test_epi2_mod = 'HRWLIGHK'
test_epi2_final = 'HRWLIGHK'

test_epi3_mod =   'HRWLIGHKCCCCC-'
test_epi3_final = 'HRWLIGHKCCCCCR'

test_deletion_before_aligned = '-LCL-S-WTTGH'

class TestEpitope:
  def test_mutations_start(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_start)
    assert epi.mod_sequence == 'TLCTTSFWFH'
    assert epi.mutation_counter == 1

  def test_mutations_middle(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_middle)
    assert epi.mod_sequence == 'ALCTTTFWFH'
    assert epi.mutation_counter == 1
 
  def test_mutations_end(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_end)
    assert epi.mod_sequence == 'ALCTTSFWFT'
    assert epi.mutation_counter == 1

  def test_mutations_wrong_protein(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_end_wrong_protein)
    assert epi.mod_sequence == 'ALCTTSFWFH'
    assert epi.mutation_counter == 0
    assert epi.has_del == False
    assert epi.has_ins == False

  def test_mutations_wrong_protein_ignore(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_end_wrong_protein,check_protein=False)
    assert epi.mod_sequence == 'ALCTTSFWFT'
    assert epi.mutation_counter == 1
    assert epi.has_del == False
    assert epi.has_ins == False

  def test_mutations_del(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_del)
    assert epi.mod_sequence == 'AL-TTSFWFH'
    assert epi.mutation_counter == 1
    assert epi.has_ins == False
    assert epi.has_del == True
  
  def test_mutations_del_and_mut(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_del_middle)
    assert epi.mod_sequence == 'AL-TTTFWFH'
    assert epi.mutation_counter == 2
    assert epi.has_ins == False
    assert epi.has_del == True

  def test_mutations_insert(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_insert)
    assert epi.mod_sequence == 'ALCTCRCTSFWFH'
    assert epi.mutation_counter == 1
    assert epi.has_ins == True
    assert epi.has_del == False


  def test_mutations_insert_and_mut(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_insert_middle)
    assert epi.mod_sequence == 'ALCTCRCTWFWFH'
    assert epi.mutation_counter == 2
    assert epi.has_ins == True
    assert epi.has_del == False

  def test_final_seq_no_indel(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_middle)
    epi.modify_indel_epitopes('AAAAAAAAA')
    assert epi.final_mutated_seq == 'ALCTTTFWFH'

  def test_final_seq_gaps_start(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_del_pos1)
    epi.modify_indel_epitopes(seq_del_start)
    assert epi.final_mutated_seq == 'GLCTTSFWFH'
    epi = Epitope(test_epitope_1)
    deletions = mutation_del_pos1 + mutation_del_pos2
    epi.apply_mutations(deletions)
    epi.modify_indel_epitopes(seq_dels_start)
    assert epi.final_mutated_seq == 'GGCTTSFWFH'

  def test_final_seq_gaps_middle1(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_del_pos1)
    epi.modify_indel_epitopes(seq_del_start)
    assert epi.final_mutated_seq == 'GLCTTSFWFH'
    epi = Epitope(test_epitope_1)
    deletions = mutation_del_pos1 + mutation_del_pos2
    epi.apply_mutations(deletions)
    epi.modify_indel_epitopes(seq_dels_start)
    assert epi.final_mutated_seq == 'GGCTTSFWFH'

  def test_final_seq_gaps_middle_front(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutations_middle_front)
    epi.modify_indel_epitopes(seq_dels_middle_front)
    assert epi.final_mutated_seq == 'GGALTSFWFH'

  def test_final_seq_gaps_middle_back(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutations_middle_back)
    epi.modify_indel_epitopes(seq_dels_middle_end)
    assert epi.final_mutated_seq == 'ALCTTSFHEE'

  def test_final_seq_insert_start(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_insert_start)
    epi.modify_indel_epitopes(seq_ins_start)
    assert epi.final_mutated_seq == 'CRCALCTTSF'

  def test_final_seq_insert_middle(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_insert)
    epi.modify_indel_epitopes(seq_ins_middle)
    assert epi.final_mutated_seq == 'TCRCTSFWFH'

  def test_final_seq_insert_back(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_insert_end)
    epi.modify_indel_epitopes(seq_ins_end)
    assert epi.final_mutated_seq == 'ALCTTSCRCF'

  def test_final_seq_insert_dels_after(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_insert_dels_after)
    epi.modify_indel_epitopes(seq_ins_del_after)
    assert epi.final_mutated_seq == 'ALCTCRCTWF'

  def test_final_seq_insert_dels_before(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_insert_dels_before)
    epi.modify_indel_epitopes(seq_ins_del_before)
    assert epi.final_mutated_seq == 'TCRCTSFWFH'

  def test_to_dict_after_mod(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_insert)
    epi.modify_indel_epitopes(seq_ins_middle)
    d = epi.to_dict()
    assert all([
      key in d.keys() for key in ['epitope_id','protein','start','end','length','original_sequence','mutations in this epitope','contains insertion','contains deletion',
            'HLA restrictions','modified not final sequence','mutated_sequence']
    ])

  def test_to_dict_before_mod(self):
    epi = Epitope(test_epitope_1)
    d = epi.to_dict()
    assert all([
      key in d.keys() for key in ['epitope_id','protein','start','end','length','original_sequence','mutations in this epitope','contains insertion','contains deletion',
            'HLA restrictions']
    ])
    assert not ('modified not final sequence' in d.keys())
    assert not ('mutated_sequence' in d.keys())

  def test_to_dict_simple_mod(self):
    epi = Epitope(test_epitope_1)
    epi.apply_mutations(mutation_start)
    epi.modify_indel_epitopes(seq_mut_start)
    d = epi.to_dict()
    assert all([
      key in d.keys() for key in ['epitope_id','protein','start','end','length','original_sequence','mutations in this epitope','contains insertion','contains deletion',
            'HLA restrictions','modified not final sequence','mutated_sequence']
    ])
    assert d['modified not final sequence'] == d['mutated_sequence']

  def translate_ORF_to_protein_start_test():
    epi = Epitope(test_epitope_2)
    epi.translate_ORF_to_protein()
    assert epi.start == test_epitope_2['start']
    assert epi.protein == 'nsp1'

  def translate_ORF_to_protein_middle_test():
    epi = Epitope(test_epitope_2)
    epi.start = 3610
    epi.end = 3619
    epi.translate_ORF_to_protein()
    assert epi.protein == 'nsp6'
    assert epi.start == 41

  def translate_ORF_to_protein_last_test():
    epi = Epitope(test_epitope_2)
    epi.start = 6900
    epi.end = 6919
    epi.translate_ORF_to_protein()
    assert epi.protein == 'nsp16'
    assert epi.start == 102
"""


def input_output_tests() -> Dict[str,pd.DataFrame]:
  print("Testing generation of epis from dictionaries.")
  epilist = [test_epitope_1,test_epitope_2, test_epitope_3]
  epis = epis_from_dicts(epilist)
  for epi in epis:
    print(epi)

  print("Testing generation of epis from input csv file.")
  f = Path('testing_data/CD8_epitopes_restrictions.csv')
  raw_CD8_epis = epilist_from_csv(f)
  # Now using as raw_epis here and then later filter out those, where the sequence is actually not contained in the ancestral wuhan sequence.
  print("Testing generation of mutations from input csv file.")
  fmut = Path('testing_data/mutations_ba286.csv')
  mutations = mutationlist_from_csv(fmut)
  print("Testing generation of first input sequence from fasta file and then generation of mutated sequence with mutations.")
  seqfilepath = Path('testing_data/spike_wt_ba286.txt')
  spike_wuhan =  read_sequence_from_fasta(seqfilepath, 'Spike-Wuhan')
  # Now filtering out those epis, where the sequence is not contained in the spike wuhan sequence.
  CD8_epis = [epi for epi in raw_CD8_epis if epi.sequence in spike_wuhan['Spike-Wuhan']]
  spike_pseudo_epi = Epitope({
    'epitope_id':'Spike',
    'protein': 'S',
    'sequence': spike_wuhan['Spike-Wuhan'],
    'start' : 1,
    'end' : len(spike_wuhan['Spike-Wuhan'])+1,
    'HLA_restrictions':'',
    })
  spike_pseudo_epi.apply_mutations(mutations)
  ba_286_seq = spike_pseudo_epi.mod_sequence.replace('-','')
  '''
  print("Testing if generated mutated sequence from WT matches Albas sequence.")
  albas_ba286_seq = read_sequence_from_fasta(seqfilepath,'spike_BA2.86')['spike_BA2.86']
  print('#'*30 + ' MY BA 286')
  print(ba_286_seq)
  print('#'*30 + ' ALBAS BA 286')
  print(albas_ba286_seq)
  print('#'*30 + ' ALIGNMENT STARTING HERE')
  aligner = PairwiseAligner()
  alignment = aligner.align(ba_286_seq,albas_ba286_seq)[0]
  #print(alignment)
  # return alignment
  # assert ba_286_seq == albas_ba286_seq, 'Generated BA2.86 sequence not the same as the one from Alba'''
  print("Now testing application of mutations from csv file to epis from csv file including alignment and autocomplete in second step.")
  for epi in CD8_epis: 
    epi.apply_mutations(mutations)
    epi.modify_indel_epitopes(mutated_sequence=ba_286_seq)
    # print(epi.to_dict())
  print("Now splitting into those with mutations and those without.")
  no_mutations = [epi.to_dict() for epi in CD8_epis if epi.mutation_counter == 0]
  with_mutations = [epi.to_dict() for epi in CD8_epis if epi.mutation_counter != 0]
  with_deletions = [epi.to_dict() for epi in CD8_epis if epi.has_del]
  with_insertions = [epi.to_dict() for epi in CD8_epis if epi.has_ins]
  all_epis = [epi.to_dict() for epi in CD8_epis ]
  print(f"Of {str(len(CD8_epis))} epis, {str(len(no_mutations))} have not mutations.")
  print(f"{str(len(with_mutations))} have any kind of mutation.")
  print(f"{str(len(with_deletions))} have at least one deletion.")
  print(f"{str(len(with_insertions))} have at least one insetion mutation.")
  print("Returning a dict of dataframes with first those without mutations, then any mutation, then deletions, then insertions.")
  no_mutations_df = pd.DataFrame.from_records(no_mutations)
  with_mutations_df = pd.DataFrame.from_records(with_mutations)
  with_deletions_df = pd.DataFrame.from_records(with_deletions)
  with_insertions_df = pd.DataFrame.from_records(with_insertions)
  all_epis_df = pd.DataFrame.from_records(all_epis)
  df_dict = {'no_mutations':no_mutations_df, 'with_mutations':with_mutations_df,'with_deletions': with_deletions_df, 'with_insertions':with_insertions_df, 'all_epis':all_epis_df}  
  folder =  Path('testing_data/testing_output_data')
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
  for name,dataframe in df_dict.items():
    # Some columns might not exist like 'protein'. So I can't just apply the list as new order, but have to keep just those, which are actually in the dataset.
    new_columns = []
    for column in column_order:
      if column in dataframe.columns:
        new_columns.append(column)
    dataframe = dataframe[new_columns]
    dataframe.to_csv(folder.joinpath(f'{name}.csv'),sep=';',index=False)
  stat_dict = {
      'Type':['All epitopes','With any mutations','With deletions','With insertions'],
      'Number': [all_epis_df.shape[0],with_mutations_df.shape[0],with_deletions_df.shape[0],with_insertions_df.shape[0]],
      'Of total': [all_epis_df.shape[0]/all_epis_df.shape[0], with_mutations_df.shape[0]/all_epis_df.shape[0],with_deletions_df.shape[0]/all_epis_df.shape[0],with_insertions_df.shape[0]/all_epis_df.shape[0]]

    }
  stats = pd.DataFrame(stat_dict)
  stats.to_csv(folder.joinpath('stats.csv'),sep=';',index=False)
  hist_data = all_epis_df.groupby(by='mutations in this epitope').count().loc[:,'start']
  hist_data.to_csv(folder.joinpath('histogramm_data.csv'),sep=';')
  df_dict['stats'] = stats
  df_dict['hist_data'] = hist_data
  return df_dict

if __name__ == "__main__":
  dfs = input_output_tests()


  # TODO:
  # DONE Filter on epitope existing in ancestral sequence
  # DONE Export stats with number of epitopes with mutation, deletions, insertions, % mutations and hist of number of mutations
  # Make sure restriction info is in there; I'll get one file with unique epitop ID HAL combo; export like that and without HLA with unique epitope ID lines and number of different HLA restrctions
  # DONE Export one csv with all
  # DONE Bring start and end into final table"""