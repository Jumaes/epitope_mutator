import pytest
from copy import deepcopy
from Bio.SeqRecord import SeqRecord
from input_output import epis_from_dicts,epilist_from_csv, mutationlist_from_csv, read_sequences_from_fasta, generate_mutated_sequences
from epitope_mutations import Epitope

test_epitope_1 = {'epitope_id': 1,'protein':'S', 'start':10, 'end':19, 'length':10, 'sequence':'ALCTTSFWFH', 'HLA_restrictions': 'HLA class II'}
test_epitope_2 = {'epitope_id': 1,'protein':'ORF1ab', 'start':1, 'end':10, 'length':10, 'sequence':'ALCTTSFWFH', 'HLA_restrictions': 'HLA class II'}
epitope_list = [test_epitope_1,test_epitope_2]

def test_epis_from_dicts():
    epis = epis_from_dicts(epitope_list)
    assert type(epis[0]) is Epitope
    assert type(epis[1]) is Epitope
    assert epis[0].protein == 'S'
    assert epis[1].protein == 'ORF1ab'


def test_epis_from_dicts_missing_value1():
    with pytest.raises(AssertionError):
        test_epitope_3 = deepcopy(test_epitope_2)
        test_epitope_3.pop('epitope_id')
        epitope_list = [test_epitope_1,test_epitope_3]
        epis = epis_from_dicts(epitope_list)
def test_epis_from_dicts_missing_value2():
    with pytest.raises(AssertionError):
        test_epitope_3 = deepcopy(test_epitope_2)
        test_epitope_3.pop('protein')
        epitope_list = [test_epitope_1,test_epitope_3]
        epis = epis_from_dicts(epitope_list)
def test_epis_from_dicts_missing_value3():
    with pytest.raises(AssertionError):
        test_epitope_3 = deepcopy(test_epitope_2)
        test_epitope_3.pop('sequence')
        epitope_list = [test_epitope_1,test_epitope_3]
        epis = epis_from_dicts(epitope_list)
def test_epis_from_dicts_missing_value4():
    with pytest.raises(AssertionError):
        test_epitope_3 = deepcopy(test_epitope_2)
        test_epitope_3.pop('start')
        epitope_list = [test_epitope_1,test_epitope_3]
        epis = epis_from_dicts(epitope_list)
def test_epis_from_dicts_missing_value5():
    with pytest.raises(AssertionError):
        test_epitope_3 = deepcopy(test_epitope_2)
        test_epitope_3.pop('end')
        epitope_list = [test_epitope_1,test_epitope_3]
        epis = epis_from_dicts(epitope_list)



def test_epilist_from_csv_import():
    # Import generates epis and right number of them
    path = 'testing_data/input_epis.csv'
    epis = epilist_from_csv(path)
    assert all([(type(epi) == Epitope ) for epi in epis])
    assert len(epis) == 14
    epi_starts = [epi.start for epi in epis]
    epi_lengths = [epi.length for epi in epis]
    # Type of start_col_name column and length column int (though input data here not any other data type...)
    assert all([(type(x) == int) for x in epi_starts])
    assert all([(type(x) == int) for x in epi_lengths])

def test_epilist_from_csv_import_insufficient_info_lines():
    # Lines with insufficient information are droppped, but only those.
    path = 'testing_data/input_epis_insufficient_lines.csv'
    epis = epilist_from_csv(path)
    assert all([(type(epi) == Epitope ) for epi in epis])
    epi_ids = [epi.epitope_id for epi in epis]
    assert all([(x in epi_ids) for x in range(1,16)])
    assert not 16 in epi_ids

def test_epilist_from_csv_import_missing_ends():
    # End column gets calculated when not present (from length??? But that is in epi calculted from start and end?!?!?!)
    path = 'testing_data/input_epis_missing_end.csv'
    epis = epilist_from_csv(path)
    epi_ends = [epi.end for epi in epis]
    assert all([(type(x) == int) for x in epi_ends])

def test_epilist_from_csv_import_start_mutstrig_error():
    # Lines with 'mut' in start column are getting dropped.
    path = 'testing_data/input_epis_start_mutstr_error.csv'
    epis = epilist_from_csv(path)
    epi_ids = [epi.epitope_id for epi in epis]
    assert all([(type(epi) == Epitope ) for epi in epis])
    assert len(epis) == 13
    assert not 4 in epi_ids

def test_epilist_from_csv_import_start_lengths_str_input():
    # Import generates epis and right number of them
    path = 'testing_data/input_epis_start_length_strings.csv'
    with pytest.raises(ValueError):
        epis = epilist_from_csv(path)

def test_epilist_from_csv_import_start_lengths_faux_str_input():
    # Import generates epis and right number of them
    path = 'testing_data/input_epis_start_length_faux_strings.csv'
    epis = epilist_from_csv(path)
    assert all([(type(epi) == Epitope ) for epi in epis])
    assert len(epis) == 14
    epi_starts = [epi.start for epi in epis]
    epi_lengths = [epi.length for epi in epis]
    # Type of start_col_name column and length column int (though input data here not any other data type...)
    assert all([(type(x) == int) for x in epi_starts])
    assert all([(type(x) == int) for x in epi_lengths])


def test_mutationlist_from_csv():
    l = mutationlist_from_csv('testing_data/mutationlist.csv', mutation_col_name= 'Residue')
    assert len(l) == 8
    assert l[4]['protein'] == 'M'
    assert l[4]['position'] == 63
    assert l[4]['new'] == 'T'
    assert l[7]['protein'] == 'N'
    assert l[7]['position'] == 31
    assert l[7]['new'] == '-'

# TODO: This test does't work as intended.. Somehow gets messed up with the Error Handling
def test_mutationlist_from_csv_wrong_colnames():
    with pytest.raises(KeyError):
        l = mutationlist_from_csv('testing_data/mutationlist.csv')


def test_read_sequences_from_fasta():
    seqdict = read_sequences_from_fasta('testing_data/sequences.fasta')
    assert type(seqdict)  == dict
    assert type(seqdict['S']) == SeqRecord
    assert all([(n in seqdict.keys()) for n in ['S','H','ORF1ab']])
    assert str(seqdict['S'].seq) ==  'AADDFGGGHSTSR'


def test_generate_mutated_sequences():
    seqdict = read_sequences_from_fasta('testing_data/sequences.fasta')
    mutatonlist  = [{'protein':'S', 'position': 3, 'new': 'CRC'},
                          {'protein':'S', 'position': 10, 'new': 'W'},
                          {'protein':'H', 'position': 3, 'new': 'W'},
                          {'protein':'H', 'position': 5, 'new': '-'}]
    output_dict = generate_mutated_sequences(seqdict,mutatonlist)
    assert len(output_dict.keys()) == 3
    assert all([(type(item)==str) for key,item in output_dict.items()])
    assert output_dict['ORF1ab'] == str(seqdict['ORF1ab'].seq)
    assert output_dict['S'] == 'AACRCDDFGGGHWTSR'
    assert output_dict['H'] == 'AAWDGGGHCRSR'
    # # All sequences actually turned into dict entries
    # assert sequences mutated or not according to mutations and according to protein




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
  """