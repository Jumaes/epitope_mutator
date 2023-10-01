from epitope_mutations import Epitope

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



  # TODO:
  # DONE Filter on epitope existing in ancestral sequence
  # DONE Export stats with number of epitopes with mutation, deletions, insertions, % mutations and hist of number of mutations
  # Make sure restriction info is in there; I'll get one file with unique epitop ID HAL combo; export like that and without HLA with unique epitope ID lines and number of different HLA restrctions
  # DONE Export one csv with all
  # DONE Bring start and end into final table"""