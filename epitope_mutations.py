import logging
from dataclasses import dataclass
from Bio.Align import PairwiseAligner as pa
from typing import Tuple, Dict, List

l = logging.getLogger('epitope_mutator.epitope')
# sh = logging.StreamHandler()
# l.addHandler(sh)
# l.setLevel(logging.DEBUG)

DEFAULT_TRANSLATION_DICT= {
    'ORF1ab': [
        [1, 181,819,2764,3264,3570,3860,3943,4141,4254,4393,5325,5926,6453,6799,7096],
        ["nsp1","nsp2","nsp3","nsp4","nsp5","nsp6","nsp7","nsp8","nsp9","nsp10","nsp12","nsp13","nsp14","nsp15","nsp16"]
    ]
}

@dataclass
class Epitope:
    epitopeinfo: dict
    mod_sequence : str = None
    mutation_counter: int = 0
    final_mutated_seq: str = None
    has_del: bool = False
    has_ins: bool = False
    epitope_id: str = ''
    protein: str = ''
    sequence: str = ''
    start: int = 0
    HLA_restrictions : str = ''
    end: int = 0
    
    def __post_init__(self):
        self.epitope_id = self.epitopeinfo['epitope_id']
        self.protein = self.epitopeinfo['protein']
        self.sequence = self.epitopeinfo['sequence']
        self.start = int(self.epitopeinfo['start'])
        self.HLA_restrictions = self.epitopeinfo['HLA_restrictions']
        if self.epitopeinfo.get('end'):
            self.end = int(self.epitopeinfo['end'])
        else:
            self.end = self.start + self.length - 1  
        l.debug(f'Epitope initialized: {self}')
    
    @property
    def length(self):
        l = len(self.sequence)
        return l

    def apply_mutations(self, mutations:list, check_protein:bool=True):
        l.debug(f'Applying mutations to epitope startin at {self.start} in protein {self.protein} with sequence {self.sequence}.')
        self.mod_sequence = self.sequence
        self.mutation_counter = 0
        # NB: This assumes the 'end' is indeed the last position of the epitope aka it is still part of the epitope.
        epirange = range(self.start, self.end +1 )
        mod_seq_list = list(self.mod_sequence)
        for mutation in mutations:
            # Double checking that the protein of the epitope and the mutation are the same. 
            # If not skip to next mutation.
            # This is not mandatory since funtion might be applied to a set of only one protein, where no protein names are given.
            if check_protein:
                # NB using 'lower()' comparison here makes this case insensitive. This is under the assumption that no two protein names are identical except for case and that case misspellings (nsp3 vs NSP3; ORF1a vs orf1a) are frequent.
                if mutation['protein'].lower() != self.protein.lower():
                    l.debug(f'\t\t Mutation {mutation} not in protein of this epitope.')
                    continue
            if mutation['position'] in epirange:
                    l.debug(f'\t Found mutation {mutation} in protein and correct range for epitope. Applying. ')
                    self.mutation_counter += 1
                    try:
                        if len(mutation['new']) == 1:
                            if mutation['new'] == '-':
                                self.has_del = True
                            mod_seq_list[mutation['position']-self.start] = mutation['new']
                        else:
                            mod_seq_list[mutation['position']-self.start] = mutation['new'] + mod_seq_list[mutation['position']-self.start]
                            self.has_ins = True
                    except(IndexError):
                        l.critical(f"Index Error occuring at mutation: {mutation['position']} and epitope: {self.epitope_id} of protein {self.protein} and position {str(self.start)}.")
                        quit()
                    self.mod_sequence = ''.join(mod_seq_list)
                    l.debug(f'\t Original sequence {self.sequence} turned into {self.mod_sequence} by mutation {mutation}.')
                    # In case there were only mutations and no indels.
                    if not '-' in self.mod_sequence and len(self.mod_sequence) == self.length:
                        l.debug(f'\t Mutation {mutation} seems simple mutation without deletions or insertions. No further modification of mutated epitope needed.')
                        self.mutated_sequence = self.mod_sequence 
                    else: 
                        l.debug(f'\t\tMutation {mutation} seems to contain deletions and/or insertions. Further modificationn of mutated epitope might be necessary.')

    def modify_indel_epitopes(self, mutated_sequence:str):
        #This is only needed, if there are insertions or deletions in the epitope (to my knowledge)
        if not (self.has_del or self.has_ins):
            self.final_mutated_seq = self.mod_sequence
            return
        l.debug(f'Need to modify final sequence of epitope {self.sequence} from modified sequence {self.mod_sequence} because of indels.')
        l.debug(f'Aligning mod epitope to original epitope to find largest unmodified core.')
        mod_start, mod_end, original_start, original_end  = self.get_unchanged_epitope_core_borders()
        len_padding_needed_start = original_start
        len_padding_needed_end = self.length - original_end
        l.debug(f'Unchanged core identified as {self.unchanged_core}.')
        l.debug(f'Padding needed in front of unchanged core {len_padding_needed_start} and after unchanged core {len_padding_needed_end}.')
        l.debug(f'Now checking position of unchanged core in mutated sequence.')
        epitope_start, epitope_end, sequence_start, sequence_end = self.align_unchanged_core_mutated_sequence(mutated_sequence)
        l.debug(f'Found unchanged core in mutated sequence starting at {sequence_start} and ending at {sequence_end}.')
        self.final_mutated_seq = mutated_sequence[sequence_start-len_padding_needed_start:sequence_end+len_padding_needed_end]
        l.debug(f'Final sequence is {self.final_mutated_seq}.')

    def get_unchanged_epitope_core_borders(self) -> Tuple[int,int,int,int]:
        modified_epitope = self.mod_sequence.replace('-','')
        aligner = pa()
        # We are looking to get the one best matching piece of the epitope aligned without any gaps.
        # Therefore making sure it's a local alignment and gap score through the roof. 
        aligner.mode = 'local'
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -10
        aligner.mismatch_score = -5
        alignments  = aligner.align(self.sequence,modified_epitope)
        best_align = alignments[0]
        l.debug(best_align)
        # Following should return tuples of start and end positions for all aligned parts.
        # Because of high gap scores, there should be only one aligned part.
        starts, ends = best_align.path
        original_start, mod_start = [int(x) for x in starts]
        original_end, mod_end = [int(x) for x in ends]
        # NB: The "ends" are exclusive like in typical python slicing and starting with 0. Alignment of only amino acid 1 would be "(0,1)".
        self.unchanged_core = self.sequence[original_start:original_end]
        return mod_start, mod_end, original_start, original_end 

    def align_unchanged_core_mutated_sequence(self,mutated_sequence:str) -> Tuple[int, int, int, int]:
        aligner = pa()
        # We are looking to get the one best matching piece of the epitope aligned without any gaps.
        # Therefore making sure it's a local alignment and gap score through the roof. 
        aligner.mode = 'local'
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -10
        aligner.mismatch_score = -5
        alignments  = aligner.align(mutated_sequence,self.unchanged_core)
        best_align = alignments[0]
        l.debug(best_align)
        # Following should return tuples of start and end positions for all aligned parts.
        # Because of high gap scores, there should be only one aligned part.
        starts, ends = best_align.path
        sequence_start, epitope_start = [int(x) for x in starts]
        sequence_end, epitope_end = [int(x) for x in ends]
        # NB: The "ends" are exclusive like in typical python slicing. Alignment of only amino acid 1 would be "(0,1)".
        return epitope_start, epitope_end, sequence_start, sequence_end 
    
    @staticmethod
    def gaps_before(start=int, sequence=str)->int:
        gaps = sequence[:start+1].count('-')
        if not gaps:
            return start
        else:
            newgaps = gaps
            while not newgaps == 0:
                newstart = start + gaps
                newgaps = sequence[:newstart+1].count('-') - gaps
                gaps = gaps + newgaps
            return start + gaps
        

    def to_dict(self) -> dict:
        #assert self.final_mutated_seq, 'This method generates epitops with their final mutated sequence. Current epitope does not (yet) contain that'
        d = {
            'epitope_id':self.epitope_id,
            'protein':self.protein, 
            'start':self.start, 
            'end':self.end, 
            'length':self.length, 
            'original_sequence':self.sequence, 
            'mutations in this epitope' : self.mutation_counter,
            'contains insertion': self.has_ins,
            'contains deletion' : self.has_del,
            'HLA restrictions' : self.HLA_restrictions}
        if self.mod_sequence:
            d['modified not final sequence'] = self.mod_sequence
        if self.final_mutated_seq:
            d['mutated_sequence'] = self.final_mutated_seq

        return d
    
    def translate_ORF_to_protein(self, translation_dict:Dict[str,List[list]]=DEFAULT_TRANSLATION_DICT):
        l.debug(f'Translating orf {self.protein} into protein.')
        if self.protein in translation_dict.keys():
            ranges = translation_dict.get(self.protein)[0]
            names = translation_dict.get(self.protein)[1]
            for num,r in enumerate(ranges[:-1]):
                if r <= self.start  < ranges[num+1]:
                    new_name = names[num]
                    new_start = self.start - r + 1
                    new_end = new_start + self.length -1
                    l.debug(f"Translated ORF {self.protein} to protein {new_name} with original start {self.start} to new start {new_start}.")
                    self.start = new_start
                    self.protein = new_name
                    self.end = new_end
                    return
