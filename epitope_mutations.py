from Bio.Align import PairwiseAligner as pa
from typing import Tuple, Dict, List




DEFAULT_TRANSLATION_DICT= {
    'ORF1ab': [
        [1, 181,819,2764,3264,3570,3860,3943,4141,4254,4393,5325,5926,6453,6799,7096],
        ["nsp1","nsp2","nsp3","nsp4","nsp5","nsp6","nsp7","nsp8","nsp9","nsp10","nsp12","nsp13","nsp14","nsp15","nsp16"]
    ]
}

class Epitope:
    def __init__(self, epitopteinfo:dict):
        self.epitope_id = epitopteinfo['epitope_id']
        self.protein = epitopteinfo['protein']
        self.sequence = epitopteinfo['sequence']
        self.start = int(epitopteinfo['start'])
        self.HLA_restrictions = epitopteinfo['HLA_restrictions']
        self.mod_sequence = None
        self.mutation_counter = 0
        self.final_mutated_seq = None
        self.has_del = False
        self.has_ins = False
        if epitopteinfo.get('end'):
            self.end = int(epitopteinfo['end'])
        else:
            self.end = self.start + self.length - 1  

    def __repr__(self) -> str:
        #print ('Epitope: '+self.__name__)
        rep = 'Epitope ID: ' + str(self.epitope_id) + '\n'
        rep += 'Start: ' + str(self.start) + '\n'
        rep += 'End: '+ str(self.end) + '\n'
        rep +=  'Length: ' +str(self.length) + '\n'
        rep += 'Sequence: ' + self.sequence + '\n'
        if self.protein:
            rep += 'Protein: ' +self.protein + '\n'
        if self.mod_sequence:
            rep += 'Modified sequence: ' + self.mod_sequence + '\n'
        if self.mutation_counter:
            rep += 'Number of mutations: ' + str(self.mutation_counter) + '\n'
        if self.final_mutated_seq:
            rep += 'Final mutated sequence: ' + str(self.final_mutated_seq) + '\n'
        return rep
    
    @property
    def length(self):
        l = len(self.sequence)
        return l

    def apply_mutations(self, mutations:list, check_protein:bool=True):
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
                    continue
            if mutation['position'] in epirange:
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
                        print (f"Index Error occuring at mutation: {mutation['position']} and epitope: {self.epitope_id} of protein {self.protein} and position {str(self.start)}.")
                        quit()
                    self.mod_sequence = ''.join(mod_seq_list)
                    # In case there were only mutations and no indels.
                    if not '-' in self.mod_sequence and len(self.mod_sequence) == self.length:
                        self.mutated_sequence = self.mod_sequence  

    def modify_indel_epitopes(self, mutated_sequence:str):
        #This is only needed, if there are insertions or deletions in the epitope (to my knowledge)
        if not (self.has_del or self.has_ins):
            self.final_mutated_seq = self.mod_sequence
            return
        # Making sure there are no gaps in the input mutated sequence.
        mutated_sequence = mutated_sequence.replace('-','')
        # In case the epitope underwent indels, we need to bring the epitope back to the original length
        # This is done by aligning it to the mutated full sequence such that only perfect matches align without gaps.
        # Then adding amino acids at beginning and end as applicable to reach original length while keeping amino acids in register.
        epitope_start, epitope_end, sequence_start, sequence_end = self.realign_mut_epitope(mutated_sequence)
        aligned_epitope_core = self.mod_sequence.replace('-','')[epitope_start:epitope_end+1]
        deletions_before_aligned = self.gaps_before(epitope_start,self.mod_sequence)
        len_padding_needed_start = epitope_start + deletions_before_aligned
        if len_padding_needed_start >= 0:
            padding_start = mutated_sequence[sequence_start-len_padding_needed_start:sequence_start]
            final_start_core = padding_start + aligned_epitope_core
        else:
            final_start_core = aligned_epitope_core[abs(len_padding_needed_start):]
        len_padding_needed_end = self.length - (epitope_end +deletions_before_aligned)
        if len_padding_needed_end >= 0:
            padding_end = mutated_sequence[sequence_end+1:sequence_end+1+len_padding_needed_end]
            self.final_mutated_seq = final_start_core + padding_end
        else:
            self.final_mutated_seq = final_start_core[:len_padding_needed_end]


    def realign_mut_epitope(self,mutated_sequence:str) -> Tuple[int, int, int, int]:
        # TODO: This needs major overhaul.
        # Should be first create dict from position in original epitope to position in ungapped new epitope.
        # Then align unagapped new to old sequence and make sure we only align the identical part.
        # Get the borders of that align on the epitope, translate back to original epitope positions.
        # Then using those borders decide on the padding left and right.
        # Finally look for the exact sequence of the ungapped mutated (which aligned to the original) in the mutated full sequence.
        # Then apply the padding left and right.
        modified_epitope = self.mod_sequence.replace('-','')
        aligner = pa()
        # We are looking to get the one best matching piece of the epitope aligned without any gaps.
        # Therefore making sure it's a local alignment and gap score through the roof. 
        aligner.mode = 'local'
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -10
        aligner.mismatch_score = -5
        alignments  = aligner.align(mutated_sequence,modified_epitope)
        best_align = alignments[0]
        #print (best_align)
        # Following should return tuples of start and end positions for all aligned parts.
        # Because of high gap scores, there should be only one aligned part.
        starts, ends = best_align.path
        sequence_start, epitope_start = [int(x) for x in starts]
        sequence_end, epitope_end = [int(x) for x in ends]
        '''target_start_ends, query_start_ends = best_align.path
        assert len(query_start_ends) == 1 , print('Alignment of eiptope lead to gaps. Something is wrong here.')
        epitope_start, epitope_end = [int(x) for x in query_start_ends]
        sequence_start, sequence_end = [int(x) for x in query_start_ends]
        '''# NB: The "ends" are exclusive like in typical python slicing. Alignment of only amino acid 1 would be "(0,1)".
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
        if self.protein in translation_dict.keys():
            ranges = translation_dict.get(self.protein)[0]
            names = translation_dict.get(self.protein)[1]
            for num,r in enumerate(ranges[:-1]):
                if r <= self.start  < ranges[num+1]:
                    new_name = names[num]
                    new_start = self.start - r + 1
                    new_end = new_start + self.length -1
                    self.start = new_start
                    self.protein = new_name
                    self.end = new_end
                    return
