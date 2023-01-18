
from .DynamicProgramming import DynamicProgramming


class NWAlignment(DynamicProgramming):
    """ implementation of Needleman-Wunsch"""

    space = -2
    MATCH = 1
    MISS_MATCH = -1

    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        super(NWAlignment, self).__init__(len(seq2) + 1, len(seq1) + 1)

    def get_initial_score(self, cell):
        if cell.row == 0 and cell.col != 0:
            return self.matrix[cell.row][cell.col - 1]
        elif cell.row != 0 and cell.col == 0:
            return self.matrix[cell.row - 1][cell.col]
        else:
            return None

    def get_initial_pointer(self, cell):
        if cell.row == 0 and cell.col != 0:
            return cell.col * NWAlignment.space
        elif cell.row != 0 and cell.col == 0:
            return cell.row * NWAlignment.space
        else:
            return 0

    def fill_cell(self, current_cell, top_cell, left_cell, top_left_cell):
        row_space_score = top_cell.score + NWAlignment.space
        col_space_score = left_cell.score + NWAlignment.space
        match_missmatch_score = top_left_cell.score
        # give match/mismatch score
        match_missmatch_score += BLOSUM62.values(self.seq2[current_cell.row - 1], self.seq1[current_cell.col - 1])

        if row_space_score >= col_space_score:
            if match_missmatch_score >= row_space_score:
                current_cell.score = match_missmatch_score
                current_cell.prev_cell = top_left_cell
            else:
                current_cell.score = row_space_score
                current_cell.prev_cell = top_cell
        else:
            if match_missmatch_score >= col_space_score:
                current_cell.score = match_missmatch_score
                current_cell.prev_cell = top_left_cell
            else:
                current_cell.score = col_space_score
                current_cell.prev_cell = left_cell

    def get_traceback(self):
        align1 = []
        align2 = []
        cell = self.get_traceback_starting_cell()
        while not self.traceback_done(cell):
            prev_cell = cell.prev_cell
            if cell.row - prev_cell.row == 1:
                align2.insert(0, self.seq2[cell.row - 1])
            else:
                align2.insert(0, '-')

            if cell.col - prev_cell.col == 1:
                align1.insert(0, self.seq1[cell.col - 1])
            else:
                align1.insert(0, '-')
            cell = prev_cell
        return align1, align2

    def get_traceback_starting_cell(self):
        matrix = self.matrix
        return matrix[len(matrix) - 1][len(matrix[0]) - 1]

    def traceback_done(self, cell):
        return cell.prev_cell is None


class BLOSUM62:
    __instance = None

    ResiduesCodes = {
        'ALA': 'A',
        'ARG': 'R',
        'ASN': 'N',
        'ASP': 'D',
        'CYS': 'C',
        'GLU': 'E',
        'GLN': 'Q',
        'GLY': 'G',
        'HIS': 'H',
        'ILE': 'I',
        'LEU': 'L',
        'LYS': 'K',
        'MET': 'M',
        'PHE': 'F',
        'PRO': 'P',
        'SER': 'S',
        'THR': 'T',
        'TRP': 'W',
        'TYR': 'Y',
        'VAL': 'V'
    }

    def __init__(self):
        with open('BLOSUM62.csv') as f:
            headers = f.readline().strip().split(',')[1:]
            self.subValues = dict()
            for l in f.readlines():
                values = l.split(',')
                resDest = values[0]
                for resSrc, v in zip(headers, values[1:]):
                    self.subValues[(resSrc, resDest)] = int(v)

    @staticmethod
    def values(resA, resB):
        return BLOSUM62.instance().subValues[(resA, resB)]

    @staticmethod
    def instance():
        if not BLOSUM62.__instance:
            BLOSUM62.__instance = BLOSUM62()
        return BLOSUM62.__instance
