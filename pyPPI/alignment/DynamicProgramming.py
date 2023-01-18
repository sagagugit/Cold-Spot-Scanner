"""Sequence alignment"""
import abc


class Cell(object):
    def __init__(self, r, c):
        self.prev_cell = None
        self.score = 0
        self.row = r
        self.col = c

    def set_prev_cell(self, prev):
        self.prev_cell = prev


class DynamicProgramming(object):
    def __init__(self, rows, cols):
        self.matrix = []
        for r in range(0, rows):
            row = list([Cell(r, c) for c in range(0, cols)])
            self.matrix.append(row)

        self.init_scores()
        self.init_pointers()

    def init_scores(self):
        for row in self.matrix:
            for col in row:
                col.score(self.get_initial_score(col))

    def init_pointers(self):
        for row in self.matrix:
            for col in row:
                col.set_prev_cell(self.get_initial_pointer(col))

    @abc.abstractmethod
    def get_initial_score(self, cell):
        return

    @abc.abstractmethod
    def get_initial_pointer(self, cell):
        return

    def fillIn(self):
        for row in self.matrix[1:]:
            for cell in row[1:]:
                cellAbove = self.matrix[cell.row - 1][cell.col]
                cellToLeft = self.matrix[cell.row][cell.col - 1]
                cellAboveLeft = self.matrix[cell.row - 1][cell.col - 1]
                self.fill_cell(cell, cellAbove, cellToLeft, cellAboveLeft)

    @abc.abstractmethod
    def fill_cell(self, current_cell, top_cell, left_cell, top_left_cell):
        return

    @abc.abstractmethod
    def get_traceback(self):
        return
