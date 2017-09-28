'''
Monte Carlo code
written by: Thien-Phuc Tu-Nguyen
Last modified: 2017
File: atom
'''


class Atom:
    def __init__(self):
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.a_type = ''
        self.p_index = 0
        self.charge = 0.0
        
    def copy(self):
        new_atom = Atom()
        new_atom.x = self.x
        new_atom.y = self.y
        new_atom.z = self.z
        new_atom.a_type = self.a_type
        new_atom.p_index = self.p_index
        new_atom.charge = self.charge
        return new_atom

