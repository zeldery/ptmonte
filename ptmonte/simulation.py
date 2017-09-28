'''
Monte Carlo code
written by: Thien-Phuc Tu-Nguyen
Last modified: 2017
File: simulation
'''

import math
import numpy as np
import random
from atom import Atom
from structure import Lattice,Adsorbent,Box
from step import StepTranslation, StepAdd, StepRemove
from forcefield import ForceField

class Simulation:
    def __init__(self):
        self.lattice = None
        self.adsorbent = None
        self.box = None
        self.ff = None
        self.record_en = []
        self.record_adsorb = []
        self.steps = []
        self.p_step = []
        
    def init(self):
        pass
        
    def run(self,n_step):
        for i in range(n_step):
            r = random.choices(range(len(self.steps)),self.p_step)[0]
            self.steps[r].run(lattice = self.lattice,adsorbent =  self.adsorbent,box = self.box)
            
    def to_csv(self):
        pass


class GrandCanonicalSimulation(Simulation):
    def __init__(self):
        Simulation.__init__(self)
        self.temperature = 273.15
        self.pressure = 1.0
        self.d_max = 1.0
        self.mass = 16.0
        
    def init(self,lattice_file,ff_file):
        self.lattice = Lattice()
        self.lattice.read_cif(lattice_file)
        self.lattice.init()
        self.ff = ForceField()
        self.ff.read_raspa_def(ff_file)
        self.ff.init()
        self.adsorbent = Adsorbent()
        self.adsorbent.copy_lattice(self.lattice)
        a = StepTranslation()
        a.ff = self.ff
        a.init(self.d_max,self.temperature)
        self.steps.append(a)
        a = StepAdd()
        a.ff = self.ff
        a.init(self.mass,self.pressure,self.temperature)
        self.steps.append(a)
        a = StepRemove()
        a.ff = self.ff
        a.init(self.mass,self.pressure,self.temperature)
        self.steps.append(a)
        

    

