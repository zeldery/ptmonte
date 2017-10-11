'''
Monte Carlo code
written by: Thien-Phuc Tu-Nguyen
Last modified: 2017
File: simulation
'''

import math
import numpy as np
import random
from .atom import Atom
from .structure import Lattice,Adsorbent,Box
from .step import StepTranslation, StepAdd, StepRemove
from .forcefield import ForceField

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
        
    def single_run(self): # Change this function if you want to record the potential
        r = random.choices(range(len(self.steps)),self.p_step)[0]
        self.steps[r].run(lattice = self.lattice,adsorbent = self.adsorbent, box = self.box)
        
    def run(self,n_step):
        for i in range(n_step):
            self.single_run()
            
    def calibrate(self,n_step):
        pass
    
    def reset(self):
        for step in self.steps:
            step.reset()
        self.record_en = []
        self.record_adsorb = []
            
    def to_csv(self):
        pass


class GrandCanonicalSimulation(Simulation):
    def __init__(self):
        Simulation.__init__(self)
        self.temperature = 273.15
        self.pressure = 1.0
        self.d_max = 1.0
        self.mass = 16.0
        self.p_step = [0.4,0.3,0.3]
        
        
    def init(self,lattice_file,ff_file,a_type):
        # Prepare the lattice
        self.lattice = Lattice()
        self.lattice.read_cif(lattice_file)
        self.lattice.init()
        # Prepare the force field
        self.ff = ForceField()
        self.ff.read_raspa_def(ff_file)
        self.ff.init()
        # Generate the adsorbent
        self.adsorbent = Adsorbent()
        self.adsorbent.copy_lattice(self.lattice)
        # Set index for the lattice
        self.ff.set_index(self.lattice)
        # Add translation step
        a = StepTranslation()
        a.ff = self.ff
        a.init(self.d_max,self.temperature)
        self.steps.append(a)
        # Add addition step
        atom = Atom()
        atom.a_type = a_type
        self.ff.set_atom(atom)
        a = StepAdd()
        a.ff = self.ff
        a.init(atom,self.mass,self.pressure,self.temperature)
        self.steps.append(a)
        a = StepRemove()
        a.ff = self.ff
        a.init(self.mass,self.pressure,self.temperature)
        self.steps.append(a)
        
    def single_run(self):
        Simulation.single_run(self)
        # self.record_en.append(self.ff.interaction(self.adsorbent,self.lattice)+self.ff.box(self.adsorbent) )
        self.record_adsorb.append(len(self.adsorbent.atoms))
    
class GibbsEnsembleSimulation(Simulation):
    def __init__(self):
        Simulation.__init__(self)
        self.temperature = 273.15
        self.pressure = 1.0
        self.d_max = 1.0
        self.d_logV = 0.1
        self.mass = 1.0
        self.p_step = []
        
    def init(self, lattice_file, ff_file, a_type):
        pass
