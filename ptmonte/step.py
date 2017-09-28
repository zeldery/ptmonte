'''
Monte Carlo code
written by: Thien-Phuc Tu-Nguyen
Last modified: 2017
File: step
'''

import random
import numpy as np
import math
from atom import Atom
from constants import *

# General class for Step
class Step:
    def __init__(self):
        self.acceptance = 0
        self.total = 0
        self.temperature = 275.15
        self.ff = None # Reference to ForceField
        
    def init(self):
        pass
    
    def run(self,adsorbent = None, lattice = None, box = None):
        pass
        
    def accept_rate(self):
        return self.acceptance / self.total

# Translation step
class StepTranslation(Step):
    def __init__(self):
        Step.__init__(self)
        self.d_max = 0.0
        
    def init(self,d_max, temperature = 273.15):
        self.d_max = d_max
        self.temperature = temperature
        
    def run(self,adsorbent = None, lattice = None, box = None):
        if box == None:
            container = adsorbent
        else:
            container = box
        if len(container.atoms) == 0:
            return
        atom = container.atoms.pop(random.randrange(len(container.atoms)))
        old_en = self.ff.one_atom(atom,container)
        if lattice != None:
            old_en += self.ff.one_atom(atom,lattice)
        # Move atom
        new_atom = atom.copy()
        new_atom.x += self.d_max * random.random()
        new_atom.y += self.d_max * random.random()
        new_atom.z += self.d_max * random.random()
        container.check(new_atom)
        new_en = self.ff.one_atom(new_atom,container)
        if lattice != None:
            new_en += self.ff.one_atom(new_atom,lattice)
        self.total += 1
        if random.random() < math.exp((old_en - new_en) / self.temperature): # Energy conversion
            container.atoms.append(new_atom)
            self.acceptance += 1
        else:
            container.atoms.append(atom)

# Add a particle to a box / adsorbent class
class StepAdd(Step):
    def __init__(self):
        Step.__init__(self)
        self.pressure = 0.0
        self.mu = 0.0
        self.lamb3 = 0.0 # De broglie wavelength
    
    def init(self,mass,pressure,temperature = 273.15):
        self.pressure = pressure
        self.temperature = temperature
        self.lamb3 = PLANK / math.sqrt(2*math.pi*mass*BOLTZMANN*temperature)
        self.lamb3 = self.lamb3**3
        self.mu = BOLTZMANN * temperature * math.log(self.lamb3 * pressure / BOLTZMANN / temperature)
        
        
    def run(self,adsorbent = None, lattice = None, box = None):
        x = random.random()
        y = random.random()
        z = random.random()
        coord = np.dot(lattice.to_cartesian,[x,y,z])
        
        atom = Atom()
        atom.x = coord[0]
        atom.y = coord[1]
        atom.z = coord[2]
        en = self.ff.one_atom(atom,adsorbent) + self.ff.one_atom(atom,lattice)
        prop = lattice.volume/self.lamb3/(len(adsorbent.atoms)+1)*math.exp((self.mu-en)/self.temperature) # Energy conversion
        self.total += 1
        if random.random() < prop:
            self.acceptance += 1
            adsorbent.atoms.append(atom)

# Remove a particle from a box/adsorbent
class StepRemove(Step):
    def __init__(self):
        Step.__init__(self)
        self.pressure = 0.0
        self.mu = 0.0
        self.lamb3 = 0.0
    
    def init(self,mass,pressure,temperature = 273.15):
        self.pressure = pressure
        self.temperature = temperature
        self.lamb3 = PLANK / math.sqrt(2*math.pi*mass*BOLTZMANN*temperature)
        self.lamb3 = self.lamb3**3
        self.mu = BOLTZMANN * temperature * math.log(self.lamb3 * pressure / BOLTZMANN / temperature)
        
    def run(self,adsorbent = None, lattice = None, box = None):
        if len(adsorbent.atoms) > 0 :
            atom = adsorbent.atoms.pop(random.randrange(len(adsorbent.atoms)))
            en = self.ff.one_atom(atom,adsorbent) + self.ff.one_atom(atom,lattice)
            prop = self.lamb3*(len(adsorbent.atoms)+1) / lattice.volume * math.exp((en-self.mu)/self.temperature) # Energy conversion
            self.total += 1
            if random.random() < prop:
                self.acceptance += 1
            else:
                adsorbent.atoms.append(atom)
    

