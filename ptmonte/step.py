'''
Monte Carlo code
written by: Thien-Phuc Tu-Nguyen
Last modified: 2017
File: step
'''

import random
import numpy as np
import math
from .atom import Atom
from .constants import *

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
        
    def reset(self):
        self.total = 0
        self.acceptance = 0

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
        new_atom.x += self.d_max * (random.random()-0.5)
        new_atom.y += self.d_max * (random.random()-0.5)
        new_atom.z += self.d_max * (random.random()-0.5)
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
        self.atom = Atom()
    
    def init(self, atom, mass, pressure, temperature = 273.15):
        self.pressure = pressure
        self.temperature = temperature
        self.lamb3 = DE_BROGLIE / math.sqrt(mass*temperature)
        self.lamb3 = self.lamb3**3
        self.mu = temperature * math.log(self.lamb3 * pressure / BOLTZMANN_ANGSTROM / temperature) # Unit mu/kB in K
        self.atom = atom.copy()
        
        
    def run(self,adsorbent = None, lattice = None, box = None):
        x = random.random()
        y = random.random()
        z = random.random()
        coord = np.dot(lattice.to_cartesian,[x,y,z])
        
        atom = self.atom.copy()
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
    
    def init(self, mass, pressure, temperature = 273.15):
        self.pressure = pressure
        self.temperature = temperature
        self.lamb3 = DE_BROGLIE / math.sqrt(mass*temperature)
        self.lamb3 = self.lamb3**3
        self.mu = temperature * math.log(self.lamb3 * pressure / BOLTZMANN_ANGSTROM / temperature) # Unit mu/kB in K
        
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
    

# Change volume step
class StepVolume(Step):
    def __init__(self):
        Step.__init__(self)
        self.pressure = 0.0
        
    def init(self,pressure,d_logV,temperature = 273.15):
        self.pressure = pressure
        self.d_logV = d_logV
        self.temperature = temperature
        
    def run(self,adsorbent = None, lattice = None, box = None):
        en_old = self.ff.box(box)
        side_old = box.side
        logV = 3 * math.log(side_old)
        logV += self.d_logV * (random.random() - 0.5)
        side_new = math.exp(logV/3)
        k = side_new / side_old
        for atom in box.atoms:
            atom.x *= k
            atom.y *= k
            atom.z *= k
        en_new = self.ff.box(box)
        n = len(box.atoms)
        prop = k **(3*n + 3) * math.exp( 
            (en_old - en_new + self.pressure*(side_old ** 3 - side_new **3) / BOLTZMANN_ANGSTROM ) / self.temperature )
        self.total += 1
        if random.random() < prop:
            self.acceptance += 1
            box.side = side_new
            box.volume = side_new ** 3
        else:
            for atom in box.atoms:
                atom.x /= k
                atom.y /= k
                atom.z /= k

# Swap the particle between box
class StepSwap(Step):
    def __init__(self):
        Step.__init__(self)
        self.total_to_box = 0
        self.total_to_adsorbent = 0
        self.acceptance_to_box = 0
        self.acceptance_to_adsorbent = 0
    
    def init(self, temperature = 273.15):
        self.temperature = temperature
        
    def run(self,adsorbent = None, lattice = None, box = None):
        self.total += 1
        if random.random() < 0.5:
            # Change from the box to the adsorbent
            r = random.randrange(len(box.atoms))
            coord = np.dot(adsorbent.to_cartesian,[random.random(),random.random(),random.random()])
            atom = box.atoms[r].copy()
            atom.x = coord[0]
            atom.y = coord[1]
            atom.z = coord[2]
            en_old = ff.one_atom(box.atom[r],box)
            en_new = ff.one_atom(atom,adsorbent) + ff.one_atom(atom,lattice)
            prop = adsorbent.volume / (len(adsorbent.atoms)+1) * len(box.atoms) / box.volume * math.exp((en_old - en_new) / self.temperature)
            self.total_to_adsorbent += 1
            if random.random() < prop:
                self.acceptance += 1
                self.acceptance_to_adsorbent += 1
                adsorbent.append(atom)
                del box.atoms[r]
            
        else:
            r = random.randrange(len(adsorbent.atoms))
            atom = adsorbent.atoms[r].copy()
            atom.x = random.random() * box.side
            atom.y = random.random() * box.side
            atom.z = random.random() * box.side
            en_old = ff.one_atom(box.atom[r],adsorbent) + ff.one_atom(box.atom[r],lattice)
            en_new = ff.one_atom(atom,box)
            prop = box.volume / (len(box.atoms)+1) * len(adsorbent.atoms) / adsorbent.volume * math.exp((en_old - en_new) / self.temperature)
            self.total_to_box += 1
            if random.random() < prop:
                self.acceptance += 1
                self.acceptance_to_box += 1
                box.append(atom)
                del adsorbent.atoms[r]
    

    
    
