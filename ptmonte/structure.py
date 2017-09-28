'''
Monte Carlo code
written by: Thien-Phuc Tu-Nguyen
Last modified: 2017
File: structure
'''

import random
import numpy as np
import math
from atom import Atom

# Class for Lattice
class Lattice:
    def __init__(self):
        self.atoms = []
        self.internal_atoms = []
        self.a = 0.0
        self.b = 0.0
        self.c = 0.0
        self.alpha = 0.0
        self.beta = 0.0
        self.gamma = 0.0
        self.volume = 0.0
        self.symmetry_x = []
        self.symmetry_y = []
        self.symmetry_z = []
        
    def init(self):
        k = math.pi/180.0
        self.volume = self.a*self.b*self.c*math.sqrt(
            1 - math.cos(k*self.alpha)**2- math.cos(k*self.beta)**2 - math.cos(k*self.gamma)**2
            + 2 * math.cos(k*self.alpha) * math.cos(k*self.beta) * math.cos(k*self.gamma))
        self.to_cartesian = np.zeros((3,3))
        self.to_cartesian[0,0] = self.a
        self.to_cartesian[0,1] = self.b*math.cos(k*self.gamma)
        self.to_cartesian[0,2] = self.c*math.cos(k*self.beta)
        self.to_cartesian[1,1] = self.b*math.sin(k*self.gamma)
        self.to_cartesian[1,2] = self.c*(math.cos(k*self.alpha)-math.cos(k*self.beta)*math.cos(k*self.gamma))/math.sin(k*self.gamma)
        self.to_cartesian[2,2] = self.volume / (self.a*self.b*math.sin(k*self.gamma))
        self.to_internal = np.linalg.inv(self.to_cartesian)
        # Generate the atoms coordinate
        temp = []
        for atom in self.internal_atoms:
            for i in range(len(self.symmetry_x)):
                x = atom.copy()
                x.x = eval(self.symmetry_x[i].replace('x',str(atom.x)).replace('y',str(atom.y)).replace('z',str(atom.z)))
                if x.x < 0:
                    x.x += 1
                if x.x > 1:
                    x.x -= 1
                x.y = eval(self.symmetry_y[i].replace('x',str(atom.x)).replace('y',str(atom.y)).replace('z',str(atom.z)))
                if x.y < 0:
                    x.y += 1
                if x.y > 1:
                    x.y -= 1
                x.z = eval(self.symmetry_z[i].replace('x',str(atom.x)).replace('y',str(atom.y)).replace('z',str(atom.z)))
                if x.z < 0:
                    x.z += 1
                if x.z > 1:
                    x.z -= 1
                temp.append(x)
        # Check the repetation in temp
        self.atoms = []
        for i in range(len(temp)-1,-1,-1):
            rep = False
            atom = temp[i]
            for j in range(i-1):
                r2 = (atom.x - temp[j].x)**2 + (atom.y - temp[j].y)**2 + (atom.z - temp[j].z)**2
                if r2 < 0.0001:
                    rep = True
                    break
            if not rep:
                self.atoms.append(atom)
            del temp[i]
        
    def read_cif(self,file_name):
        f = open(file_name,'r')
        temp = f.readline()
        n = 0
        while temp != '':
            if temp.find('_cell_length_a') != -1:
                self.a = float(temp.split()[1])
            if temp.find('_cell_length_b') != -1:
                self.b = float(temp.split()[1])
            if temp.find('_cell_length_c') != -1:
                self.c = float(temp.split()[1])
            if temp.find('_cell_angle_alpha') != -1:
                self.alpha = float(temp.split()[1])
            if temp.find('_cell_angle_beta') != -1:
                self.beta = float(temp.split()[1])
            if temp.find('_cell_angle_gamma') != -1:
                self.gamma = float(temp.split()[1])
            if temp.find('loop_') != -1:
                n += 1
                if n == 1:
                    temp = f.readline()
                    temp = f.readline()
                    while temp != '\n':
                        tmp = temp[(temp.find("'")+1):temp.rfind("'")]
                        tmp = tmp.split(',')
                        self.symmetry_x.append(tmp[0])
                        self.symmetry_y.append(tmp[1])
                        self.symmetry_z.append(tmp[2])
                        temp = f.readline()
                if n == 2:
                    temp = f.readline()
                    while temp[:1] == '_':
                        temp = f.readline()
                    while temp != '\n':
                        tmp = temp.split()
                        atom = Atom()
                        atom.a_type = tmp[0]
                        atom.x = float(tmp[2])
                        atom.y = float(tmp[3])
                        atom.z = float(tmp[4])
                        atom.charge = float(tmp[5])
                        temp = f.readline()
                        self.internal_atoms.append(atom)
            temp = f.readline()
        f.close()
        
    def check(self,atom):
        coord = [atom.x,atom.y,atom.z]
        new_coord = np.dot(self.to_internal,coord)
        for i in range(3):
            while new_coord[i] < 0:
                new_coord[i] += 1
            while new_coord[i] > 1:
                new_coord[i] -= 1
        new_coord = np.dot(self.to_cartesian,new_coord)
        atom.x = new_coord[0]
        atom.y = new_coord[1]
        atom.z = new_coord[2]


# Class for adsorbents in lattice box
class Adsorbent:
    def __init__(self):
        self.atoms = []
    
    def copy_lattice(self,lattice):
        self.a = lattice.a
        self.b = lattice.b
        self.c = lattice.c
        self.alpha = lattice.alpha
        self.beta = lattice.beta
        self.gamma = lattice.gamma
        self.volume = lattice.volume
        self.to_cartesian = lattice.to_cartesian # Copy reference
        self.to_internal = lattice.to_internal # Copy reference
    
    def check(self,atom):
        coord = [atom.x,atom.y,atom.z]
        new_coord = np.dot(self.to_internal,coord)
        for i in range(3):
            if new_coord[i] < 0:
                new_coord[i] += 1
            if new_coord[i] > 1:
                new_coord[i] -= 1
        new_coord = np.dot(self.to_cartesian,new_coord)
        atom.x = new_coord[0]
        atom.y = new_coord[1]
        atom.z = new_coord[2]


# Lattice for adsorbent in another box
class Box:
    def __init__(self):
        self.atoms = []
        self.side = 0.0
        self.volume = 0.0
    
    def change_volume(self,d_volume = 0.0, d_log_volume = 0.0):
        pass
        
    def check(self,atom):
        while atom.x > self.side:
            atom.x -= self.side
        while atom.x < 0:
            atom.x += self.side
        while atom.y > self.side:
            atom.y -= self.side
        while atom.y < 0:
            atom.y += self.side
        while atom.z > self.side:
            atom.z -= self.side
        while atom.z < 0:
            atom.z += self.side
        
