'''
Monte Carlo code
written by: Thien-Phuc Tu-Nguyen
Last modified: 2017
File: forcefield
'''

import math
import numpy as np

class ForceField:
    def __init__(self):
        self.raw_sigma = []
        self.raw_epsilon = []
        self.sigma2 = None
        self.epsilon4 = None
        self.p_type = []
        self.p_type_missing = []
        
    def init(self):
        n = len(self.p_type)
        self.sigma2 = np.zeros((n,n))
        self.epsilon4 = np.zeros((n,n))
        for i in range(n-1):
            for j in range(i+1,n):
                self.sigma2[i,j] = (self.raw_sigma[i] + self.raw_sigma[j])**2/4
                self.sigma2[j,i] = self.sigma2[i,j]
                self.epsilon4[i,j] = 4*math.sqrt(self.raw_epsilon[i]*self.raw_epsilon[j])
                self.epsilon4[j,i] = self.epsilon4[i,j]
        for i in range(n):
            self.sigma2[i,i] = self.raw_sigma[i] ** 2
            self.epsilon4[i,i] = 4 * self.raw_epsilon[i]
        
    def read_raspa_def(self,file_name):
        f = open(file_name,'r')
        temp = None
        for i in range(5):
            temp = f.readline()
        n = int(f.readline())
        temp = f.readline()
        for i in range(n):
            temp = f.readline()
            temp = temp.split()
            if len(temp) > 2:
                self.p_type.append(temp[0])
                self.raw_epsilon.append(float(temp[2]))
                self.raw_sigma.append(float(temp[3]))
            else:
                self.p_type_missing.append(temp[0])
        f.close()
        
    def set_index(self,container): # Set a_index
        for a in container.atoms:
            a.p_index = self.p_type.index(a.a_type) # This way does not work out
        
    def pair(self,a,b):
        r2 = (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2
        if r2 > 0.00001:
            r2 = self.sigma2[a.p_index,b.p_index] / r2
            en = self.epsilon4[a.p_index,b.p_index] * (r2 ** 6 - r2 ** 3)
            en += a.charge * b.charge / r2 ** 0.5 # Energy conversion missing
            return en
        else:
            return 0.0
        
    def one_atom(self,atom, container):
        en = 0.0
        for x in container.atoms:
            en += self.pair(atom,x)
        return en # Did not count periodic
        
    def box(self,box):
        en = 0.0
        n = len(box.atoms)
        for i in range(n-1):
            for j in range(i+1,n):
                en += self.pair(box.atoms[i],box.atoms[j])
        return en # Did not count periodic


