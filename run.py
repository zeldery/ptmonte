''' Running program for adsorption
'''

from ptmonte import GrandCanonicalSimulation

if __name__ == '__main__':
    x = GrandCanonicalSimulation()
    x.temperature = 273.15+25
    x.pressure = 3e6
    x.d_max = 2.0
    x.mass = 16
    x.p_step = [0.5,0.25,0.25]
    x.init('IRMOF-1.cif','force_field_mixing_rules.def','CH4_sp3')
    x.run(50000)
    x.reset()
    x.run(50000)
