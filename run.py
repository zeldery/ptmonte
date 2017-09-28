''' Running program for adsorption
'''

from ptmonte import GrandCanonicalSimulation

if __name__ == '__main__':
    x = GrandCanonicalSimulation()
    x.temperature = 273.15+25
    x.pressure = 1e5
    x.d_max = 5
    x.mass = 16
    x.p_step = [0.4,0.3,0.3]
    x.init('IRMOF-1.cif','force_field_mixing_rules.def')
    x.run(500)
