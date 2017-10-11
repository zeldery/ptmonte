''' Running program for adsorption
'''

from ptmonte import GrandCanonicalSimulation
import pandas as pd

if __name__ == '__main__':
    x = GrandCanonicalSimulation()
    x.temperature = 273.15+25
    x.pressure = 5e6
    x.d_max = 0.4
    x.mass = 16
    x.p_step = [0.5,0.25,0.25]
    x.init('IRMOF-1.cif','force_field_mixing_rules.def','CH4_sp3')
    x.run(1000)
    x.reset()
    x.run(1000)
    y = pd.Series(x.record_adsorb)
    print('Mean :{:10.5} , STD :{:10.5}'.format(y.mean(),y.std()/1000**0.5))
