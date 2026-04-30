# Need to be able to get on Klone to do this
from astropy.io import ascii 
import json 
import os

planets = ['b', 'c', 'd', 'e', 'f', 'g', 'h']

b_suites = ['Bco2', 'Bco2wL', 'Bco2wL2', 'Bco2wL3', 'Bz1Re', 'Bzoom1', 'Bzoom2', 'AdjSO2/Bco2', 'AdjSO2/Bz1Re', 'AdjSO2/Bzoom1', 'AdjSO2/Bzoom2']
c_suites = ['Cinit', 'Czoom1', 'Czoom2', 'CzoomRT', 'Cco2', 'Cco2wL', 'Cco2wL2', 'AdjSO2/Cinit', 'AdjSO2/Czoom1', 'AdjSO2/Czoom2', 'AdjSO2/CzoomRT', 'AdjSO2/Cco2', 'AdjSO2/Cco2wL', 'AdjSO2/Cco2wL2']
d_suites = ['Dinit', 'Dzoom1', 'Dzoom2', 'Dco2', 'Dco2rego', 'Dco2wL', 'Dco2wL2', 'AdjSO2/Dinit', 'AdjSO2/Dzoom1', 'AdjSO2/Dzoom2', 'AdjSO2/Dco2', 'AdjSO2/Dco2rego', 'AdjSO2/Dco2wL', 'AdjSO2/Dco2wL2']
e_suites = ['Einit', 'Ezoom1', 'Ezoom2', 'Ezoom3', 'Eco2', 'Eco2lo', 'Eco2wL', 'Eco2wL2', 'AdjSO2/Einit', 'AdjSO2/Ezoom1', 'AdjSO2/Ezoom2', 'AdjSO2/Ezoom3', 'AdjSO2/Eco2', 'AdjSO2/Eco2wL', 'AdjSO2/Eco2wL2']
f_suites = ['Fco2', 'Fco2lo', 'Fco2wL', 'Fco2wL2', 'AdjSO2/Fco2', 'AdjSO2/Fco2wL', 'AdjSO2/Fco2wL2']
g_suites = ['Gco2', 'Gco2lo', 'Gco2wL', 'Gco2wL2', 'AdjSO2/Gco2lo', 'AdjSO2/Gco2wL', 'AdjSO2/Gco2wL2']
h_suites = ['Hco2', 'Hco2lo', 'Hco2wL', 'AdjSO2/Hco2', 'AdjSO2/Hco2lo', 'AdjSO2/Hco2wL']

suites = [b_suites, c_suites, d_suites, e_suites, f_suites, g_suites, h_suites]

path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'
tablefile = 'ParameterSweep_RunStats_failedrun.dat'

output = {}
for sdat, p in zip(suites, planets):

    output[p] = {}
    output[p]['FailPressures'] = []
    output[p]['Unknown'] = 0

    for s in sdat:
        try:
            tab = ascii.read(path+s+'/'+tablefile)
        except:
            tab = ascii.read(path+s+'/'+tablefile, delimiter=' ', format='fixed_width')
        
        for i in range(len(tab)):
            if tab[i]['FinalState'] != 'Converged':
                rn = tab[i]['ModelNumber']
                if rn[0] == 'u':
                    rn = 'R'+rn
                for ds, sds, fs in os.walk(path+s+'/'+rn):
                    break

                if 'FINAL_PTZ_mixingratios_out_FAILED.dist' in fs:

                    ptz = ascii.read(path+s+'/'+rn+'/FINAL_PTZ_mixingratios_out_FAILED.dist')
                    output[p]['FailPressures'].append(ptz['PRESS'][0])

                elif 'PT_profile_'+rn+'.pt' in fs:

                    pt = ascii.read(path+s+'/'+rn+'/PT_profile_'+rn+'.pt')
                    output[p]['FailPressures'].append(pt['Press'][0])

                else:
                    output[p]['Unknown'] += 1


fout = open('/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/FailFinalPressures.json', 'w')
dh = json.dumps(output)
json.dump(dh, fout)
fout.close()