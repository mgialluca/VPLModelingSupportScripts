import os
from astropy.io import ascii
import numpy as np 
import astropy.units as u
import json

def get_conden_esc_flux(sweepdir, planet='T1b'):

    if planet == 'T1b':
        R_p = 1.116*u.Rearth
    elif planet == 'T1c':
        R_p = 1.097*u.Rearth # Currently radius of Trappist-1c 
    elif planet == 'T1d':
        R_p = 0.788*u.Rearth
    elif planet == 'T1e':
        R_p = 0.920*u.Rearth
    elif planet == 'T1f':
        R_p = 1.045*u.Rearth
    elif planet == 'T1g':
        R_p = 1.129*u.Rearth
    elif planet == 'T1h':
        R_p = 0.755*u.Rearth

    d = {}

    for ds, sds, fis in os.walk(sweepdir):
        break 

    names = ['Z', 'T', 'EDD', 'DEN', 'P', 'H2OSAT', 'H2O', 'RELH', 'CONDEN', 'HFLUX']

    for run in sds:
        for dsh, sdsh, fish in os.walk(sweepdir+run):
            break
        if 'FINAL_out.out' in fish:
            d[run] = {}

            f = open(sweepdir+run+'/FINAL_out.out', 'r')
            lines = f.readlines()

            for l in range(len(lines)):
                if 'CONDEN' in lines[l].split():
                    tab = ascii.read(sweepdir+run+'/FINAL_out.out', data_start=l-83, data_end=l+17, names=names)
                    d[run]['P'] = list(tab['P'])
                    d[run]['CONDEN'] = list(tab['CONDEN'])
                    d[run]['H2O'] = list(tab['H2O'])
                
                if 'FLUXES' in lines[l].split() and 'ENERGY' not in lines[l].split():
                    hold = lines[l+105].split()
                    oesc = float(hold[1])*(1/(u.cm**2 * u.s))
                    o2esc = float(hold[2])*(1/(u.cm**2 * u.s))
                    oesc = oesc*(4*np.pi*(R_p**2))
                    o2esc = o2esc*(4*np.pi*(R_p**2))
                    oesc = oesc.to(1/u.s).value
                    o2esc = o2esc.to(1/u.s).value

                    d[run]['Oesc'] = oesc
                    d[run]['O2esc'] = o2esc

    f = open(sweepdir+planet+'_conden_escapefluxes.json', 'w')
    dh = json.dumps(d)
    json.dump(dh, f)
    f.close()
    
    return d