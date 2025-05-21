import os
import re
import astropy.units as u
from astropy.io import ascii
import numpy as np

# Runs considered to be 'modes'
runnums = [36895,
 54978,
 54978,
 16234,
 38978,
 38978,
 69466,
 18919,
 18919,
 14068,
 14068,
 72100,
 72100,
 66281,
 27597,
 77364,
 26233,
 26233,
 1233,
 1233,
 386,
 386,
 71294,
 57224,
 57224,
 45921,
 45921,
 44417,
 44417,
 77187,
 73564,
 73564,
 74692,
 74692,
 78726,
 84413,
 30426,
 30426,
 8568,
 49150,
 49150,
 6383,
 6383,
 31690,
 92899,
 92899,
 47582,
 12351,
 12351,
 68472,
 68472,
 66576,
 66576,
 65694,
 65694,
 99650,
 8264,
 50438,
 50438,
 57218,
 57218,
 24768,
 24768,
 49182,
 49182,
 74776,
 63869,
 63869,
 44693,
 76877,
 60421,
 60421,
 71427,
 81575,
 81575,
 74862,
 74862,
 56735,
 56735,
 49233,
 49233,
 74640,
 74640,
 23132,
 23132,
 50724,
 78643,
 87055,
 58363,
 24637,
 7119,
 22448,
 22149,
 22149,
 64873,
 64873,
 6472,
 6472,
 44339,
 44339]

for r in runnums:
    root_dir =  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/RestrOMultiN/RunNumber'+str(r)+'/' 

    # Regex pattern to extract numbers
    pattern = re.compile(r'outout_Psurfsubtry(\d+)_Innertry_(\d+)\.out$')
    otherpattern = re.compile(r'outout_subtry(\d+)\.out$')

    # Store matches in a dict: {psurfsubtry_number: list of (innertry_number, filepath)}
    psurf_map = {}
    fis = []

    for dirpath, _, filenames in os.walk(root_dir):
        for fname in filenames:
            match = pattern.search(fname)
            if match:
                psurf_num = int(match.group(1))
                inner_num = int(match.group(2))
                full_path = os.path.join(dirpath, fname)
                psurf_map.setdefault(psurf_num, []).append((inner_num, full_path))

    # Now find the highest Psurfsubtry
    if psurf_map:
        max_psurf = max(psurf_map)
        # Within that, find the highest Innertry
        max_inner, max_file = max(psurf_map[max_psurf], key=lambda x: x[0])
        #print(f"Highest file: {max_file}")
        #print(f"Psurfsubtry: {max_psurf}, Innertry: {max_inner}")
        #print('\n')
        fis.append(max_file)
    else:
        print("No matching files found.")
        print(root_dir)
        print('\n')

R_p = 1.097*u.Rearth
fac = 4*np.pi*(R_p**2)

oflux = []
o2flux = []
for fi in fis:

    fio = open(fi, 'r')
    lines = fio.readlines()
    fio.close()

    for l in range(len(lines)):
        if len(lines[l].split('FLUXES OF')) > 1:
            tab = ascii.read(fi, header_start=l-45, data_end=l+58) 

    flo2 = tab['O2'][len(tab['O2'])-1]
    flo = tab['O'][len(tab['O'])-1]

    flo2 = ((flo2*(u.cm**-2 * u.s**-1))*fac).to(u.s**-1).value
    flo = ((flo*(u.cm**-2 * u.s**-1))*fac).to(u.s**-1).value

    oflux.append(flo)
    o2flux.append(flo2)

oflux = np.array(oflux)
o2flux = np.array(o2flux)

np.save('/gscratch/vsm/gialluca/VPLModelingTools_Dev/RestrOMultiN/OxyFluxesTOA.npy', np.array(oflux, o2flux))




        

'''
Think it's this to get the fluxes tables: 

for l in range(len(lines)):
    if len(lines[l].split('FLUXES OF')) > 1:
        tab = ascii.read('outout_subtry1.out', header_start=l-45, data_end=l+58)

fl = tab['O2'][len(tab['O2'])-1]
'''