import os
import re
from astropy.io import ascii
import json

# Define the directory containing the files
directory = '../TestPress2/RunNumber7/'

# Get all files in the directory
files = os.listdir(directory)

# Define a regex pattern to extract the numbers after 'Psurfsubtry' and 'Innertry'
pattern = r"PTZMix_Psurfsubtry(\d+)_Innertry_(\d+).run"

# Create a list to store tuples of the form (Psurfsubtry, Innertry, file_path)
file_data = []

# Loop through the files and apply the regex to extract Psurfsubtry and Innertry
for file in files:
    match = re.search(pattern, file)
    if match:
        Psurfsubtry = int(match.group(1))  # Extract the Psurfsubtry number
        Innertry = int(match.group(2))    # Extract the Innertry number
        file_data.append((Psurfsubtry, Innertry, file))

# Sort the files first by Psurfsubtry, then by Innertry
file_data.sort(key=lambda x: (x[0], x[1]))

# Now loop through the sorted files and perform your operation
d = {}
currPsurfSubtry = 1
for Psurfsubtry, Innertry, file in file_data:
    # Perform the operation on each file
    # print(f"Processing file: {file} (Psurfsubtry: {Psurfsubtry}, Innertry: {Innertry})")
    
    if Psurfsubtry > currPsurfSubtry:
        d['PTry'+str(currPsurfSubtry)]['toa_press'] = toa_press
        d['PTry'+str(currPsurfSubtry)]['toa_alt'] = toa_alt
        d['PTry'+str(currPsurfSubtry)]['toa_o'] = toa_o
        d['PTry'+str(currPsurfSubtry)]['toa_o2'] = toa_o2
        d['PTry'+str(currPsurfSubtry)]['toa_o3'] = toa_o3
        d['PTry'+str(currPsurfSubtry)]['toa_h2o'] = toa_h2o
        d['PTry'+str(currPsurfSubtry)]['toa_h'] = toa_h

        d['PTry'+str(currPsurfSubtry)]['boa_press'] = boa_press
        d['PTry'+str(currPsurfSubtry)]['boa_alt'] = boa_alt
        d['PTry'+str(currPsurfSubtry)]['boa_o'] = boa_o
        d['PTry'+str(currPsurfSubtry)]['boa_o2'] = boa_o2
        d['PTry'+str(currPsurfSubtry)]['boa_o3'] = boa_o3
        d['PTry'+str(currPsurfSubtry)]['boa_h2o'] = boa_h2o
        d['PTry'+str(currPsurfSubtry)]['boa_h'] = boa_h
        
        print('Convergence found, moved to next pressure iteration\n')
    
    currPsurfSubtry = Psurfsubtry

    ptz = ascii.read(file)
    boa = len(ptz)-1

    print(f"Psurfsubtry: {Psurfsubtry}, Innertry: {Innertry}")
    print('--------------------')
    print('Top of Atmosphere: ')
    print('Pressure - ', ptz['PRESS'][0])
    print('Altitude - ', ptz['ALT'][0])
    print('O Mixing - ', ptz['O'][0])
    print('O2 Mixing - ', ptz['O2'][0])
    print('O3 Mixing - ', ptz['O3'][0])
    print('H2O Mixing - ', ptz['H2O'][0])
    print('H Mixing - ', ptz['H'][0])
    print('--------------------')
    print('Bottom of Atmosphere: ')
    print('Pressure - ', ptz['PRESS'][boa])
    print('Altitude - ', ptz['ALT'][boa])
    print('O Mixing - ', ptz['O'][boa])
    print('O2 Mixing - ', ptz['O2'][boa])
    print('O3 Mixing - ', ptz['O3'][boa])
    print('H2O Mixing - ', ptz['H2O'][boa])
    print('H Mixing - ', ptz['H'][boa])
    print('\n')

    toa_press = ptz['PRESS'][0]
    toa_alt = ptz['ALT'][0]
    toa_o = ptz['O'][0]
    toa_o2 = ptz['O2'][0]
    toa_o3 = ptz['O3'][0]
    toa_h2o = ptz['H2O'][0]
    toa_h = ptz['H'][0]
    boa_press = ptz['PRESS'][boa]
    boa_alt = ptz['ALT'][boa]
    boa_o = ptz['O'][boa]
    boa_o2 = ptz['O2'][boa]
    boa_o3 = ptz['O3'][boa]
    boa_h2o = ptz['H2O'][boa]
    boa_h = ptz['H'][boa]


d['PTry'+str(currPsurfSubtry)]['toa_press'] = toa_press
d['PTry'+str(currPsurfSubtry)]['toa_alt'] = toa_alt
d['PTry'+str(currPsurfSubtry)]['toa_o'] = toa_o
d['PTry'+str(currPsurfSubtry)]['toa_o2'] = toa_o2
d['PTry'+str(currPsurfSubtry)]['toa_o3'] = toa_o3
d['PTry'+str(currPsurfSubtry)]['toa_h2o'] = toa_h2o
d['PTry'+str(currPsurfSubtry)]['toa_h'] = toa_h

d['PTry'+str(currPsurfSubtry)]['boa_press'] = boa_press
d['PTry'+str(currPsurfSubtry)]['boa_alt'] = boa_alt
d['PTry'+str(currPsurfSubtry)]['boa_o'] = boa_o
d['PTry'+str(currPsurfSubtry)]['boa_o2'] = boa_o2
d['PTry'+str(currPsurfSubtry)]['boa_o3'] = boa_o3
d['PTry'+str(currPsurfSubtry)]['boa_h2o'] = boa_h2o
d['PTry'+str(currPsurfSubtry)]['boa_h'] = boa_h

fout = open('PressureIterations.json', 'w')
dh = json.dumps(d)
json.dump(dh, fout)
fout.close()