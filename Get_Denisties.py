import json
import numpy as np
from astropy.io import ascii

def ingest_outdist(photochem_InputsDir, originalpth):

        # Read in species.dat file
        species = photochem_InputsDir+'species.dat'
        nsp = open(species, 'r')

        # initialize output dictionary
        d = {}

        # Loop through species and find all gas names for naming in dictionary
        done = False
        new_gases = []
        for l in nsp.readlines():
            if l.split()[0][0] == '*':
                if done == True:
                    break
            else:
                new_gases.append(l.split()[0])
                done = True

        # Total number of LL gases (will have mixing ratios in out.dist)
        NQ = len(new_gases)

        # Define relavant trackers for looping through out.dist
        blockstart = 0
        blockend = 200

        NQblocks = np.ceil(NQ/10)
        species_ind = 0

        # Loop through out.dist file NQ blocks and append species mixing ratios to dictionary
        for i in range(int(NQblocks)):
            curr_nq_block = ascii.read(originalpth+'FINAL_out.dist', data_start=blockstart, data_end=blockend, header_start=None)
            blockstart = blockend
            blockend = blockend+200

            for i in curr_nq_block.columns:
                d[new_gases[species_ind]] = list(curr_nq_block[i])
                species_ind += 1

        # Add Temp, Edd, and number density to output dictionary 
        T_edd_block = ascii.read(originalpth+'FINAL_out.dist', data_start=blockstart, data_end=blockend)
        blockstart = blockend

        d['Temp'] = list(T_edd_block['col1'])
        d['Edd'] = list(T_edd_block['col2'])
        d['NDens'] = list(T_edd_block['col3'])

        return d

def get_true_number_densities(indistdic):

    d = {}
    for i in indistdic.keys():
        if i not in ['Temp', 'Edd', 'NDens']:
            newndens = [indistdic['NDens'][lvl]*indistdic[i][lvl] for lvl in range(len(indistdic[i]))]
            d[i] = newndens

    return d

def get_numdens(planet):

    if planet == 'T1b' or planet == 'TRAPPIST-1b':
        finaldb = '../VPLModelingSupportScripts/T1b_FinalSpectra_Database.json'
    elif planet == 'T1c' or planet == 'TRAPPIST-1c':
        finaldb = '../VPLModelingSupportScripts/T1c_FinalSpectra_Database.json'
    elif planet == 'T1d' or planet == 'TRAPPIST-1d':
        finaldb = '../VPLModelingSupportScripts/T1d_FinalSpectra_Database.json'
    elif planet == 'T1e' or planet == 'TRAPPIST-1e':
        finaldb = '../VPLModelingSupportScripts/T1e_FinalSpectra_Database.json'
    elif planet == 'T1f' or planet == 'TRAPPIST-1f':
        finaldb = '../VPLModelingSupportScripts/T1f_FinalSpectra_Database.json'
    elif planet == 'T1g' or planet == 'TRAPPIST-1g':
        finaldb = '../VPLModelingSupportScripts/T1g_FinalSpectra_Database.json'
    elif planet == 'T1h' or planet == 'TRAPPIST-1h':
        finaldb = '../VPLModelingSupportScripts/T1h_FinalSpectra_Database.json'

    # Read in the final DB
    f = open(finaldb, 'r')
    fj = json.load(f)
    fj = json.loads(fj)
    fh2o = fj['H2O-O2']
    fco2 = fj['CO2']
    if planet not in ['T1f', 'T1g', 'T1h']:
        fso2h2o = fj['SO2-H2O']
    fso2co2 = fj['SO2-CO2']
    f.close()

    d = {}
    d['H2O-O2'] = {}
    d['CO2'] = {}
    d['SO2-H2O'] = {}
    d['SO2-CO2'] = {}

    for atm in fh2o['AtmIDs']:

        path = fh2o['Atm'+str(atm)]['OriginalPath']
        photochempath = path+'PhotochemInputs/'
        outdist = ingest_outdist(photochempath, path)

        speciesndens = get_true_number_densities(outdist)

        d['H2O-O2']['Atm'+str(atm)] = speciesndens

    for atm in fco2['AtmIDs']:

        path = fco2['Atm'+str(atm)]['OriginalPath']
        photochempath = path+'PhotochemInputs/'
        outdist = ingest_outdist(photochempath, path)

        speciesndens = get_true_number_densities(outdist)

        d['CO2']['Atm'+str(atm)] = speciesndens

    if planet not in ['T1f', 'T1g', 'T1h']:
        for atm in fso2h2o['AtmIDs']:

            path = fso2h2o['Atm'+str(atm)]['OriginalPath']
            photochempath = path+'PhotochemInputs/'
            outdist = ingest_outdist(photochempath, path)

            speciesndens = get_true_number_densities(outdist)

            d['SO2-H2O']['Atm'+str(atm)] = speciesndens

    for atm in fso2co2['AtmIDs']:

        path = fso2co2['Atm'+str(atm)]['OriginalPath']
        photochempath = path+'PhotochemInputs/'
        outdist = ingest_outdist(photochempath, path)

        speciesndens = get_true_number_densities(outdist)

        d['SO2-CO2']['Atm'+str(atm)] = speciesndens

    f = open('/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/'+planet+'_SpeciesNumberDensities.json', 'w')
    dh = json.dumps(d)
    json.dump(dh, f)
    f.close()