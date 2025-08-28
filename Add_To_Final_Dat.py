import numpy as np
import astropy.io as ascii
import os
import json

### This is to compile all atmospheric spectra from a planets converged set in a database file 
# atm_type refers to the type of atmospheric composition ('H2O-O2', 'H2O-O2-CO2', 'H2O-O2-SO2', 'H2O-O2-CO2-SO2')
def add_spectra(planet='T1b', atm_type='H2O-O2', sweep_dir=None):

    if os.path.exists('/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/'+planet+'_FinalSpectra_Database.json'):
        f = open('/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/'+planet+'_FinalSpectra_Database.json', 'r')
        fj = json.load(f)
        fj = json.loads(fj)
        f.close()
    
    else:
        fj = {}
    
    # Make sure atm type already exists or add it
    # Set the first ID number to add (arbitrarily assigned)
    if atm_type not in fj.keys():
        fj[atm_type] = {}
        fj[atm_type]['AtmIDs'] = [] # ID numbers, arbitrarily assigned
        curr_id = 1
    
    else:
        curr_id = fj[atm_type]['AtmIDs'][len(fj[atm_type]['AtmIDs'])-1] + 1

    # Read in sweep dir compile table
    partab = ascii.read(sweep_dir+'ParameterSweep_RunStats_failedrun.dat', format='fixed_width', delimiter=' ')
    
    for i, status in enumerate(partab['FinalState']):
        if status == 'Converged':

            # For the particular atmosphere type, metadata is an array with all the source/sink options
            if atm_type == 'H2O-O2':
                metadat = [float(partab['H2O_OutgassRate'][i]), float(partab['O_EscapeRate'][i]), 
                           float(partab['O2_EscapeRate'][i]), float(partab['O3_EscapeRate'][i]),
                           float(partab['H2O2_EscapeRate'][i])]
                
            # Check if this atmosphere is already in the database
            add_to_db = True
            for key in fj[atm_type].keys():
                if key != 'AtmIDs':
                    if fj[atm_type][key]['MetaData'] == metadat:
                        add_to_db = False

            # If you need to add it to DB, do so
            if add_to_db == True:

                # Append current ID to the AtmIDs list
                fj[atm_type]['AtmIDs'].append(curr_id)

                # Add in metadata
                fj[atm_type]['Atm'+str(curr_id)] = {}
                fj[atm_type]['Atm'+str(curr_id)]['MetaData'] = metadat

                # Set path to take data from
                currpath = sweep_dir+partab['ModelNumber']+'/'

                # List surface pressure
                fj[atm_type]['Atm'+str(curr_id)]['SurfPress'] = float(partab['LastPsurf'][i])

                # List tiered convergence
                fj[atm_type]['Atm'+str(curr_id)]['ConvTier'] = partab['2colConv'][i]

                # Read in transmission data from limb 
                trnst = ascii.read(currpath+partab['ModelNumber']+'_SMART.trnst')
                fj[atm_type]['Atm'+str(curr_id)]['Trnst_Wav'] = list(trnst['col1'])
                fj[atm_type]['Atm'+str(curr_id)]['Trnst_Depth'] = list(trnst['col4'])

                # If the planet is interior, need to get the emission data too 
                if planet in ['T1b', 'T1c', 'T1d']:
                    dayside = ascii.read(currpath+partab['ModelNumber']+'_dayside_SMART_toa.rad')
                    fj[atm_type]['Atm'+str(curr_id)]['Dayside_Wav'] = list(dayside['col1'])
                    fj[atm_type]['Atm'+str(curr_id)]['Dayside_Fp'] = list(dayside['col4'])
                    fj[atm_type]['Atm'+str(curr_id)]['Dayside_Fstar'] = list(dayside['col3'])

                    nightside = ascii.read(currpath+partab['ModelNumber']+'_nightside_SMART_toa.rad')
                    fj[atm_type]['Atm'+str(curr_id)]['Nightside_Wav'] = list(nightside['col1'])
                    fj[atm_type]['Atm'+str(curr_id)]['Nightside_Fp'] = list(nightside['col4'])
                    fj[atm_type]['Atm'+str(curr_id)]['Nightside_Fstar'] = list(nightside['col3'])

                # Iterate the current id for the next atmosphere if more are added 
                curr_id += 1

    f = open('/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/'+planet+'_FinalSpectra_Database.json', 'w')
    dh = json.dumps(fj)
    json.dump(dh, f)
    f.close()