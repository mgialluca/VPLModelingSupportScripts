import numpy as np
from astropy.io import ascii
import os
import json
import astropy.units as u

### This is to compile all atmospheric spectra from a planets converged set in a database file 
# atm_type refers to the type of atmospheric composition ('H2O-O2', 'H2O-O2-CO2', 'H2O-O2-SO2', 'H2O-O2-CO2-SO2')
def add_spectra(planet='T1b', atm_type='H2O-O2', sweep_dir=None):

    if planet == 'T1b':
        Rp = 1.116*u.Rearth
    elif planet == 'T1c':
        Rp = 1.097*u.Rearth
    elif planet == 'T1d':
        Rp = 0.788*u.Rearth
    elif planet == 'T1e':
        Rp = 0.920*u.Rearth
    elif planet == 'T1f':
        Rp = 1.045*u.Rearth
    elif planet == 'T1g':
        Rp = 1.129*u.Rearth
    elif planet == 'T1h':
        Rp = 0.755*u.Rearth

    if atm_type in ['SO2-H2O', 'SO2-CO2']:
        fso2 = open('/gscratch/vsm/gialluca/VPLModelingTools_Dev/AddSO2/'+planet+'_Tracking.json', 'r')
        fjso2 = json.load(fso2)
        fjso2 = json.loads(fjso2)
        fso2.close()

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
                
            elif atm_type == 'CO2':
                metadat = [float(partab['H2O_OutgassRate'][i]), float(partab['CO2_OutgassRate'][i]), float(partab['O_EscapeRate'][i]), 
                           float(partab['O2_EscapeRate'][i]), float(partab['O3_EscapeRate'][i]),
                           float(partab['H2O2_EscapeRate'][i]), float(partab['CO_EscapeRate'][i]), float(partab['CO2_EscapeRate'][i])]
                
            elif atm_type == 'SO2-H2O':
                metadat = [float(partab['H2O_OutgassRate'][i]), float(partab['SO2_OutgassRate']), float(partab['O_EscapeRate'][i]), 
                           float(partab['O2_EscapeRate'][i]), float(partab['O3_EscapeRate'][i]),
                           float(partab['H2O2_EscapeRate'][i])]
            
            elif atm_type == 'SO2-CO2':
                metadat = [float(partab['H2O_OutgassRate'][i]), float(partab['CO2_OutgassRate'][i]), float(partab['SO2_OutgassRate']), 
                           float(partab['O_EscapeRate'][i]), float(partab['O2_EscapeRate'][i]), float(partab['O3_EscapeRate'][i]),
                           float(partab['H2O2_EscapeRate'][i]), float(partab['CO_EscapeRate'][i]), float(partab['CO2_EscapeRate'][i])]
                
            # Check if this atmosphere is already in the database
            add_to_db = True
            for key in fj[atm_type].keys():
                if key != 'AtmIDs':
                    if fj[atm_type][key]['MetaData'] == metadat:
                        add_to_db = False

            # Need to exclude bad atmospheres:
            if add_to_db == True:
                if partab['ModelNumber'][i][0] == 'u':
                    modnum = 'R'+partab['ModelNumber'][i]
                else:
                    modnum = partab['ModelNumber'][i]
                
                currpath = sweep_dir+modnum+'/'
                ptz = ascii.read(currpath+'FINAL_PTZ_mixingratios_out.dist')
                if atm_type == 'H2O-O2':
                    if (ptz['O'][0] + ptz['O2'][0] + ptz['H2O'][0] + ptz['O3'][0]) < 0.9:
                        add_to_db = False
                elif atm_type == 'CO2':
                    if ptz['C'][0] > 0.05:
                        add_to_db = False


            # If you need to add it to DB, do so
            if add_to_db == True:

                # Append current ID to the AtmIDs list
                fj[atm_type]['AtmIDs'].append(curr_id)

                # Add in metadata
                fj[atm_type]['Atm'+str(curr_id)] = {}
                fj[atm_type]['Atm'+str(curr_id)]['MetaData'] = metadat

                # Set path to take data from
                currpath = sweep_dir+modnum+'/'

                # Save path in case
                fj[atm_type]['Atm'+str(curr_id)]['OriginalPath'] = currpath

                # If this is an SO2 atmosphere, need to sort out the original stable atmosphere
                if atm_type in ['SO2-H2O', 'SO2-CO2']:
                    ogsweep = sweep_dir.split('/')[len(sweep_dir.split('/'))-2]
                    modnum_num = modnum.split('RunNumber')[1]
                    copiedatm = fjso2[ogsweep]['Atm'+str(modnum_num)]['Inputs'][0]
                    so2amount = fjso2[ogsweep]['Atm'+str(modnum_num)]['SO2Amount']
                    fj[atm_type]['Atm'+str(curr_id)]['CopiedStableAtm'] = copiedatm
                    fj[atm_type]['Atm'+str(curr_id)]['SO2Amount'] = so2amount

                # List surface pressure
                fj[atm_type]['Atm'+str(curr_id)]['SurfPress'] = float(partab['LastPsurf'][i])

                # List tiered convergence
                fj[atm_type]['Atm'+str(curr_id)]['ConvTier'] = partab['2colConv'][i]

                # Read in transmission data from limb 
                trnst = ascii.read(currpath+modnum+'_SMART.trnst')
                fj[atm_type]['Atm'+str(curr_id)]['Trnst_Wav'] = list(trnst['col1'])
                fj[atm_type]['Atm'+str(curr_id)]['Trnst_Depth'] = list(trnst['col4'])

                # If the planet is interior, need to get the emission data too 
                if planet in ['T1b', 'T1c', 'T1d']:
                    dayside = ascii.read(currpath+modnum+'_dayside_SMART_toa.rad')
                    fj[atm_type]['Atm'+str(curr_id)]['Dayside_Wav'] = list(dayside['col1'])
                    fj[atm_type]['Atm'+str(curr_id)]['Dayside_Fp'] = list(dayside['col4'])
                    fj[atm_type]['Atm'+str(curr_id)]['Dayside_Fstar'] = list(dayside['col3'])

                    nightside = ascii.read(currpath+modnum+'_nightside_SMART_toa.rad')
                    fj[atm_type]['Atm'+str(curr_id)]['Nightside_Wav'] = list(nightside['col1'])
                    fj[atm_type]['Atm'+str(curr_id)]['Nightside_Fp'] = list(nightside['col4'])
                    fj[atm_type]['Atm'+str(curr_id)]['Nightside_Fstar'] = list(nightside['col3'])

                # Get the PTZ mixing ratios file in the database 
                #ptz = ascii.read(currpath+'FINAL_PTZ_mixingratios_out.dist') # Now read in above
                fj[atm_type]['Atm'+str(curr_id)]['PTZOut'] = {}
                for col in ptz.colnames:
                    fj[atm_type]['Atm'+str(curr_id)]['PTZOut'][col] = list(ptz[col])

                # Get the condensation and escape rates (fluxes)
                fout = open(currpath+'FINAL_out.out', 'r')
                lines = fout.readlines()
                names = ['Z', 'T', 'EDD', 'DEN', 'P', 'H2OSAT', 'H2O', 'RELH', 'CONDEN', 'HFLUX']

                for l in range(len(lines)):
                    '''if 'CONDEN' in lines[l].split():
                        tab = ascii.read(currpath+'FINAL_out.out', data_start=l-83, data_end=l+17, names=names)
                        fj[atm_type]['Atm'+str(curr_id)]['P_CONDEN'] = list(tab['P'])
                        fj[atm_type]['Atm'+str(curr_id)]['CONDEN'] = list(tab['CONDEN'])
                        fj[atm_type]['Atm'+str(curr_id)]['H2O_CONDEN'] = list(tab['H2O'])'''
                    
                    if 'FLUXES' in lines[l].split() and 'ENERGY' not in lines[l].split():
                        hold = lines[l+105].split()
                        oesc = float(hold[1])*(1/(u.cm**2 * u.s))
                        o2esc = float(hold[2])*(1/(u.cm**2 * u.s))
                        oesc = oesc*(4*np.pi*(Rp**2))
                        o2esc = o2esc*(4*np.pi*(Rp**2))
                        oesc = oesc.to(1/u.s).value
                        o2esc = o2esc.to(1/u.s).value

                        fj[atm_type]['Atm'+str(curr_id)]['Oesc_per-s'] = oesc
                        fj[atm_type]['Atm'+str(curr_id)]['O2esc_per-s'] = o2esc

                        '''if atm_type == 'CO2':
                            assert lines[l+211].split()[9] == 'CO2'
                            hold = lines[l+313].split()
                            co2esc = float(hold[9])*(1/(u.cm**2 * u.s))
                            co2esc = co2esc*(4*np.pi*(Rp**2))
                            co2esc = co2esc.to(1/u.s).value
                            fj[atm_type]['Atm'+str(curr_id)]['CO2esc_per-s'] = co2esc'''

                # Iterate the current id for the next atmosphere if more are added 
                curr_id += 1

    f = open('/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/'+planet+'_FinalSpectra_Database.json', 'w')
    dh = json.dumps(fj)
    json.dump(dh, f)
    f.close()