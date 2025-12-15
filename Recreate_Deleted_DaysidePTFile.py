import numpy as np
import pandas as pd
import shutil
from astropy.io import ascii

def get_final_2column_climate_output_temp_profile(outputfile):
    # Define python dictionary to compile data in
    dat = {}
    dat['Atm_Levels'] = 50

    # Column names used in VPL Climate output run
    colnames = ['P[Pa]', 'Alt[km]', 'T[K]', 'Q_s[K/day]', 'Q_t[K/day]', 'Q_c[K/day]', 
                'Q_ad[K/day]', 'Q_net[K/day]', 'fs_net[W/m/m]', 'ft_net[W/m/m]', 'fc[W/m/m]', 
                'f_ad[W/m/m]', 'pc[Pas]', 'Altc[km]', 'Tc[K]', 'dt[s]', 'lr[K/km]', 
                'aid_lr[K/km]', 'Km[m2/s]', 'rmix[kg/kg]']

    # Open a simple text instance of the Climate output to use for parsing
    dat['FileName'] = outputfile
    fop = open(outputfile, 'r')

    flines = fop.readlines()
    fop.close()

    # Loop through text instance of output file to retrieve atmospheric profiles
    # Start at bottom of file to extract only the last profile
    nightside_done = False
    dat['Nightside'] = {}
    dat['Dayside'] = {}
    for i in reversed(range(len(flines))):
        curr_line = flines[i].split()
        if len(curr_line) > 0:
            if curr_line[0] == '(Pas)': # This checks if you're at a profile
                # This reads in that profile beautifully as pandas data frame
                curr_step = pd.read_csv(dat['FileName'], delimiter=' ', skipinitialspace=True, header=0, 
                                        names=colnames, skiprows=i, nrows=50)
                
                if nightside_done == False:
                    # Add that profile to the dictionary
                    for k in colnames:
                        dat['Nightside'][k] = np.array(curr_step[k])

                    # Find net flux at each level
                    Fnet = np.zeros(len(curr_step['fs_net[W/m/m]']))
                    for lvl in range(len(Fnet)):
                        Fnet[lvl] = curr_step['fs_net[W/m/m]'][lvl] - curr_step['ft_net[W/m/m]'][lvl] - curr_step['fc[W/m/m]'][lvl]
                    dat['Nightside']['f_net[W/m/m]'] = Fnet

                    nightside_done = True
                
                else:
                    # Add that profile to the dictionary
                    for k in colnames:
                        dat['Dayside'][k] = np.array(curr_step[k])

                    # Find net flux at each level
                    Fnet = np.zeros(len(curr_step['fs_net[W/m/m]']))
                    for lvl in range(len(Fnet)):
                        Fnet[lvl] = curr_step['fs_net[W/m/m]'][lvl] - curr_step['ft_net[W/m/m]'][lvl] - curr_step['fc[W/m/m]'][lvl]
                    dat['Dayside']['f_net[W/m/m]'] = Fnet
                    break

    return dat

def recreate_PT_Dayside_File():

    # Read in the file of input options that need extant PT dayside profiles:
    options = pd.read_csv('/gscratch/vsm/gialluca/VPLModelingTools_Dev/Cinit/VMRSSurfP_RatesInSweep_ForFutureInputOptions.dat', delimiter=' ', index_col=['ModelNumber'])
    models = list(options.index)

    for m in models:

        # pick 2 column output file 
        runnum = m.split('/')
        runnum = runnum[len(runnum)-1]
        curr_outputfi = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+m+'/vpl_2col_climate_output_'+str(runnum)+'.run'

        # Retrieve profiles:
        dat = get_final_2column_climate_output_temp_profile(curr_outputfi)

        # Copy nightside PT to replace temp col:
        shutil.copyfile('/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+m+'/PT_profile_nightside_'+str(runnum)+'.pt', '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+m+'/PT_profile_dayside_'+str(runnum)+'.pt')
        daysidept = ascii.read('/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+m+'/PT_profile_dayside_'+str(runnum)+'.pt')

        daysidept['Temp'] = dat['Dayside']['T[K]']
        ascii.write(daysidept, '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+m+'/PT_profile_dayside_'+str(runnum)+'.pt', overwrite=True)