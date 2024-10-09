import numpy as np
import os,sys,shutil,subprocess
from astropy.io import ascii
import sys
import astropy.units as u
import astropy.constants as const
from astropy.table import Table
import json
import pandas as pd

################################
##
## Working on a full semi-automatic VPL modeling pipeline
## Author: Megan Gialluca
##
################################

# Set Global:
photochemDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/'
atmosDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/'
OutPath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/'
AtmProfPath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/AtmProfiles/'

nlevel_fine = 200
nlevel_coarse = 70

######################################################### Support Functions

### Run 1 photochem run
##
## Inputs:
# casename - name of case you're runnin (to name make and run output files)
# CleanMake - do a clean make of photochem or no
# InputCopy - the path to the input photochem files, if false assumes they are already in the PHOTOCHEM/INPUTS directory
# trynum - the iteration number you're on for the specific case
# photochemDir - Path to the PHOTOCHEM dir in atmos (so you can use specific instances)
# atmosDir - Path to the atmos/ dir
# OutPath - path to write model make and run outputs to
##

#def run_photochem_1instance(casename, CleanMake=True, InputCopy=False, trynum=1, 
#                            photochemDir='/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/', 
#                            atmosDir='/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/', 
#                            OutPath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/'):
def run_photochem_1instance(casename, CleanMake=True, 
                            #InputCopy='/home/mgialluca/Nextcloud/VPL_Modeling/Atmos_Dev/atmos/PHOTOCHEM/INPUTFILES/TEMPLATES/ModernEarth/', 
                            #InputCopy='/home/mgialluca/Nextcloud/VPL_Modeling/Atmos_Dev/atmos/PHOTOCHEM/INPUTFILES/TEMPLATES/Venus/',
                            InputCopy='/home/mgialluca/Nextcloud/VPL_Modeling/VPLModelingSupportScripts/Bodies/T1c/T1cOutgas_Testing/',
                            trynum=1, 
                            photochemDir='/home/mgialluca/Nextcloud/VPL_Modeling/Atmos_Dev/atmos/PHOTOCHEM/', 
                            atmosDir='/home/mgialluca/Nextcloud/VPL_Modeling/Atmos_Dev/atmos/', 
                            OutPath='/home/mgialluca/Nextcloud/VPL_Modeling/JWSTObserving_TestOutgassing/'):

    # If you have new input files to use, give 'InputCopy' the dir path
    if InputCopy != False:
        if CleanMake != True:
            choice = input('Youve selected new inputs but do not want a clean make, this is typically incorrect. Do you still want to continue? [y/n] \n')
            if choice == 'n' or choice == 'N':
                sys.exit('Try to Clean Make with new inputs')
            elif choice == 'y' or choice == 'Y':
                pass
            else:
                choice = input('I SAID continue? - """ y """ or """ n """"\n')
                if choice == 'n' or choice == 'N':
                    sys.exit('Try to Clean Make with new inputs')
                elif choice == 'y' or choice == 'Y':
                    pass
                else:
                    sys.exit('User doesnt follow instructions, terminating')

        # Remove old input files if they exist
        subprocess.run('rm '+photochemDir+'INPUTFILES/input_photchem.dat', shell=True)
        subprocess.run('rm '+photochemDir+'INPUTFILES/parameters.inc', shell=True)
        #subprocess.run('rm '+photochemDir+'INPUTFILES/params.dat', shell=True)
        subprocess.run('rm '+photochemDir+'INPUTFILES/PLANET.dat', shell=True)
        subprocess.run('rm '+photochemDir+'INPUTFILES/reactions.rx', shell=True)
        subprocess.run('rm '+photochemDir+'INPUTFILES/species.dat', shell=True)
        subprocess.run('rm '+photochemDir+'in.dist', shell=True)

        # Copy new input files to the right places
        subprocess.run('cp '+InputCopy+'input_photchem.dat '+photochemDir+'INPUTFILES/', shell=True)
        subprocess.run('cp '+InputCopy+'parameters.inc '+photochemDir+'INPUTFILES/', shell=True)
        #subprocess.run('cp '+InputCopy+'params.dat '+photochemDir+'INPUTFILES/', shell=True)
        subprocess.run('cp '+InputCopy+'PLANET.dat '+photochemDir+'INPUTFILES/', shell=True)
        subprocess.run('cp '+InputCopy+'reactions.rx '+photochemDir+'INPUTFILES/', shell=True)
        subprocess.run('cp '+InputCopy+'species.dat '+photochemDir+'INPUTFILES/', shell=True)
        subprocess.run('cp '+InputCopy+'in.dist '+photochemDir, shell=True)

    # Clear the outputs
    subprocess.run('rm -rf '+photochemDir+'OUTPUT/*', shell=True)
    subprocess.run('rm -rf '+photochemDir+'PTZ_mixingratios_in.dist', shell=True)

    # Clean make, if requested
    if CleanMake:
        if trynum == 1:
            fmake = open(OutPath+'photochem_make_output_'+casename+'.txt', 'w')
        else:
            fmake = open(OutPath+'photochem_make_output_'+casename+'_Try'+str(trynum)+'.txt', 'w')
        workdir = os.getcwd()
        os.chdir(atmosDir)
        subprocess.run('make -f ./PhotoMake clean', shell=True)
        subprocess.run('make -f ./PhotoMake', shell=True, stdout=fmake)
        os.chdir(workdir)

    # Run photochem
    if trynum == 1:
        f = open(OutPath+'photochem_run_output_'+casename+'.run', 'w')
    else:
        f = open(OutPath+'photochem_run_output_'+casename+'_Try'+str(trynum)+'.run', 'w')
    workdir = os.getcwd()
    os.chdir(atmosDir)
    subprocess.run('./Photo.run', shell=True, stdout=f)
    os.chdir(workdir)

### Run VPL Climate for a given runscript/executable
##
## Inputs:
# casename - name of case you're running (to name output file)
# runscript - name of VPL Climate runscript WITH PATH
# exec - the VPL Climate executable to use WITH PATH
# trynum - the iteration number youre on for the specific case
# OutPath - path to put run output in
##
def run_climate_1instance(casename, runscript, exec, trynum=1, 
                          OutPath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/'):
    if trynum == 1:
        subprocess.run(exec+' < '+runscript+' > '+OutPath+'vpl_climate_output_'+casename+'.run', shell=True)
    else:
        subprocess.run(exec+' < '+runscript+' > '+OutPath+'vpl_climate_output_'+casename+'_Try'+str(trynum)+'.run', shell=True)

### Run LBLABC for a given runscript (just useful to get the output w/e)
##
## Inputs:
# casename - name of case you're running (to name output file)
# runscript - name of LBLABC runscript WITH PATH
# OutPath - path to put run output in
##
def run_lblabc_1instance(casename, runscript, molecule, 
                         OutPath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/'):
    f = open(OutPath+'lblabc_run_output_'+molecule+'_'+casename+'.run', 'w')
    workdir = os.getcwd()
    os.chdir('/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/')
    subprocess.run('/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/lblabc < '+runscript, shell=True, stdout=f)
    os.chdir(workdir)

### Take the PT profile output from photochem and degrade it to a specified number of layers ...
### ... to create PT profile for LBLABC and SMART
##
## Inputs:
# casename - name of case on your grid you're doing
# nlevel_new - number of layers you want in your degraded atmosphere
# Prof - PTZ_mixingratios_out WITH PATH output by photochem
# PressUnits - can do Bar or Pa, Megan stick with bar for the forseeable forever
# AtmProfPath - path to output new PT profile to
##
def degrade_PT(casename, nlevel_new, photochemDir='/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/', 
               PressUnits='Bar', AtmProfPath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/AtmProfiles/'):
    atm = ascii.read(photochemDir+'OUTPUT/PTZ_mixingratios_out.dist', delimiter=' ')
    alt = atm['ALT']
    pres = atm['PRESS']
    temp = atm['TEMP']

    new_grid = np.linspace(alt[0], alt[len(alt)-1], nlevel_new)
    new_temp = np.interp(new_grid, alt, temp)
    new_pres = np.interp(new_grid, alt, pres)

    if PressUnits == 'Pa':
        new_pres = new_pres*u.bar.to(u.Pa)

    dat = Table([new_pres[::-1], new_temp[::-1]], names=('Press', 'Temp'))
    ascii.write(dat, AtmProfPath+'PT_profile_'+casename+'.pt', overwrite=True)

# MEGAN NOTE YOU SHOULD CREATE THE MIXING RATIO COLS TO BE THE SAME AS PTZ OUT (OR JUST USE THAT FILE IF THEY DONT NEED TO BE DEGRADED)


### Purpose: To compile all of the steps in a VPL Climate Run in a quick
#### and easy way to a fast-readable python dictionary using pandas dataframe.
#
## Inputs:
# casename - name of case to find output file naming scheme
# nlevel - number of atmospheric layers (use nlevel_coarse)
# trynum - the iteration number youre on for the specific case
# save_dic - False or output file name, option to save the python dict as json 
#       Note, output file name HAS TO HAVE a '.json' extension
# OutPath - path where output climate file is
# DataOutPath - path to write a dictionary data file to, maybe same as OutPath?
##
def explore_climate_out(casename, nlevel, trynum=1, save_dic=False, 
                        OutPath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/', 
                        DataOutPath='./'):
    # Define python dictionary to compile data in
    dat = {}
    dat['Atm_Levels'] = nlevel

    # Column names used in VPL Climate output run
    colnames = ['P[Pa]', 'Alt[km]', 'T[K]', 'Q_s[K/day]', 'Q_t[K/day]', 'Q_c[K/day]', 
                'Q_ad[K/day]', 'Q_net[K/day]', 'fs_net[W/m/m]', 'ft_net[W/m/m]', 'fc[W/m/m]', 
                'f_ad[W/m/m]', 'pc[Pas]', 'Altc[km]', 'Tc[K]', 'dt[s]', 'lr[K/km]', 
                'aid_lr[K/km]', 'Km[m2/s]', 'rmix[kg/kg]']

    # Open a simple text instance of the Climate output to use for parsing
    if trynum == 1:
        dat['FileName'] = OutPath+'vpl_climate_output_'+casename+'.run'
        fop = open(OutPath+'vpl_climate_output_'+casename+'.run', 'r')
    else:
        dat['FileName'] = OutPath+'vpl_climate_output_'+casename+'_Try'+str(trynum)+'.run'
        fop = open(OutPath+'vpl_climate_output_'+casename+'_Try'+str(trynum)+'.run', 'r')
    flines = fop.readlines()
    fop.close()

    IDcount = 1 # Keep track of how many atmospheric profiles there are for ID and plotting
    # Loop through text instance of output file to retrieve atmospheric profiles
    for i in range(len(flines)):
        curr_line = flines[i].split()
        if len(curr_line) > 0:
            if curr_line[0] == '(Pas)': # This checks if you're at a profile
                # This reads in that profile beautifully as pandas data frame
                curr_step = pd.read_csv(dat['FileName'], delimiter=' ', skipinitialspace=True, header=0, 
                                        names=colnames, skiprows=i, nrows=nlevel)
                
                # Add that profile to the dictionary
                dat[str(IDcount)] = {} # step is the ID count number
                for k in colnames:
                    dat[str(IDcount)][k] = np.array(curr_step[k])

                # Find net flux at each level
                Fnet = np.zeros(len(curr_step['fs_net[W/m/m]']))
                for lvl in range(len(Fnet)):
                    Fnet[lvl] = curr_step['fs_net[W/m/m]'][lvl] - curr_step['ft_net[W/m/m]'][lvl] - curr_step['fc[W/m/m]'][lvl]
                dat[str(IDcount)]['f_net[W/m/m]'] = Fnet

                IDcount += 1

    # Find the dT from step to step
    for i in range(1,IDcount):
        dtemp = (dat[str(i)]['T[K]'] - dat['1']['T[K]'])/dat['1']['T[K]']
        dat[str(i)]['dT'] = dtemp

    # Save the number of atmospheric profiles that were in the output file
    dat['Num_Steps'] = IDcount-1

    # If save_dic is not false, save the dictionary as a json with the specified output file name
    if save_dic != False:
        dichold = json.dumps(dat)
        final = open(DataOutPath+save_dic, 'w')
        json.dump(dichold, final)
        final.close()

    return dat

### Check the photochem output for convergence, might need to be played with
##
## Inputs:
# casename - name of case you're running (to find output file of photochem)
# trynum - the iteration number youre on for the specific case
# NormGrossTolerance - Convergence check, abs value of normalized gross error from photochem run
#          must be <= this tolerance to be converged
# L2Tolerance - Convergence check, abs value of L2 Norm error from photochem run must be <= this
#          tolerance to be converged
# TimeTolerance - Convergence check, time of last step must be >= this tolerance to be converged
# OutPath - the path where model run outputs have been written
##
def check_photochem_conv(casename, trynum=1, NormGrossTolerance=5, L2Tolerance=5, TimeTolerance=1e16, 
                         OutPath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/'):
    # Set the output flag of converged or not (boolean)
    # Guilty until proven innocent
    HasItConverged = False

    # Flags for each tolerance check
    NormGrossConverged = False
    L2Converged = False
    TimeConverged = False 

    # Read in the output from the photochem run, try number defines naming scheme for automatic pipeline
    if trynum == 1:
        fi = open(OutPath+'photochem_run_output_'+casename+'.run', 'r')
    else:
        fi = open(OutPath+'photochem_run_output_'+casename+'_Try'+str(trynum)+'.run', 'r')
    
    lines = fi.readlines()
    fi.close()

    # Want to get the normalized gross error, l2 error, and last outputed time
    # These are all at the end of the file, loop through backwards and break iteration to convserve efficiency
    for i in reversed(range(len(lines))):
        if len(lines[i].split()) > 0: # Make sure you dont get length errors
            # Check if you're at the lines of Normalized Gross and L2 error
            if lines[i].split()[0] == 'Normalized':
                vals = lines[i+1].split() # Get Normalized Gross and L2 error
                NormGrosserr = float(vals[0])
                L2err = float(vals[1])
            # Check if youre at the last timestep line
            elif lines[i].split()[0] == 'N':
                hold = lines[i].split()
                for k in range(len(hold)):
                    if hold[k] == 'TIME': # Find the final time from the last timestep
                        FinalTime = float(hold[k+2])
                        break
                # If you have all desired values, break the loop
                break

    # Do convergence checking:
    if NormGrosserr <= NormGrossTolerance:
        NormGrossConverged = True

    if L2err <= L2Tolerance:
        L2Converged = True

    if FinalTime >= TimeTolerance:
        TimeConverged = True

    # Overall convergence check:
    if NormGrossConverged == True and L2Converged == True and TimeConverged == True:
        HasItConverged = True

    # Print messages:
    if HasItConverged:
        print('Photochem run '+casename+' Try '+str(trynum)+' has converged!')
    else:
        print('Photochem run '+casename+' Try '+str(trynum)+' has NOT converged.')
    print('Normalized Gross error: '+str(NormGrosserr))
    print('L2 Error: '+str(L2err))
    print('Time of final timestep: '+str(FinalTime)+'\n')

    return HasItConverged, NormGrosserr, L2err, FinalTime
# usage should be 'convergence, grosserr, l2err, finaltime = check_photochem_conv()

### Check the vpl climate output for convergence, might need more work currently pretty unconstrained
##
## Inputs:
# casename - name of case you're running (to find output file of climate)
# trynum - the iteration number youre on for the specific case
# TropHeatingTolerance - Convergence check, last output value of avg trop heatin rate magnitude must be
#          <= this tolerance [K/day] to be converged
# AvgFluxTolerance - Convergence check, last output value of avg flux must be <= this tolerance
#          [W/m^2] to be converged
# OutPath - the path where model run outputs have been written
##
def check_vplclimate_conv(casename, trynum=1, TropHeatingTolerance=1e-4, AvgFluxTolerance=0.1, 
                          OutPath='/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/'):
    # Set the output flag of converged or not (boolean)
    # Guilty until proven innocent
    HasItConverged = False

    # Flags for each tolerance check
    TropHeatingConverged = False
    AvgFluxConverged = False

    # Read in the output from the climate run, try number defines naming scheme for automatic pipeline
    if trynum == 1:
        fi = open(OutPath+'vpl_climate_output_'+casename+'.run', 'r')
    else:
        fi = open(OutPath+'vpl_climate_output_'+casename+'_Try'+str(trynum)+'.run', 'r')
    
    lines = fi.readlines()
    fi.close()

    # Want to get the last output trop heating rate and avg flux, should be last two lines
    # so loop in reversed order, break loop after to conserve efficiency
    for i in reversed(range(len(lines))):
        hold = lines[i].split()
        if hold[0] == 'avg' and hold[1] == 'flux:':
            AvgFlux = float(hold[2])
        elif hold[0] == 'avg' and hold[1] == 'trop':
            TropHeating = float(hold[5])
            # After retrieving tropospheric heating rate will have all values, break
            break

    # Do Convergence checking
    if TropHeating <= TropHeatingTolerance:
        TropHeatingConverged = True

    if AvgFlux <= AvgFluxTolerance:
        AvgFluxConverged = True

    # Overall convergence check:
    if TropHeatingConverged == True and AvgFluxConverged == True:
        HasItConverged = True

    # Print messages:
    if HasItConverged:
        print('VPL Climate run '+casename+' Try '+str(trynum)+' has converged!')
    else:
        print('VPL Climate run '+casename+' Try '+str(trynum)+' has NOT converged.')
    print('Avg Tropospheric Heating Rate Magnitude: '+str(TropHeating)+' K/day')
    print('Avg Flux: '+str(AvgFlux)+' W/m**2\n')

    return HasItConverged, TropHeating, AvgFlux
# usage should be 'convergence, tropheating, avgflux = check_vplclimate_conv()
