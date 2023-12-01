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

# Molecule dict example:
# molecules_I_want = {'O2':7, 'O3':3}
# where 7 and 3 are the HITRAN gas codes for O2 and O3, respectively

class VPLModelingPipeline:

    # Set Global and initialize atmosphere object:
    def __init__(self, casename, photochemInitial, vplclimateInitial, verbose, molecule_dict, hitran_year='2020') -> None:
        # Set any and all needed paths
        self.photochemDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/' # path to PHOTOCHEM/ dir
        self.atmosDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/' # path to atmos/ dir
        self.lblabcDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/' # path to lblabc/ dir (such that lblabcDir/lblabc is the executable to call)
        self.OutPath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/' # path for the raw model run outputs (NOT for created data products like dictionaries)
        self.DataOutPath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/' # path for created data products like dictionaries
        self.AtmProfPath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/AtmProfiles/' # path to put atmospheric profile files (.pt files really)
        self.photochemBackupDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/Save_Photochem_Output/'+casename+'/' # path to save output from each photochem run
        self.LBLABC_AbsFilesDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/LinebyLine_absFiles/'+casename+'/' # path to put the created lbl .abs files in 
        self.lblabc_RunScriptDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/LBLABC/'+casename+'/' # path to put lbl runscripts in
        self.vplclimate_RunScriptDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/VPLClimate/'+casename+'/' # path to put vpl climate runscripts in
        self.photochem_InputsDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/Bodies/'+casename+'/' # The path to create new photochem inputs in

        # The climate executable:
        self.vplclimate_executable = 'something_Beta17?' # The VPL Climate executable you want to use WITH FULL PATH

        # Set the appropriate HITRAN variables
        if hitran_year == '2020':
            self.HITRAN_FundamentalFile = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/fundamental2020.dat'
            self.HITRAN_parFile = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/HITRAN_Data/HITRAN2020.par'
            self.HITRAN_parFileOptionNumber = '10'
            self.lblabc_qtxt_dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/hitranQtips2020/' # For the hitran distribution you want
        elif hitran_year == '2016':
            self.HITRAN_FundamentalFile = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/fundamntl2016.dat'
            self.HITRAN_parFile = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/HITRAN_Data/HITRAN2016.par'
            self.HITRAN_parFileOptionNumber = '9'
            self.lblabc_qtxt_dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/hitranQtips/'

        # Set the number of levels for the fine (photochem) and coarse (everything else) grids
        self.nlevel_fine = 200 # Number of atm layers for fine grids (i.e., for photochem)
        self.nlevel_coarse = 32 # Number of atm layers for coarse grids (i.e., for all models besides photochem)

        # Start counters to track how many times each model has been ran
        self.num_photochem_runs = 0
        self.num_climate_runs = 0
        self.num_lblabc_runs = 0

        # Initialize convergence logic gates
        self.photochem_global_converge = False
        self.climate_global_converge = False
        self.global_convergence = False

        # User defined inputs
        self.casename = casename # Case name youre running, user defined
        self.photochemInitial = photochemInitial # path to users initial guess for photochem inputs
        self.vplclimateInitial = vplclimateInitial # path to vpl climate template file
        # User needs to create a mostly filled in vpl climate template run script
        # Pipeline will automatically point to the right PT prof/mixing ratio files and iterate the MMW
        self.verbose = verbose # Boolean, whether or not you want print statements (False for computational efficiency)
        self.molecule_dict = molecule_dict # key-value pairs of molecules of interest (keys, str) and their hitran codes (value, int)
        # MOLECULES MUST BE ALL CAPITAL LETTERS AS THEY WILL PRINT OUT FROM PHOTOCHEM

    ######################################################### Support Functions

    ### Prepare the hyak environment with ifort, python, and the HITRAN you want to use
    ##
    ## Attribute Dependencies:
    # lblabc_qtxt_dir - the qtxt dir for the hitran distribution you want to use with lblabc
    #
    ## Fxn-Specific Inputs: NONE
    ##
    def prepare_hyak_env(self):
        subprocess.run('module load intel/oneAPI/2021.1.1', shell=True)
        subprocess.run("alias python='/gscratch/vsm/gialluca/anaconda3/bin/python'", shell=True)
        subprocess.run('export LBLABC_QTXT_DIR='+self.lblabc_qtxt_dir, shell=True)

    ### Run 1 photochem run
    ##
    ## Attributes Dependencies:
    # casename - name of case you're runnin (to name make and run output files)
    # photochemDir - Path to the PHOTOCHEM dir in atmos (so you can use specific instances)
    # atmosDir - Path to the atmos/ dir
    # OutPath - path to write model make and run outputs to
    #
    ## Fxn-specific Inputs:
    # CleanMake - do a clean make of photochem or no
    # InputCopy - the path to the input photochem files, if false assumes they are already in the PHOTOCHEM/INPUTS directory
    # trynum - the iteration number you're on for the specific case, defined by self.num_photochem_runs
    ##
    def run_photochem_1instance(self, CleanMake=True, InputCopy=False, trynum=1):

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
            subprocess.run('rm '+self.photochemDir+'INPUTFILES/input_photchem.dat', shell=True)
            subprocess.run('rm '+self.photochemDir+'INPUTFILES/parameters.inc', shell=True)
            #subprocess.run('rm '+self.photochemDir+'INPUTFILES/params.dat', shell=True)
            subprocess.run('rm '+self.photochemDir+'INPUTFILES/PLANET.dat', shell=True)
            subprocess.run('rm '+self.photochemDir+'INPUTFILES/reactions.rx', shell=True)
            subprocess.run('rm '+self.photochemDir+'INPUTFILES/species.dat', shell=True)
            subprocess.run('rm '+self.photochemDir+'in.dist', shell=True)

            # Copy new input files to the right places
            subprocess.run('cp '+InputCopy+'input_photchem.dat '+self.photochemDir+'INPUTFILES/', shell=True)
            subprocess.run('cp '+InputCopy+'parameters.inc '+self.photochemDir+'INPUTFILES/', shell=True)
            #subprocess.run('cp '+InputCopy+'params.dat '+self.photochemDir+'INPUTFILES/', shell=True)
            subprocess.run('cp '+InputCopy+'PLANET.dat '+self.photochemDir+'INPUTFILES/', shell=True)
            subprocess.run('cp '+InputCopy+'reactions.rx '+self.photochemDir+'INPUTFILES/', shell=True)
            subprocess.run('cp '+InputCopy+'species.dat '+self.photochemDir+'INPUTFILES/', shell=True)
            subprocess.run('cp '+InputCopy+'in.dist '+self.photochemDir, shell=True)

        # Clear the outputs
        subprocess.run('rm -rf '+self.photochemDir+'OUTPUT/*', shell=True)
        subprocess.run('rm -rf '+self.photochemDir+'PTZ_mixingratios_in.dist', shell=True)

        # Clean make, if requested
        if CleanMake:
            if trynum == 1:
                fmake = open(self.OutPath+'photochem_make_output_'+self.casename+'.txt', 'w')
            else:
                fmake = open(self.OutPath+'photochem_make_output_'+self.casename+'_Try'+str(trynum)+'.txt', 'w')
            workdir = os.getcwd()
            os.chdir(self.atmosDir)
            subprocess.run('make -f ./PhotoMake clean', shell=True)
            subprocess.run('make -f ./PhotoMake', shell=True, stdout=fmake)
            os.chdir(workdir)

        # Run photochem
        if trynum == 1:
            f = open(self.OutPath+'photochem_run_output_'+self.casename+'.run', 'w')
        else:
            f = open(self.OutPath+'photochem_run_output_'+self.casename+'_Try'+str(trynum)+'.run', 'w')
        workdir = os.getcwd()
        os.chdir(self.atmosDir)
        subprocess.run('./Photo.run', shell=True, stdout=f)
        os.chdir(workdir)

    ### Run VPL Climate for a given runscript/executable
    ##
    ## Attribute Dependencies:
    # casename - name of case you're running (to name output file)
    # OutPath - path to put run output in
    #
    ## Fxn-specific Inputs:
    # runscript - name of VPL Climate runscript WITH PATH
    # exec - the VPL Climate executable to use WITH PATH
    # trynum - the iteration number youre on for the specific case
    ##
    def run_climate_1instance(self, runscript, exec, trynum=1):
        if trynum == 1:
            subprocess.run(exec+' < '+runscript+' > '+self.OutPath+'vpl_climate_output_'+self.casename+'.run', shell=True)
        else:
            subprocess.run(exec+' < '+runscript+' > '+self.OutPath+'vpl_climate_output_'+self.casename+'_Try'+str(trynum)+'.run', shell=True)

    ### Run LBLABC for a given runscript (just useful to get the output w/e)
    ##
    ## Attribute Dependencies:
    # casename - name of case you're running (to name output file)
    # OutPath - path to put run output in
    #
    # Fxn-specific Inputs:
    # runscript - name of LBLABC runscript WITH PATH
    ##
    def run_lblabc_1instance(self, runscript, molecule):
        f = open(self.OutPath+'lblabc_run_output_'+molecule+'_'+self.casename+'.run', 'w')
        workdir = os.getcwd()
        os.chdir(self.lblabcDir)
        subprocess.run(self.lblabcDir+'lblabc < '+runscript, shell=True, stdout=f)
        f.close()
        os.chdir(workdir)

    ### Take the PT profile output from photochem and degrade it to a specified number of layers ...
    ### ... to create PT profile for LBLABC and SMART
    ##
    ## Attribute Dependencies:
    # casename - name of case on your grid you're doing
    # nlevel_coarse - number of layers you want in your degraded atmosphere
    # AtmProfPath - path to output new PT profile to
    #
    ## Fxn-specific Inputs:
    # PressUnits - can do Bar or Pa, Megan stick with bar for the forseeable forever
    ##
    def degrade_PT(self, PressUnits='Bar'):
        atm = ascii.read(self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist', delimiter=' ')
        alt = atm['ALT']
        pres = atm['PRESS']
        temp = atm['TEMP']

        new_grid = np.linspace(alt[0], alt[len(alt)-1], self.nlevel_coarse)
        new_temp = np.interp(new_grid, alt, temp)
        new_pres = np.interp(new_grid, alt, pres)

        if PressUnits == 'Pa':
            new_pres = new_pres*u.bar.to(u.Pa)

        dat = Table([new_pres[::-1], new_temp[::-1]], names=('Press', 'Temp'))
        ascii.write(dat, self.AtmProfPath+'PT_profile_'+self.casename+'.pt', overwrite=True)

    # MEGAN NOTE YOU SHOULD CREATE THE MIXING RATIO COLS TO BE THE SAME AS PTZ OUT (OR JUST USE THAT FILE IF THEY DONT NEED TO BE DEGRADED)


    ### Purpose: To compile all of the steps in a VPL Climate Run in a quick
    #### and easy way to a fast-readable python dictionary using pandas dataframe.
    #
    ## Attribute Dependencies:
    # casename - name of case to find output file naming scheme
    # nlevel_coarse - number of atmospheric layers (use nlevel_coarse)
    # OutPath - path where output climate file is
    # DataOutPath - path to write a dictionary data file to, maybe same as OutPath?
    # 
    ## Fxn-specific Inputs:
    # trynum - the iteration number youre on for the specific case, defined by self.num_climate_runs
    # save_dic - False or output file name, option to save the python dict as json 
    #       Note, output file name HAS TO HAVE a '.json' extension
    ##
    def explore_climate_out(self, trynum=1, save_dic=False):
        # Define python dictionary to compile data in
        dat = {}
        dat['Atm_Levels'] = self.nlevel_coarse

        # Column names used in VPL Climate output run
        colnames = ['P[Pa]', 'Alt[km]', 'T[K]', 'Q_s[K/day]', 'Q_t[K/day]', 'Q_c[K/day]', 
                    'Q_ad[K/day]', 'Q_net[K/day]', 'fs_net[W/m/m]', 'ft_net[W/m/m]', 'fc[W/m/m]', 
                    'f_ad[W/m/m]', 'pc[Pas]', 'Altc[km]', 'Tc[K]', 'dt[s]', 'lr[K/km]', 
                    'aid_lr[K/km]', 'Km[m2/s]', 'rmix[kg/kg]']

        # Open a simple text instance of the Climate output to use for parsing
        if trynum == 1:
            dat['FileName'] = self.OutPath+'vpl_climate_output_'+self.casename+'.run'
            fop = open(self.OutPath+'vpl_climate_output_'+self.casename+'.run', 'r')
        else:
            dat['FileName'] = self.OutPath+'vpl_climate_output_'+self.casename+'_Try'+str(trynum)+'.run'
            fop = open(self.OutPath+'vpl_climate_output_'+self.casename+'_Try'+str(trynum)+'.run', 'r')
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
                                            names=colnames, skiprows=i, nrows=self.nlevel_coarse)
                    
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
            final = open(self.DataOutPath+save_dic, 'w')
            json.dump(dichold, final)
            final.close()

        return dat

    ### Check the photochem output for convergence, might need to be played with
    ##
    ## Attribute Dependencies:
    # casename - name of case you're running (to find output file of photochem)
    # OutPath - the path where model run outputs have been written
    #
    ## Fxn-specific Inputs:
    # trynum - the iteration number youre on for the specific case, defined by self.num_photochem_runs
    # NormGrossTolerance - Convergence check, abs value of normalized gross error from photochem run
    #          must be <= this tolerance to be converged
    # L2Tolerance - Convergence check, abs value of L2 Norm error from photochem run must be <= this
    #          tolerance to be converged
    # TimeTolerance - Convergence check, time of last step must be >= this tolerance to be converged
    ##
    def check_photochem_conv(self, trynum=1, NormGrossTolerance=5, L2Tolerance=5, TimeTolerance=1e16):
        # Set the output flag of converged or not (boolean)
        # Guilty until proven innocent
        HasItConverged = False

        # Flags for each tolerance check
        NormGrossConverged = False
        L2Converged = False
        TimeConverged = False 

        # Read in the output from the photochem run, try number defines naming scheme for automatic pipeline
        if trynum == 1:
            fi = open(self.OutPath+'photochem_run_output_'+self.casename+'.run', 'r')
        else:
            fi = open(self.OutPath+'photochem_run_output_'+self.casename+'_Try'+str(trynum)+'.run', 'r')
        
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
        if self.verbose == True:
            if HasItConverged:
                print('Photochem run '+self.casename+' Try '+str(trynum)+' has converged!')
            else:
                print('Photochem run '+self.casename+' Try '+str(trynum)+' has NOT converged.')
            print('Normalized Gross error: '+str(NormGrosserr))
            print('L2 Error: '+str(L2err))
            print('Time of final timestep: '+str(FinalTime)+'\n')

        return HasItConverged, NormGrosserr, L2err, FinalTime
    # usage should be 'convergence, grosserr, l2err, finaltime = check_photochem_conv()

    ### Check the vpl climate output for convergence, might need more work currently pretty unconstrained
    ##
    ## Attribute Dependencies:
    # casename - name of case you're running (to find output file of climate)
    # OutPath - the path where model run outputs have been written
    #
    ## Fxn-specific Inputs:
    # trynum - the iteration number youre on for the specific case, defined by self.num_climate_runs
    # TropHeatingTolerance - Convergence check, last output value of avg trop heatin rate magnitude must be
    #          <= this tolerance [K/day] to be converged
    # AvgFluxTolerance - Convergence check, last output value of avg flux must be <= this tolerance
    #          [W/m^2] to be converged
    ##
    def check_vplclimate_conv(self, trynum=1, TropHeatingTolerance=1e-4, AvgFluxTolerance=0.1):
        # Set the output flag of converged or not (boolean)
        # Guilty until proven innocent
        HasItConverged = False

        # Flags for each tolerance check
        TropHeatingConverged = False
        AvgFluxConverged = False

        # Read in the output from the climate run, try number defines naming scheme for automatic pipeline
        if trynum == 1:
            fi = open(self.OutPath+'vpl_climate_output_'+self.casename+'.run', 'r')
        else:
            fi = open(self.OutPath+'vpl_climate_output_'+self.casename+'_Try'+str(trynum)+'.run', 'r')
        
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
        if self.verbose == True:
            if HasItConverged:
                print('VPL Climate run '+self.casename+' Try '+str(trynum)+' has converged!')
            else:
                print('VPL Climate run '+self.casename+' Try '+str(trynum)+' has NOT converged.')
            print('Avg Tropospheric Heating Rate Magnitude: '+str(TropHeating)+' K/day')
            print('Avg Flux: '+str(AvgFlux)+' W/m**2\n')

        return HasItConverged, TropHeating, AvgFlux
    # usage should be 'convergence, tropheating, avgflux = check_vplclimate_conv()

    ### Save the photochem run outputs to the backup directory
    ##
    ## Attribute Dependencies:
    # photochemBackupDir - The directory created for the casename to make backups (this function will create subdirs in that)
    #
    ## Fxn-specific Inputs:
    # trynum - the iteration number youre on for the specific case, defined by self.num_photochem_runs
    ##
    def backup_photochem_run(self, trynum=1):
        os.mkdir(self.photochemBackupDir+'/RunNumber'+str(trynum)+'/')
        subprocess.run('cp '+self.photochemDir+'OUTPUT/* '+self.photochemBackupDir+'/RunNumber'+str(trynum)+'/', shell=True)
        if self.verbose == True:
            print('Photochem Run Number '+str(trynum)+' Output Backup Created')

    ### Create the lblabc runscripts, this creates them from scratch, but this function could easily be re-written to use a
    ### ... template, as long as it writes an lbl file to the self.LBLABC_AbsFilesDir
    ##
    ## Attribute Dependencies:
    # casename - name of the case youre running
    # LBLABC_AbsFilesDir - directory path to write new lblabc .abs files in
    # lblabc_RunScriptDir - directory path to write runscripts to
    # molecule_dict - dictionary with the molecules of interest (keys), their hitran codes (value), and rmix columns
    #
    ## Fxn-specific Inputs: NONE
    ##
    def make_lblabc_runscripts(self):
        for i in self.molecule_dict['Gas_names']:
            f = open(self.lblabc_RunScriptDir+'RunLBLABC_'+i+'_'+self.casename+'.script', 'w')
            f.write('3                                       format list directed\n')
            f.write(self.AtmProfPath+'PT_profile_'+self.casename+'.pt\n')
            f.write('1                                       # records to skip at TOF\n')
            f.write('1,2                                     columns of p and t\n')
            f.write('100000                                  pressure scaling factor\n')
            f.write('1                                       number of t profiles\n')
            f.write(str(self.molecule_dict[i])+'                                       hitran gas index\n')
            f.write('0                                       # of foreign broadening gases\n')
            f.write('3                                       format  (list directed)\n')
            f.write(self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist\n')
            f.write('1                                       # records to skip at TOF\n')
            f.write('1                                       vert coordinate (pres)\n')
            f.write('1,'+str(self.molecule_dict[i+'_RmixCol'])+'                                     columns of p and rmix\n')
            f.write('1                                       type - volume mixing ratio (mass = 2)\n')
            f.write('100000,1.0                              p, rmix scaling factors\n')
            f.write(self.HITRAN_parFileOptionNumber+'                                       read HITRAN .par formatted file\n')
            f.write(self.HITRAN_parFile+'\n')
            f.write(str(self.planetary_gravity)+'                                   gravitational acceleration\n') # grav accel
            f.write(str(self.planetary_radius)+'                               planet radius\n') # planet radius
            f.write(str(self.MMW)+'                      mol. wgt. of atmosphere (kg/kmole)\n')
            f.write('50,100000                               min, max wavenumber\n')
            f.write('200.                                    maximum line width\n')
            f.write('1.e-5                                   minimum column optical depth\n')
            f.write(self.HITRAN_FundamentalFile+'\n')
            f.write(self.LBLABC_AbsFilesDir+i+'_'+self.casename+'.abs\n')
            f.write('2                                       overwrite existing files\n')
            f.close()

    # THIS IS THE FUNCTION TO BE CAREFUL OF
    # IT IS *HIGHLY* DEPENDENT ON HOW YOU WRITE YOUR CLIMATE TEMPLATE FILE
    # A NEW USER SHOULD TEST THIS WORKS BEFORE RUNNING THE AUTOMATIC PIPELINE
    ##
    ### Take user defined vpl climate run file, update profile, .abs file paths, and MMW
    ##
    ## Attribute Dependencies:
    #
    #
    ## Fxn-specific Inputs:
    ##
    def make_first_climate_runscript(self):
        # Start a new runscript to create, and open/read the template user defined climate script
        f = open(self.vplclimate_RunScriptDir+'RunVPLClimate_'+self.casename+'_FirstTry.script', 'w')
        template = open(self.vplclimateInitial, 'r')
        lines = template.readlines()
        template.close()

        # Iterate through the template runscript, only change necessary paths and MMW
        PTwritten = False # Flag if PT profile has been written yet
        RmixWritten = True # Flag to know if you are in a gas-specific block that needs a rmix profile written
        for i in range(len(lines)):
            # If its the mean mol weight line, write the MMW:
            if len(lines[i].split('atm mean mol wgt')) > 1:
                f.write(self.MMW+'			atm mean mol wgt [g\mol]\n')
            # If the PT profile has not been written yet, and the line 2 away is for the columns of P,T, write the pt profile
            elif PTwritten == False and len(lines[i+2].split('columns of P,T')) > 1:
                f.write(self.AtmProfPath+'PT_profile_'+self.casename+'.pt\n')
                PTwritten == True
            # If the line is the columns of P,T for the pt profile, they will be '1,2'
            elif len(lines[i].split('columns of P,T')) > 1:
                f.write('1,2			columns of P,T\n')
            # If youve reached a hitran gas code line, decide what the gas will be based on the index
            elif len(lines[i].split('HITRAN GAS CODE')) > 1:
                curr_gas_index = int(lines[i].split()[0])
                for g in self.molecule_dict['Gas_names']:
                    if self.molecule_dict[g] == curr_gas_index:
                        curr_gas = g
                RmixWritten = False
                f.write(lines[i])
            ## You should be in gas-specific info block by this point:
            # If the line youre on is a .abs HITRAN line absorber file, write the appropriate on for the gas youre on
            elif RmixWritten == False and len(lines[i-2].split('HITRAN Line Absorbers')) > 1:
                f.write(self.LBLABC_AbsFilesDir+curr_gas+'_'+self.casename+'.abs\n')
            # If the line youre on is for the mixing ratios, write the atmospheric profile
            elif RmixWritten == False and len(lines[i+3].split('columns of P,rmix')) > 1:
                f.write(self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist\n')
            # If the line youre on is for the columns of P,rmix, give appropriate column number
            elif RmixWritten == False and len(lines[i].split('columns of P,rmix')) > 1:
                f.write('1,'+str(self.molecule_dict[curr_gas+'_RmixCol'])+'			columns of P,rmix\n')
                RmixWritten = True

    # Where you left off: Just wrote the correct Rmix column for the gas youre on
    # Need to go through rest of runclimate runscript to make sure nothing else needs to be written explicitly
    # also whats the deal with pressure scaling? Does it need to convert to Pa from bars? I think sooooo 






    ### AUTOMATIC PIPELINE, calls all above functions sequentially defined by flow chart, updates necessary object parameters
    ##
    ## Dependent on all Attributes, no fxn-specific inputs
    ##
    def run_automatic(self):
        # Prepare your backup directory for photochem data
        if not os.path.exists(self.photochemBackupDir):
            os.mkdir(self.photochemBackupDir)
        # Prepare your directory for storing lblabc .abs files
        if not os.path.exists(self.LBLABC_AbsFilesDir):
            os.mkdir(self.LBLABC_AbsFilesDir)
        # Prepare directory for storing lblabc run script files
        if not os.path.exists(self.lblabc_RunScriptDir):
            os.mkdir(self.lblabc_RunScriptDir)
        # Prepare directory for storing vpl climate run script files
        if not os.path.exists(self.vplclimate_RunScriptDir):
            os.mkdir(self.vplclimate_RunScriptDir)
        # Prepare directory for storing new photochem inputs
        if not os.path.exists(self.photochem_InputsDir):
            os.mkdir(self.photochem_InputsDir)
        # Prepare the Hyak environment
        self.prepare_hyak_env()

        # Start loop to find global convergence with photochem + lblabc + vpl climate
        while self.global_convergence == False:


            ### Run photochem section start ---------------------

            self.num_photochem_runs += 1
            if self.verbose == True:
                print('----> Beginning photochem run Try number '+str(self.num_photochem_runs))
            if self.num_photochem_runs == 1: # If this is the first time, user defined input path, dont need to check global convergence
                self.run_photochem_1instance(CleanMake=True, InputCopy=self.photochemInitial, trynum=self.num_photochem_runs)
            else:
                self.run_photochem_1instance(CleanMake=True, InputCopy='WHEREVER NEW INPUTS ARE CREATED', trynum=self.num_photochem_runs)
            
            ### If this is the first run, Retrieve the surf gravity and radius of the planet from PLANET.dat
            if self.num_photochem_runs == 1:
                planet = open(self.photochemDir+'INPUTFILES/PLANET.dat', 'r')
                planet_lines = planet.readlines()
                planet.close()
                grav = float(planet_lines[0].split()[0])*1e-2 # Get the gravity from the first line of PLANET.dat and convert to m/s**2 (should be originally cm/s**2)
                # Get the radius
                for i in planet_lines:
                    hold = i.split()
                    if len(hold) > 4:
                        if hold[3] == 'radius':
                            rad = float(hold[0])*1e-5 # Get radius from PLANET.dat and convert from cm to km
                            break
                # Set object values
                self.planetary_gravity = grav
                self.planetary_radius = rad

            photochem_subtries = 1
            # Currently, the trynum will only refer to the converged case (should we save every single try no matter what?)
            local_photochem_conv, grosserr, l2err, finaltime = self.check_photochem_conv()

            # If photochem did not converge, try try again
            while local_photochem_conv == False:
                subprocess.run('rm -rf '+self.photochemDir+'in.dist', shell=True)
                subprocess.run('cp '+self.photochemDir+'OUTPUTS/out.dist '+self.photochemDir+'in.dist', shell=True)
                self.run_photochem_1instance(CleanMake=False, InputCopy=False, trynum=self.num_photochem_runs)
                photochem_subtries += 1
                local_photochem_conv, grosserr, l2err, finaltime = self.check_photochem_conv()

            if self.verbose == True:
                print('Photochem local convergence found with '+str(photochem_subtries)+' subtries')

            # Save backup of photochem output
            self.backup_photochem_run(trynum=self.num_photochem_runs)

            # Retrieve atmosphere MMW for use moving forward:
            if self.num_photochem_runs == 1:
                fi = open(self.OutPath+'photochem_run_output_'+self.casename+'.run', 'r')
            else:
                fi = open(self.OutPath+'photochem_run_output_'+self.casename+'_Try'+str(self.num_photochem_runs)+'.run', 'r')
            lines = fi.readlines()
            fi.close()
            for i in lines:
                if len(i.split('Molecular weight of atmosphere')) > 1:
                    self.MMW = float(i.split()[len(i.split())-1])
                    break

            # If you only had 1 subtry, and its not the first run, check for global convergence
            if photochem_subtries == 1 and self.num_photochem_runs != 1:
                print('THIS IS WHERE YOU PUT THE GLOBAL CONVERGENCE PHOTOCHEM CHECK MEG')

            ### Run photochem section end --------------------------

            ### Create degraded atmospheric profile to prepare for LBLABC and Climate -------------
            self.degrade_PT()
            if self.verbose == True:
                print('Degraded PT profile created from photochem run '+str(self.num_photochem_runs))
            ### Degraded PT Profile finished ------------

            ### Run LBLABC section start ----------------------------

            # If LBLABC hasnt run yet, Get the columns of rmix for the molecules of interest, save molecule names though
            if self.num_lblabc_runs == 0:
                mixingRs = ascii.read(self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist')
                self.molecule_dict['Gas_names'] = list(self.molecule_dict.keys())
                for gas in self.molecule_dict['Gas_names']:
                    for i in range(len(mixingRs.keys())):
                        if mixingRs.keys()[i] == gas:
                            self.molecule_dict[gas+'_RmixCol'] = i+1
                            break

            # Create the LBLABC run scripts
            # Wont change much after first run, but needs to be done everytime to update MMW
            # Not much computational expense saved by only replacing MMW after each run ...
            # ... So cleaner just to recreate everytime
            self.make_lblabc_runscripts()

            # Now Run LBLABC for all the gases of interest
            for gas in self.molecule_dict['Gas_names']:
                self.run_lblabc_1instance(self.lblabc_RunScriptDir+'RunLBLABC_'+gas+'_'+self.casename+'.script', gas)
                if self.verbose == True:
                    print('LBLABC run for '+gas+', LBLABC iteration '+str(self.num_lblabc_runs+1))
            self.num_lblabc_runs += 1

            ### Run LBLABC section end -------------------------------------

            ### Run VPL Climate section start ------------------------------




