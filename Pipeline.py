import numpy as np
import os,sys,shutil,subprocess
from astropy.io import ascii
import sys
import astropy.units as u
import astropy.constants as const
from astropy.table import Table
import json
import pandas as pd
import re
from multiprocessing import Pool
from NewPressure_HelperFunctions import get_true_number_densities, sum_mixing_ratios, new_total_Ndens, new_mixing_rats, find_tot_column_mass_dens

################################
##
## Working on a full semi-automatic VPL modeling pipeline
## Author: Megan Gialluca
##
################################

class VPLModelingPipeline:

    # Set Global and initialize atmosphere object:
    def __init__(self, casename, photochemInitial, verbose, find_molecules_of_interest=False, hitran_year='2020', planet='T1c', force_dz_adjust=False) -> None:
        # Set any and all needed paths
        self.photochemDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/megan_atmos/atmos/PHOTOCHEM/' #'/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/' # path to PHOTOCHEM/ dir
        self.atmosDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/megan_atmos/atmos/'# '/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/' # path to atmos/ dir
        self.lblabcDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/' # path to lblabc/ dir (such that lblabcDir/lblabc is the executable to call)
        self.OutPath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/'+casename+'/' # path for the raw model run outputs (NOT for created data products like dictionaries)
        self.DataOutPath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/'+casename+'/' # path for created data products like dictionaries
        self.AtmProfPath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/AtmProfiles/' # path to put atmospheric profile files (.pt files really)
        self.BackupPhotochemRuns = False # Make backups of individual photochem runs
        self.photochemBackupDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/Save_Photochem_Output/'+casename+'/' # path to save output from each photochem run
        self.LBLABC_AbsFilesDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/LinebyLine_absFiles/'+casename+'/' # path to put the created lbl .abs files in 
        self.lblabc_RunScriptDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/LBLABC/'+casename+'/' # path to put lbl runscripts in
        self.vplclimate_RunScriptDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/VPLClimate/'+casename+'/' # path to put vpl climate runscripts in
        self.photochem_InputsDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/'+casename+'/' # The path to create new photochem inputs in
        self.xsec_Path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/xsec/' # The path where cross section files can be found
        self.SMARTDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/smart/' # path to smart/ dir (such that smart/smart_spectra is the executable) ??????????? may need to fix these when klones online again 
        self.SMART_RunScriptDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/SMART/'+casename+'/' # path to put SMART runscripts in
        self.photochemInitialInput = photochemInitial

        # For T1-c, Mp = 1.308*u.Mearth.to(u.kg)
        self.c_NumberSolarZeniths = 4
        self.planet = planet

        if self.planet == 'T1b':
            self.planetary_mass = 1.374*u.Mearth.to(u.kg)
        elif self.planet == 'T1c':
            self.planetary_mass = 1.308*u.Mearth.to(u.kg)
        elif self.planet == 'T1d':
            self.planetary_mass = 0.388*u.Mearth.to(u.kg)
        elif self.planet == 'T1e':
            self.planetary_mass = 0.692*u.Mearth.to(u.kg)
        elif self.planet == 'T1f':
            self.planetary_mass = 1.039*u.Mearth.to(u.kg)
        elif self.planet == 'T1g':
            self.planetary_mass = 1.321*u.Mearth.to(u.kg)
        elif self.planet == 'T1h':
            self.planetary_mass = 0.326*u.Mearth.to(u.kg)
        elif self.planet == 'Earth':
            self.planetary_mass = 1*u.Mearth.to(u.kg)

        elif self.planet == 'GJ12b':
            self.planetary_mass = 0.75*u.Mearth.to(u.kg)
        
        # Force dz adjust, added for T1b high steam atmospheres to force a larger dzgrid
        self.force_dz_adjust = force_dz_adjust

        # The climate executable:
        self.vplclimate_executable = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimateModel/vpl_climate_exec' # The VPL Climate executable you want to use WITH FULL PATH

        self.HITRAN_year = hitran_year
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
        self.nlevel_coarse = 50 # Number of atm layers for coarse grids (i.e., for all models besides photochem)

        # Start counters to track how many times each model has been ran
        self.num_photochem_runs = 0
        self.num_climate_runs = 0
        self.num_lblabc_runs = 0
        self.num_2col_climate_runs = 0

        # Initialize convergence logic gates
        self.photochem_global_converge = False
        self.climate_global_converge = False
        self.global_convergence = False
        self.max_iterations_master = 300 # Never do anything more than 100x
        self.max_iterations_climate = 200 # Never run climate more than 15x
        self.suppress_IOerrors = False # if convergence fails, raise IO errors if False, or just break running function if True
        self.run_spectra = True # If true, finished a converged run with smart 
        self.rerun_smart_for_2col = True # If true, rerun smart for day and night side after climate 2 column mode
        self.fixsgbsl = False # If true, try to fix sgbsl error
        self.MCMC_pressure_only = False # If true, remove everything after photochem for an MCMC pressure fitting walker
        self.MultiNest_DataFit = False # If true, skip any climate runs, parallelize over smart calls and LBLABC
        self.clim2col_restarting = False # If true, restarting a run that was in the middle of running 2 column

        self.adjust_N2_amount = False # If true, adjust inert N2 to a specific pressure given by self.N2_fixed_pressure
        self.N2_fixed_pressure = 0.05 # pressure to adjust inert N2 to, in bars

        self.include_2column_climate = True

        # Option to start with a different 2 column dayside/nightside PT than what is found through 1D 
        self.dayside_starting_PT = None
        self.nightside_starting_PT = None

        # User defined inputs
        self.casename = casename # Case name youre running, user defined
        #self.vplclimateInitial = vplclimateInitial # path to vpl climate template file
        # User needs to create a mostly filled in vpl climate template run script
        # Pipeline will automatically point to the right PT prof/mixing ratio files and iterate the MMW
        self.verbose = verbose # Boolean, whether or not you want print statements (False for computational efficiency)
        # MOLECULES MUST BE ALL CAPITAL LETTERS AS THEY WILL PRINT OUT FROM PHOTOCHEM

        # If the atmospheric pressure needs to be checked and adjusted outside of photochem, set this flag to True
        self.adjust_atmospheric_pressure = True
        # when adjusting pressure, the number density of each level much change by less than this percentage (as a decimal) on each level to achieve convergence:
        self.NewPressure_Ndens_tolerance = 15 # DEPRECEATED
        self.NewPressure_Psurf_tolerance = 0.035 # new surface pressure must change by <= to this (x100 percent)

        # lookup table to connect Hitran gas codes to molecule names
        self.hitran_lookup = pd.read_csv('HitranTable.csv', index_col='Molecule')

        if find_molecules_of_interest == False:
            self.molecule_dict = {} # key-value pairs of molecules of interest (keys, str) and their hitran codes (value, int)
            gas_names = ['O2', 'H2O', 'O3', 'CO2', 'CO']#, 'CO', 'CO2', 'HNO3', 'N2O', 'NO2', 'SO2']
            self.molecule_dict['Gas_names'] = gas_names
            for m in range(len(gas_names)):
                self.molecule_dict[gas_names[m]] = self.hitran_lookup.loc[gas_names[m]]['HitranNumber']
                self.molecule_dict[gas_names[m]+'_RmixCol'] = m+2


    ######################################################### Support Functions

    def setup_intial_photochem_dir(self):

        # put photochem initial inputs into the inputs dir if they aren't there already
        if not os.path.exists(self.photochem_InputsDir):
            os.mkdir(self.photochem_InputsDir)
        if self.photochemInitialInput != self.photochem_InputsDir:
            subprocess.run('cp '+self.photochemInitialInput+'input_photchem.dat '+self.photochem_InputsDir, shell=True)
            subprocess.run('cp '+self.photochemInitialInput+'parameters.inc '+self.photochem_InputsDir, shell=True)
            #subprocess.run('cp '+self.photochemInitialInput+'params.dat '+self.photochem_InputsDir, shell=True)
            subprocess.run('cp '+self.photochemInitialInput+'PLANET.dat '+self.photochem_InputsDir, shell=True)
            subprocess.run('cp '+self.photochemInitialInput+'reactions.rx '+self.photochem_InputsDir, shell=True)
            subprocess.run('cp '+self.photochemInitialInput+'species.dat '+self.photochem_InputsDir, shell=True)
            subprocess.run('cp '+self.photochemInitialInput+'in.dist '+self.photochem_InputsDir, shell=True)

    ### Switch paths to Megan's local dev environment if off of hyak
    ##
    ## Attribute Dependencies: NONE
    #
    ## Fxn-Specific Inputs: NONE
    ##
    def switch_to_local_dev(self):
        self.photochemDir = '/home/mgialluca/Nextcloud/VPL_Modeling/Atmos_Dev/atmos/PHOTOCHEM/' # path to PHOTOCHEM/ dir
        self.atmosDir = '/home/mgialluca/Nextcloud/VPL_Modeling/Atmos_Dev/atmos/' # path to atmos/ dir
        self.lblabcDir = '/home/mgialluca/Nextcloud/VPL_Modeling/LBLABC_Dev/lblabc/' # path to lblabc/ dir (such that lblabcDir/lblabc is the executable to call)
        self.OutPath = '/home/mgialluca/Nextcloud/VPL_Modeling/RunOutputs/' # path for the raw model run outputs (NOT for created data products like dictionaries)
        self.DataOutPath = '/home/mgialluca/Nextcloud/VPL_Modeling/RunOutputs/' # path for created data products like dictionaries
        self.AtmProfPath = '/home/mgialluca/Nextcloud/VPL_Modeling/AtmProfiles/' # path to put atmospheric profile files (.pt files really)
        self.photochemBackupDir = '/home/mgialluca/Nextcloud/VPL_Modeling/RunOutputs/Photochem_Backup_dir/'+self.casename+'/' # path to save output from each photochem run
        self.LBLABC_AbsFilesDir = '/home/mgialluca/Nextcloud/VPL_Modeling/RunOutputs/LBL_Abs_Files/'+self.casename+'/' # path to put the created lbl .abs files in 
        self.lblabc_RunScriptDir = '/home/mgialluca/Nextcloud/VPL_Modeling/VPLModelingSupportScripts/RunFiles/LBLABC/'+self.casename+'/' # path to put lbl runscripts in
        self.vplclimate_RunScriptDir = '/home/mgialluca/Nextcloud/VPL_Modeling/VPLModelingSupportScripts/RunFiles/VPLClimate/'+self.casename+'/' # path to put vpl climate runscripts in
        self.photochem_InputsDir = '/home/mgialluca/Nextcloud/VPL_Modeling/VPLModelingSupportScripts/Bodies/'+self.casename+'/' # The path to create new photochem inputs in

        # The climate executable:
        self.vplclimate_executable = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimateModel/vpl_climate/vpl_climate' # The VPL Climate executable you want to use WITH FULL PATH

        # Set the appropriate HITRAN variables
        if self.HITRAN_year == '2020':
            self.HITRAN_FundamentalFile = '/home/mgialluca/Nextcloud/VPL_Modeling/LBLABC_Dev/lblabc/fundamental2020.dat'
            self.HITRAN_parFile = '/home/mgialluca/Nextcloud/VPL_Modeling/LBLABC_Dev/HITRAN/HITRAN2020.par'
            self.HITRAN_parFileOptionNumber = '10'
            self.lblabc_qtxt_dir = '/home/mgialluca/Nextcloud/VPL_Modeling/LBLABC_Dev/lblabc/hitranQtips2020/' # For the hitran distribution you want
        elif self.HITRAN_year == '2016':
            self.HITRAN_FundamentalFile = '/home/mgialluca/Nextcloud/VPL_Modeling/LBLABC_Dev/lblabc/fundamntl2016.dat'
            self.HITRAN_parFile = '/home/mgialluca/Nextcloud/VPL_Modeling/LBLABC_Dev/HITRAN/HITRAN2016.par'
            self.HITRAN_parFileOptionNumber = '9'
            self.lblabc_qtxt_dir = '/home/mgialluca/Nextcloud/VPL_Modeling/LBLABC_Dev/lblabc/hitranQtips/'

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
    def run_climate_1instance(self, runscript, exec, trynum=1, twocol=False):
        if twocol == False:
            if trynum == 1:
                subprocess.run(exec+' < '+runscript+' > '+self.OutPath+'vpl_climate_output_'+self.casename+'.run', shell=True)
            else:
                subprocess.run(exec+' < '+runscript+' > '+self.OutPath+'vpl_climate_output_'+self.casename+'_Try'+str(trynum)+'.run', shell=True)
        elif twocol == True:
            if trynum == 1:
                subprocess.run(exec+' < '+runscript+' > '+self.OutPath+'vpl_2col_climate_output_'+self.casename+'.run', shell=True)
            else:
                subprocess.run(exec+' < '+runscript+' > '+self.OutPath+'vpl_2col_climate_output_'+self.casename+'_Try'+str(trynum)+'.run', shell=True)



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

    ### Run LBLABC for a given runscript (just useful to get the output w/e)
    ##
    ## Attribute Dependencies:
    # casename - name of case you're running (to name output file)
    # OutPath - path to put run output in
    #
    # Fxn-specific Inputs:
    # runscript - name of LBLABC runscript WITH PATH
    ##
    def run_lblabc_1instance_Parallel(self, input):
        runscript, molecule, col = input
        if col == 'Avg':
            f = open(self.OutPath+'lblabc_run_output_'+molecule+'_'+self.casename+'.run', 'w')
        elif col == 'dayside':
            f = open(self.OutPath+'lblabc_run_output_d_'+molecule+'_'+self.casename+'.run', 'w')
        elif col == 'nightside':
            f = open(self.OutPath+'lblabc_run_output_n_'+molecule+'_'+self.casename+'.run', 'w')
        workdir = os.getcwd()
        os.chdir(self.lblabcDir)
        subprocess.run(self.lblabcDir+'lblabc < '+runscript, shell=True, stdout=f)
        f.close()
        os.chdir(workdir)
        self.num_lblabc_runs += 1

    ### Run SMART for a given runscript (just useful to get the output or be fully in python w/e)
    ##
    ## Inputs:
    # runscript - name of SMART runscript WITH PATH
    # casename - name of case you're running (to name output file)
    # outpath - path to put run output in
    ##
    def run_smart_1instance(self, runscript, whichcol=None):
        workdir = os.getcwd()
        os.chdir(self.SMARTDir)
        if whichcol == None:
            subprocess.run(self.SMARTDir+'smart_spectra < '+runscript+' > '+self.OutPath+'smart_run_output_'+self.casename+'.run', shell=True)
        elif whichcol == 'dayside':
            subprocess.run(self.SMARTDir+'smart_spectra < '+runscript+' > '+self.OutPath+'smart_run_output_dayside_'+self.casename+'.run', shell=True)
        elif whichcol == 'nightside':
            subprocess.run(self.SMARTDir+'smart_spectra < '+runscript+' > '+self.OutPath+'smart_run_output_nightside_'+self.casename+'.run', shell=True)
        os.chdir(workdir)


    ### Run SMART all columns in a way that can be parallelized when using multinest
    ##
    ## Inputs:
    # col - 'Avg', 'dayside', or 'nightside'
    ##
    def run_multinest_smart_parallel(self, col):

        if col == 'Avg':
            self.make_smart_runscript()

            self.run_smart_1instance(self.SMART_RunScriptDir+'RunSMART_'+self.casename+'.run')

            if self.verbose == True:
                ftestingoutput = open(self.OutPath+self.casename+'_SavingInfoOut.txt', 'a')
                ftestingoutput.write('Terminator SMART run completed\n')
                ftestingoutput.close()
            
        elif col == 'dayside':

            self.make_smart_runscript(whichcol='dayside')

            self.run_smart_1instance(self.SMART_RunScriptDir+'RunSMART_dayside_'+self.casename+'.run', whichcol='dayside')

            if self.verbose == True:
                ftestingoutput = open(self.OutPath+self.casename+'_SavingInfoOut.txt', 'a')
                ftestingoutput.write('Dayside SMART run completed\n')
                ftestingoutput.close()

        elif col == 'nightside':

            self.make_smart_runscript(whichcol='nightside')

            self.run_smart_1instance(self.SMART_RunScriptDir+'RunSMART_nightside_'+self.casename+'.run', whichcol='nightside')

            if self.verbose == True:
                ftestingoutput = open(self.OutPath+self.casename+'_SavingInfoOut.txt', 'a')
                ftestingoutput.write('Nightside SMART run completed\n')
                ftestingoutput.close()

    ### Take the PT profile output from photochem and degrade it to a specified number of layers ...
    ### ... to create PT profile for LBLABC and SMART
    ##
    ## Attribute Dependencies:
    # casename - name of case on your grid you're doing
    # nlevel_coarse - number of layers you want in your degraded atmosphere
    # AtmProfPath - path to output new PT profile to
    #
    ## Fxn-specific Inputs:
    # grid_spacing - to make new degraded grid with numpy linear spacing or log spacing (for thin atmospheres)
    # PressUnits - can do Bar or Pa, Megan stick with bar for the forseeable forever
    ##
    def degrade_PT(self, grid_spacing='linear', PressUnits='Bar'):
        atm = ascii.read(self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist', delimiter=' ')
        alt = atm['ALT']
        pres = atm['PRESS']
        temp = atm['TEMP']

        # Save surface temperature for use in climate from PTZ file IF climate has not been run yet
        #     ... if climate has run, surface temp should be updated automatically from the climate run output
        if self.num_climate_runs == 0:
            self.surface_temp = temp[0] # This is NOT THE SURFACE TEMPERATURE; use the 'surface level' line with Tc column

        if grid_spacing == 'linear':
            new_grid = np.linspace(alt[0], alt[len(alt)-1], self.nlevel_coarse) # target lower level to be 0.5 km for thin atmospheres
        elif grid_spacing == 'log':
            new_grid = np.logspace(np.log10(alt[0]), np.log10(alt[len(alt)-1]), self.nlevel_coarse) # target lower level to be 0.5 km for thin atmospheres
        new_temp = np.interp(new_grid, alt, temp)
        new_pres = np.interp(new_grid, alt, pres)

        if PressUnits == 'Pa':
            new_pres = new_pres*u.bar.to(u.Pa)

        dat = Table([new_pres[::-1], new_temp[::-1]], names=('Press', 'Temp'))
        ascii.write(dat, self.AtmProfPath+'PT_profile_'+self.casename+'.pt', overwrite=True)

    ### Purpose: Replace the temperature column in the PT profile with an updated temperature column;
    #### Used when climate doesn't converge on the first try
    ##
    ## Attribute Dependencies:
    # None
    #
    ## Fxn-specific Inputs:
    # new_temp - the new temperature column to put in the PT file, usually found with get_final_climate_output_temp_profile()
    ##
    def replace_PT_tempcol(self, new_temp, whichcolumn=None):

        # Read in the old PT profile
        if whichcolumn == None:
            oldpt = ascii.read(self.AtmProfPath+'PT_profile_'+self.casename+'.pt')

        elif whichcolumn == 'dayside' or whichcolumn == 'Dayside':
            oldpt = ascii.read(self.AtmProfPath+'PT_profile_dayside_'+self.casename+'.pt')

        elif whichcolumn == 'nightside' or whichcolumn == 'Nightside':
            oldpt = ascii.read(self.AtmProfPath+'PT_profile_nightside_'+self.casename+'.pt')

        # replace the temp
        oldpt['Temp'] = new_temp

        # Update Surface Temperature
        #self.surface_temp = new_temp[len(new_temp)-1]

        # Overwrite the PT profile
        if whichcolumn == None:
            ascii.write(oldpt, self.AtmProfPath+'PT_profile_'+self.casename+'.pt', overwrite=True)
        
        elif whichcolumn == 'dayside' or whichcolumn == 'Dayside':
            ascii.write(oldpt, self.AtmProfPath+'PT_profile_dayside_'+self.casename+'.pt', overwrite=True)
        
        elif whichcolumn == 'nightside' or whichcolumn == 'Nightside':
            ascii.write(oldpt, self.AtmProfPath+'PT_profile_nightside_'+self.casename+'.pt', overwrite=True)

    ### Purpose: Create the mixing ratio file with gases of interest to lblabc/climate/smart
    #
    ## Attribute Dependencies:
    # casename - name of case to find output file naming scheme
    # AtmProfPath - path to output new mixing ratio profile to
    #
    ## Fxn-specific Inputs:
    # PTZOut - the path to the PTZOut profile with mixing ratios from photochem
    # PressUnits - can do Bar or Pa, Megan stick with bar for the forseeable forever
    #
    def prep_rmix_file(self, PTZOut, PressUnits='Bar'):

        atm = ascii.read(PTZOut)
        datfortab = [atm['PRESS'][::-1]]
        namesfortab = ['Press']
        for i in self.molecule_dict['Gas_names']:
            datfortab.append(atm[i][::-1])
            namesfortab.append(i)

        if PressUnits == 'Pa':
            datfortab[0] = datfortab[0]*u.bar.to(u.Pa)

        dat = Table(datfortab, names=namesfortab)
        ascii.write(dat, self.AtmProfPath+'MixingRs_'+self.casename+'.dat', overwrite=True)

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
    
    ### Purpose: Get only the last profile output from a climate run
    ## 
    ## Attribute Dependencies:
    # casename - name of case to find output file naming scheme
    # nlevel_coarse - number of atmospheric layers (use nlevel_coarse)
    # OutPath - path where output climate file is
    # DataOutPath - path to write a dictionary data file to, maybe same as OutPath?
    # 
    ## Fxn-specific Inputs:
    # trynum - the iteration number youre on for the specific case, defined by self.num_climate_runs
    ##
    def get_final_climate_output_temp_profile(self, trynum=1):
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

        # Loop through text instance of output file to retrieve atmospheric profiles
        # Start at bottom of file to extract only the last profile
        for i in reversed(range(len(flines))):
            curr_line = flines[i].split()
            if len(curr_line) > 0:
                if curr_line[0] == '(Pas)': # This checks if you're at a profile
                    # This reads in that profile beautifully as pandas data frame
                    curr_step = pd.read_csv(dat['FileName'], delimiter=' ', skipinitialspace=True, header=0, 
                                            names=colnames, skiprows=i, nrows=self.nlevel_coarse)
                    
                    # Add that profile to the dictionary
                    for k in colnames:
                        dat[k] = np.array(curr_step[k])

                    # Find net flux at each level
                    Fnet = np.zeros(len(curr_step['fs_net[W/m/m]']))
                    for lvl in range(len(Fnet)):
                        Fnet[lvl] = curr_step['fs_net[W/m/m]'][lvl] - curr_step['ft_net[W/m/m]'][lvl] - curr_step['fc[W/m/m]'][lvl]
                    dat['f_net[W/m/m]'] = Fnet
                    break

        if self.force_dz_adjust == False:
            self.dzg = dat['Alt[km]'][0]*u.km.to(u.cm)/200

        return dat
    
    ### Purpose: Get only the last DAYSIDE and NIGHTSIDE profile output from a TWO COLUMN climate run
    ## 
    ## Attribute Dependencies:
    # casename - name of case to find output file naming scheme
    # nlevel_coarse - number of atmospheric layers (use nlevel_coarse)
    # OutPath - path where output climate file is
    # DataOutPath - path to write a dictionary data file to, maybe same as OutPath?
    # 
    ## Fxn-specific Inputs:
    # trynum - the iteration number youre on for the specific case, defined by self.num_climate_runs
    ##
    def get_final_2column_climate_output_temp_profile(self, trynum=1):
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
            dat['FileName'] = self.OutPath+'vpl_2col_climate_output_'+self.casename+'.run'
            fop = open(self.OutPath+'vpl_2col_climate_output_'+self.casename+'.run', 'r')
        else:
            dat['FileName'] = self.OutPath+'vpl_2col_climate_output_'+self.casename+'_Subtry'+str(trynum)+'.run'
            fop = open(self.OutPath+'vpl_2col_climate_output_'+self.casename+'_Subtry'+str(trynum)+'.run', 'r')
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
                                            names=colnames, skiprows=i, nrows=self.nlevel_coarse)
                    
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
    def check_photochem_conv(self, trynum=1, NormGrossTolerance=1, L2Tolerance=np.inf, TimeTolerance=1e17, subtries=0, prevtime=0):
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

        sgbslerror = False
        # Want to get the normalized gross error, l2 error, and last outputed time
        # These are all at the end of the file, loop through backwards and break iteration to convserve efficiency
        for i in reversed(range(len(lines))):
            if len(lines[i].split()) > 0: # Make sure you dont get length errors
                # Check if you're at the lines of Normalized Gross and L2 error
                if lines[i].split()[0] == 'Normalized':
                    vals = lines[i+1].split() # Get Normalized Gross and L2 error
                    NormGrosserr = float(vals[0])
                    L2err = float(vals[1])

                elif lines[i].split()[0] == 'ERROR' and lines[i].split()[2] == 'SGBSL-RHS': # Check if there was an SGBSL error
                    sgbslerror = True
                
                # Check if youre at the last timestep line
                elif lines[i].split()[0] == 'N':
                    hold = lines[i].split()
                    Nstep_photochemrun = int(lines[i].split('=')[1].split('EMAX')[0])
                    for k in range(len(hold)):
                        if hold[k] == 'TIME': # Find the final time from the last timestep
                            FinalTime = float(hold[k+2])
                            break
                    # If you have all desired values, break the loop
                    break

        # Do convergence checking:
        if subtries < 30:
            if NormGrosserr <= NormGrossTolerance:
                NormGrossConverged = True

        # Greater than 30 subtries, might be converged with high error
        else:
            if NormGrosserr < 100:
                if prevtime >= TimeTolerance and FinalTime >= TimeTolerance:
                    NormGrossConverged = True
                else:
                    NormGrossConverged = False
            else:
                NormGrossConverged = False

        if L2err <= L2Tolerance:
            L2Converged = True

        if FinalTime >= TimeTolerance:
            TimeConverged = True

        # Overall convergence check:
        if NormGrossConverged == True and TimeConverged == True: # and L2Converged == True:
            HasItConverged = True

        '''
        # Print messages:
        if self.verbose == True:
            if HasItConverged:
                print('Photochem run '+self.casename+' Try '+str(trynum)+' has converged!')
            else:
                print('Photochem run '+self.casename+' Try '+str(trynum)+' has NOT converged.')
            print('Normalized Gross error: '+str(NormGrosserr))
            print('L2 Error: '+str(L2err))
            print('Time of final timestep: '+str(FinalTime))
        '''

        return HasItConverged, NormGrosserr, L2err, FinalTime, Nstep_photochemrun, sgbslerror
    # usage should be 'convergence, grosserr, l2err, finaltime, nstepsphoto, sgbslerror = check_photochem_conv()

    ### Check the vpl climate output for convergence, might need more work currently pretty unconstrained
    ##
    ## Attribute Dependencies:
    # casename - name of case you're running (to find output file of climate)
    # OutPath - the path where model run outputs have been written
    #
    ## Fxn-specific Inputs:
    # trynum - the iteration number youre on for the specific case, defined by self.num_climate_runs
    # TropHeatingTolerance - Convergence check, last output value of avg trop heatin rate magnitude must be
    #          <= this tolerance [K/day] to be converged -> could maybe loosen this ()
    # AvgFluxTolerance - Convergence check, last output value of avg flux must be <= this tolerance
    #          [W/m^2] to be converged
    ##
    def check_vplclimate_conv(self, trynum=1, TropHeatingTolerance=9e-2, AvgFluxTolerance=1, subtries=0):

        # Relax flux criteria a bit after 30 subtries 
        if subtries > 30:
            AvgFluxTolerance = 3
            TropHeatingTolerance = 5

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
            elif hold[0] == 'surface:':
                self.surface_temp = float(hold[8])
                # After retrieving surface temp, will have all values, break
                break

        # Do Convergence checking
        if np.abs(TropHeating) <= TropHeatingTolerance:
            TropHeatingConverged = True

        if np.abs(AvgFlux) <= AvgFluxTolerance:
            AvgFluxConverged = True

        # Overall convergence check:
        if TropHeatingConverged == True and AvgFluxConverged == True:
            HasItConverged = True

        '''
        # Print messages:
        if self.verbose == True:
            if HasItConverged:
                print('VPL Climate run '+self.casename+' Try '+str(trynum)+' has converged!')
            else:
                print('VPL Climate run '+self.casename+' Try '+str(trynum)+' has NOT converged.')
            print('Avg Tropospheric Heating Rate Magnitude: '+str(TropHeating)+' K/day')
            print('Avg Flux: '+str(AvgFlux)+' W/m**2\n')
        '''

        return HasItConverged, TropHeating, AvgFlux
    # usage should be 'convergence, tropheating, avgflux = check_vplclimate_conv()

    ### Check the vpl TWO COLUMN climate output for convergence, might need more work currently pretty unconstrained
    ##
    ## Attribute Dependencies:
    # casename - name of case you're running (to find output file of climate)
    # OutPath - the path where model run outputs have been written
    #
    ## Fxn-specific Inputs:
    # trynum - the iteration number youre on for the specific case, defined by self.num_climate_runs
    # TropHeatingTolerance - Convergence check, last output value of avg trop heatin rate magnitude must be
    #          <= this tolerance [K/day] to be converged -> could maybe loosen this ()
    # AvgFluxTolerance - Convergence check, last output value of avg flux must be <= this tolerance
    #          [W/m^2] to be converged
    ##
    def check_2column_vplclimate_conv(self, trynum=1, TropHeatingTolerance=9e-2, AvgFluxTolerance=10, climsubtry=1):
        # Set the output flag of converged or not (boolean)
        # Guilty until proven innocent
        HasItConverged = False

        # Flags for each tolerance check
        TropHeatingConverged = False
        AvgFluxConverged = False

        # Read in the output from the climate run, try number defines naming scheme for automatic pipeline
        if trynum == 1:
            fi = open(self.OutPath+'vpl_2col_climate_output_'+self.casename+'.run', 'r')
        else:
            fi = open(self.OutPath+'vpl_2col_climate_output_'+self.casename+'_Try'+str(trynum)+'.run', 'r')
        
        lines = fi.readlines()
        fi.close()

        # Want to get the last output trop heating rate and avg flux, should be last two lines
        # so loop in reversed order, break loop after to conserve efficiency

        nightside_found = False
        for i in reversed(range(len(lines))):
            hold = lines[i].split()
            if len(hold) > 2:
                if hold[0] == 'avg' and hold[1] == 'flux:':
                    if nightside_found == False:
                        AvgFlux_nightside = float(hold[2])
                    else:
                        AvgFlux_dayside = float(hold[2])

                elif hold[0] == 'avg' and hold[1] == 'trop':
                    if nightside_found == False:
                        TropHeating_nightside = float(hold[5])
                    else:
                        TropHeating_dayside = float(hold[5])

                elif hold[0] == 'surface:':
                    if nightside_found == False:
                        self.surface_temp_nightside = float(hold[8])
                        nightside_found = True
                    else:
                        self.surface_temp_dayside = float(hold[8])
                        # After retrieving surface temp for nightside, will have all values, break
                        break

        # Do Convergence checking
        if climsubtry < 50:
            cnvtype = 'Tier1'
            if np.abs(TropHeating_nightside) <= TropHeatingTolerance and np.abs(TropHeating_dayside) <= TropHeatingTolerance:
                TropHeatingConverged = True

            AvgFlux = AvgFlux_dayside + AvgFlux_nightside
            if np.abs(AvgFlux) <= AvgFluxTolerance:
                AvgFluxConverged = True

            # Overall convergence check:
            if TropHeatingConverged == True and AvgFluxConverged == True:
                HasItConverged = True
        
        else:
            AvgFlux = AvgFlux_dayside + AvgFlux_nightside
            if np.abs(AvgFlux) <= AvgFluxTolerance:
                AvgFluxConverged = True

            # Overall convergence check:
            if AvgFluxConverged == True:
                HasItConverged = True

            if np.abs(TropHeating_nightside) <= TropHeatingTolerance and np.abs(TropHeating_dayside) <= TropHeatingTolerance:
                TropHeatingConverged = True
                cnvtype = 'Tier1'
            
            elif np.abs(TropHeating_nightside) <= TropHeatingTolerance or np.abs(TropHeating_dayside) <= TropHeatingTolerance:
                cnvtype = 'Tier2'
            
            else:
                cnvtype = 'Tier3'
            
        '''
        # Print messages:
        if self.verbose == True:
            if HasItConverged:
                print('VPL Climate run '+self.casename+' Try '+str(trynum)+' has converged!')
            else:
                print('VPL Climate run '+self.casename+' Try '+str(trynum)+' has NOT converged.')
            print('Avg Dayside Tropospheric Heating Rate Magnitude: '+str(TropHeating_dayside)+' K/day')
            print('Avg Nightside Tropospheric Heating Rate Magnitude: '+str(TropHeating_nightside)+' K/day')
            print('Avg Flux: '+str(AvgFlux)+' W/m**2\n')
        '''

        return HasItConverged, TropHeating_dayside, TropHeating_nightside, AvgFlux, cnvtype
    # usage should be 'convergence, tropheating_dayside, tropheating_nightside, avgflux = check_vplclimate_conv()

    ### Save the photochem run outputs to the backup directory
    ##
    ## Attribute Dependencies:
    # photochemBackupDir - The directory created for the casename to make backups (this function will create subdirs in that)
    #
    ## Fxn-specific Inputs:
    # trynum - the iteration number youre on for the specific case, defined by self.num_photochem_runs
    ##
    def backup_photochem_run(self, trynum=1):
        os.mkdir(self.photochemBackupDir+'RunNumber'+str(trynum)+'/')
        subprocess.run('cp '+self.photochemDir+'PTZ_mixingratios_in.dist '+self.photochemBackupDir+'RunNumber'+str(trynum)+'/', shell=True)
        subprocess.run('cp '+self.photochemDir+'OUTPUT/* '+self.photochemBackupDir+'RunNumber'+str(trynum)+'/', shell=True)
        #if self.verbose == True:
        #    print('Photochem Run Number '+str(trynum)+' Output Backup Created')

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
    def make_lblabc_runscripts(self, whichcol=None):
        for i in self.molecule_dict['Gas_names']:
            if whichcol == None:
                f = open(self.lblabc_RunScriptDir+'RunLBLABC_'+i+'_'+self.casename+'.script', 'w')
            elif whichcol == 'dayside':
                f = open(self.lblabc_RunScriptDir+'RunLBLABC_d_'+i+'_'+self.casename+'.script', 'w')
            elif whichcol == 'nightside':
                f = open(self.lblabc_RunScriptDir+'RunLBLABC_n_'+i+'_'+self.casename+'.script', 'w')
            f.write('3                                       format list directed\n')
            if whichcol == None:
                f.write(self.AtmProfPath+'PT_profile_'+self.casename+'.pt\n')
            elif whichcol == 'dayside':
                f.write(self.AtmProfPath+'PT_profile_dayside_'+self.casename+'.pt\n')
            elif whichcol == 'nightside':
                f.write(self.AtmProfPath+'PT_profile_nightside_'+self.casename+'.pt\n')
            f.write('1                                       # records to skip at TOF\n')
            f.write('1,2                                     columns of p and t\n')
            f.write('100000                                  pressure scaling factor\n')
            f.write('1                                       number of t profiles\n')
            f.write(str(self.molecule_dict[i])+'                                       hitran gas index\n')
            f.write('0                                       # of foreign broadening gases\n')
            f.write('3                                       format  (list directed)\n')
            f.write(self.AtmProfPath+'MixingRs_'+self.casename+'.dat\n')
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
            if self.planet == 'Earth':
                f.write('50.,100000.                               min, max wavenumber\n')
            #elif self.planet == 'GJ12b':
            #    f.write('1000.,20000.                              min, max wavenumber\n')   
            else:
                f.write('330.,20000.                               min, max wavenumber\n')
            f.write('200.                                    maximum line width\n')
            f.write('1.e-5                                   minimum column optical depth\n')
            f.write(self.HITRAN_FundamentalFile+'\n')
            if whichcol == None:
                f.write(self.LBLABC_AbsFilesDir+i+'_'+self.casename+'.abs\n')
            elif whichcol == 'dayside':
                f.write(self.LBLABC_AbsFilesDir+i+'_d_'+self.casename+'.abs\n')
            elif whichcol == 'nightside':
                f.write(self.LBLABC_AbsFilesDir+i+'_n_'+self.casename+'.abs\n')
            f.write('2                                       overwrite existing files\n')
            f.close()

    ### Purpose: For modifying settings in the climate runscript, you should edit this function directly
    #### For all settings, the prefix 'c_' indicates these are used specifically for climate
    #
    ## Attribute Dependencies:
    # None
    #
    ## Fxn-specific Inputs:
    # None
    #
    def set_climate_settings(self):

        #self.planet = 'T1c' # Would suggest keeping all settings available for each target, this provides a quick way to switch between them

        self.c_NumberTimesteps = 10000
        self.c_TimestepLength = 86400.0 # [s]
        self.c_SubstepIntervals = 20
        self.c_TimestepOutputIntervals = 50
        self.c_TimeStepMethodIndex = 7
        self.c_TempChangeTolerance = 0.01
        self.c_DoubledRadiationGrid = True
        self.c_HRTCalcType = 3 # 3 - Global Hrt
        #self.c_NumberSolarZeniths = 4 # number of solar zenith angles used in avg
        self.c_IncludeRadiativeHrt = True
        self.c_IncludeConvectiveHrt = True
        self.c_IncludeConductiveHrt = False
        # Day length in specific planet settings
        self.c_YearLength = 1.0 # Length of year in days according to c_DayLength (1 when tidally locked)
        self.c_OrbitalCalcType = 0 # 0 - fixed orbital distance
        # Semi Major Axis in specific planet settings

        #### THIS IS SPECIFIC TO THE TYPE OF ATMOSPHERE BEING TESTED, REVISIT
        self.c_NumberMajorGases = 1 # Number of major absorbing gases?
        self.c_MajorAbsorbingGas = 4 # For O2
        ####################################

        self.c_PressureJacobians = 0 # 0 - None, 1 - Radiance, 2 - Flux
        self.c_TempJacobians = 2 # 0 - None, 1 - Radiance, 2 - Flux
        self.c_Fractional_dtemp = 0.1
        self.c_SolarTolerance = 1.0 
        self.c_ThermalTolerance = 0.03
        self.c_InternalSurfaceFlux = 0.0 # [W/m2]
        self.c_ConvectiveType = 2 # 1 - adjustment, 2 - mixing length scheme, 3 - turbulent, 4 - moist mixing
        self.c_MixingLengthType = 3 # 1 - fixed, 2 - proport to scale height, 3 - Blackadar aymptotic ML
        self.c_MixingLengthProportionality = 0.01
        self.c_MinEddyDiffusivity = 0.5 # [m2/s]
        self.c_SurfaceWindSpeed = 10.0 #[m/s]
        self.c_SurfRoughnessHeight = 0.0002 # [m]
        self.c_NumberCondensibles = 0 
        self.c_NumberAerosols = 0

        # Surface Specs
        self.c_SurfaceType = 0 # 0 - Lambertian Surface
        self.c_AlbedoJacobians = 0 # 0 - None
        self.c_SurfaceProfile = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/albedo/basalt_fresh.alb'
        self.c_SurfProfile_SkipLines = 16
        self.c_SurfProfile_wlalbedo_cols = '1,2'
        self.c_SurfProfile_wlType = 1 # 1 - Wavelength
        self.c_Convert_SurfProfilewl_microns = 1.0 # Conversion factor to microns
        self.c_ScaleAlbedo = 1.0 # Factor to scale albedo

        # Heating Specs
        self.c_NumberStreams = 4
        self.c_HRTSources = 3 # 1 - Solar, 2 - Thermal, 3 - Both

        # Host Star Specs
        if self.planet == 'GJ12b':
            self.c_StellarSpectrum = '/gscratch/vsm/gialluca/StellarSpectra/gj12.dat'
            self.c_StellarSpect_SkipLines = 5
        else:
            self.c_StellarSpectrum = '/gscratch/vsm/gialluca/StellarSpectra/TRAPPIST-1_2020.dat'
            self.c_StellarSpect_SkipLines = 10
        self.c_SolarFluxUnits = 2 
        self.c_SolarSpectralUnits = 1
        self.c_Convert_Stellar_microns = 1.0 # Conversion factor to microns
        self.c_StellarSpect_wnflux_col = '1,2'

        # Output specs
        self.c_OutputsType = 1 # 1 - Fluxes
        self.c_NumberOutputAzimuths = 1
        self.c_Azimuth = 0.0
        self.c_OutputsLevel = 1 # 1 - TOA --- UNSURE ABOUT THIS AND THE LAST ONE DIFFERENCE
        self.c_OutputsUnits = 2 # 2 - Radiance
        self.c_ThermalMinMaxWn = '40.0,5000.0'
        self.c_SolarMinMaxWn = '300.0, 50000.0'
        self.c_SlitFunction = 2 # 2 - Triangular Slit function
        self.c_SolarHWHM = 10.0
        self.c_SolarRes = 10.0
        self.c_ThermalHWHM = 1.0
        self.c_ThermalRes = 1.0
        self.c_TauError = 0.5
        self.c_SingleScatterError = 0.35
        self.c_AsymmetryParamError = 0.35
        self.c_AlbedoError = 0.25
        self.c_OutputFileType = 3 # 1 - ASCII, 3 - Binary No Header Output

        if self.planet == 'T1b':
            self.c_DayLength = 130464.0 # [s], 1.510 day
            self.c_SemiMajorAxis = 0.01154 # [AU]

        elif self.planet == 'T1c':
            self.c_DayLength = 209174.4 # [s], 2.421 days
            self.c_SemiMajorAxis = 0.0158 # [AU]

        elif self.planet == 'T1d':
            self.c_DayLength = 349833.6 # [s], 4.049 days
            self.c_SemiMajorAxis = 0.02227 # [AU]
        
        elif self.planet == 'T1e':
            self.c_DayLength = 527126.4 # [s], 6.101 days
            self.c_SemiMajorAxis = 0.02925 # [AU]

        elif self.planet == 'T1f':
            self.c_DayLength = 795484.8 # [s], 9.207 days
            self.c_SemiMajorAxis = 0.03849 # [AU]

        elif self.planet == 'T1g':
            self.c_DayLength = 1067212.8 # [s], 12.352 days
            self.c_SemiMajorAxis = 0.04683 # [AU]

        elif self.planet == 'T1h':
            self.c_DayLength = 1621900.8 # [s], 18.772 days
            self.c_SemiMajorAxis = 0.06189 # [AU]
        
        elif self.planet == 'Earth':
            self.c_DayLength = 86400.0 # [s], 1 day
            self.c_YearLength = 365.25 # Length of year in days according to c_DayLength (1 when tidally locked)
            self.c_SemiMajorAxis = 1 # [AU]

            self.c_SurfaceType = 0 # 0 - Lambertian Surface
            self.c_AlbedoJacobians = 0 # 0 - None
            self.c_SurfaceProfile = '/gscratch/vsm/alinc/fixed_input/albedo/earth1.alb'
            self.c_SurfProfile_SkipLines = 6
            self.c_SurfProfile_wlalbedo_cols = '1,2'
            self.c_SurfProfile_wlType = 1 # 1 - Wavelength
            self.c_Convert_SurfProfilewl_microns = 1.0 # Conversion factor to microns
            self.c_ScaleAlbedo = 1.15 # Factor to scale albedo

            self.c_StellarSpectrum = '/gscratch/vsm/alinc/fixed_input/specs/Kurucz1cm-1_susim_atlas2_1361.dat'
            self.c_StellarSpect_SkipLines = 204
            self.c_SolarFluxUnits = 1 
            self.c_SolarSpectralUnits = 2
            self.c_Convert_Stellar_microns = 1.0 # Conversion factor to microns
            self.c_StellarSpect_wnflux_col = '1,2'

        elif self.planet == 'GJ12b':
            self.c_DayLength = 1102586.429 # [s], 12.76 days 
            self.c_SemiMajorAxis = 0.06681 # [AU]



        elif self.planet == 'Earth_SO2':
            self.c_NumberTimesteps = 10000
            self.c_TimestepLength = 86400.0 # [s]
            self.c_SubstepIntervals = 10
            self.c_TimestepOutputIntervals = 25
            self.c_TimeStepMethodIndex = 7
            self.c_TempChangeTolerance = 0.01
            self.c_DoubledRadiationGrid = True
            self.c_HRTCalcType = 3 # 3 - Global Hrt
            self.c_NumberSolarZeniths = 4 # number of solar zenith angles used in avg
            self.c_IncludeRadiativeHrt = True
            self.c_IncludeConvectiveHrt = True
            self.c_IncludeConductiveHrt = False
            self.c_DayLength = 86400.0 # [s], 2.42 days
            self.c_YearLength = 365.25 # Length of year in days according to c_DayLength (1 when tidally locked)
            self.c_OrbitalCalcType = 0 # 0 - fixed orbital distance
            self.c_SemiMajorAxis = 1.0 # [AU]

            #### THIS IS SPECIFIC TO THE TYPE OF ATMOSPHERE BEING TESTED, REVISIT
            self.c_NumberMajorGases = 1 # Number of major absorbing gases?
            self.c_MajorAbsorbingGas = 2 # For CO2
            ####################################

            self.c_PressureJacobians = 0 # 0 - None, 1 - Radiance, 2 - Flux
            self.c_TempJacobians = 2 # 0 - None, 1 - Radiance, 2 - Flux
            self.c_Fractional_dtemp = 0.1
            self.c_SolarTolerance = 1.0 
            self.c_ThermalTolerance = 0.05
            self.c_InternalSurfaceFlux = 0.0 # [W/m2]
            self.c_ConvectiveType = 2 # 1 - adjustment, 2 - mixing length scheme, 3 - turbulent, 4 - moist mixing
            self.c_MixingLengthType = 3 # 1 - fixed, 2 - proport to scale height, 3 - Blackadar aymptotic ML
            self.c_MixingLengthProportionality = 0.085
            self.c_MinEddyDiffusivity = 0.5 # [m2/s]
            self.c_SurfaceWindSpeed = 7.0 #[m/s]
            self.c_SurfRoughnessHeight = 0.004 # [m]
            self.c_NumberCondensibles = 0 
            self.c_NumberAerosols = 0

            # Surface Specs
            self.c_SurfaceType = 0 # 0 - Lambertian Surface
            self.c_AlbedoJacobians = 0 # 0 - None
            self.c_SurfaceProfile = '/gscratch/vsm/alinc/fixed_input/albedo/earth1.alb'
            self.c_SurfProfile_SkipLines = 6
            self.c_SurfProfile_wlalbedo_cols = '1,2'
            self.c_SurfProfile_wlType = 1 # 1 - Wavelength
            self.c_Convert_SurfProfilewl_microns = 1.0 # Conversion factor to microns
            self.c_ScaleAlbedo = 1.15 # Factor to scale albedo

            # Heating Specs
            self.c_NumberStreams = 4
            self.c_HRTSources = 3 # 1 - Solar, 2 - Thermal, 3 - Both

            # Host Star Specs
            self.c_StellarSpectrum = '/gscratch/vsm/alinc/fixed_input/specs/Kurucz1cm-1_susim_atlas2_1361.dat'
            self.c_StellarSpect_SkipLines = 204
            self.c_SolarFluxUnits = 1 
            self.c_SolarSpectralUnits = 2
            self.c_Convert_Stellar_microns = 1.0 # Conversion factor to microns
            self.c_StellarSpect_wnflux_col = '1,2'

            # Output specs
            self.c_OutputsType = 1 # 1 - Fluxes
            self.c_NumberOutputAzimuths = 1
            self.c_Azimuth = 0.0
            self.c_OutputsLevel = 1 # 1 - TOA --- UNSURE ABOUT THIS AND THE LAST ONE DIFFERENCE
            self.c_OutputsUnits = 2 # 2 - Radiance
            self.c_ThermalMinMaxWn = '40.0,5000.0'
            self.c_SolarMinMaxWn = '300.0, 80000.0'
            self.c_SlitFunction = 2 # 2 - Triangular Slit function
            self.c_SolarHWHM = 10.0
            self.c_SolarRes = 10.0
            self.c_ThermalHWHM = 1.0
            self.c_ThermalRes = 1.0
            self.c_TauError = 0.5
            self.c_SingleScatterError = 0.35
            self.c_AsymmetryParamError = 0.35
            self.c_AlbedoError = 0.25
            self.c_OutputFileType = 3 # 1 - ASCII, 3 - Binary No Header Output

    ### Purpose: For modifying settings in the SMART runscript, you should edit this function directly
    #### For all settings, the prefix 's_' indicates these are used specifically for SMART
    #
    ## Attribute Dependencies:
    # None
    #
    ## Fxn-specific Inputs:
    # None
    #
    def set_smart_settings(self):

        #self.planet = 'T1c' # Would suggest keeping all settings available for each target, this provides a quick way to switch between them

        self.s_PressJacobians = 0 # 0 - No pressure jacobians
        self.s_TempJacobians = 0 # 0 - No temperature jacobians

        # Following should be same as what's passed to climate
        self.s_NumberAerosols = self.c_NumberAerosols # Number of Aerosols, should be same as climate, might be redundant
        self.s_SurfaceType = self.c_SurfaceType # e.g., 0 for Lambertian reflection
        self.s_AlbedoJacobians = self.c_AlbedoJacobians # albedo jacobians
        self.s_SurfaceProfile = self.c_SurfaceProfile # Should be same as climate
        self.s_SurfProfile_SkipLines = self.c_SurfProfile_SkipLines
        self.s_SurfProfile_wlalbedo_cols = self.c_SurfProfile_wlalbedo_cols
        self.s_SurfProfile_wlType = self.c_SurfProfile_wlType
        self.s_Convert_SurfProfilewl_microns = self.c_Convert_SurfProfilewl_microns
        self.s_ScaleAlbedo = self.c_ScaleAlbedo
        self.s_SemiMajorAxis = self.c_SemiMajorAxis # account for planet specific
        if self.planet == 'Earth':
            self.s_StellarRadius = 1
        elif self.planet == 'GJ12b':
            self.s_StellarRadius = 0.2617 # [Rsun]
        else:
            self.s_StellarRadius = 0.1192 # [Rsun]

        # Rayleigh scattering
        self.s_NumRayleigh = 1 # Number of Rayleigh scatterers
        self.s_RayleighIndex = 1 # Rayleigh scatter index
        self.s_vmrRayleigh = 1.0 # Volumn Mixing Ratio of Rayleigh scatterer (?)
        
        # Heating Sources (NOTE Climate also uses these, but may be set differently)
        self.s_NumberStreams = 4 # Number of streams
        self.s_HRTSources = 3 # 1 - Solar, 2 - Thermal, 3 - Both
        
        # Host Star Specs (NOTE Climate also uses these, should match)

        self.s_StellarSpectrum = self.c_StellarSpectrum 
        self.s_StellarSpect_SkipLines = self.c_StellarSpect_SkipLines
        self.s_SolarFluxUnits = self.c_SolarFluxUnits
        self.s_SolarSpectralUnits = self.c_SolarSpectralUnits
        self.s_Convert_Stellar_microns = self.c_Convert_Stellar_microns # Conversion factor to microns
        self.s_StellarSpect_wnflux_col = self.c_StellarSpect_wnflux_col

        # Orientation Specs
        self.s_NumZenith = 1 # number of solar zenith angles
        self.s_ZenithAzimuth_Angles = '60,0' # Zenith and Azimuth angles
        self.s_ConvergenceCriteria = 0.01 # convergence criteria, unsure exactly what it is in reference to
        self.s_OutputFormat = 9 # Can't remember exact options, could see in smart manual run
        self.s_TranistType = 2 # 2 - Ray Tracing Transit, need to verify this one for model generalizability
        self.s_ImpactParam = 0.0
        self.s_LimbDarkening = 0 # 0 - No limb darkening
        self.s_NumAzimuths = 1 # Number of Azimuth angles
        self.s_AzimuthAngles = 0.0 # Azimuth angles
        self.s_OutputLevels = 1 # 1 - Top of atmosphere only
        self.s_OutputUnits = 2 # 2 - [W/m**2/sr/um]
        if self.planet == 'Earth':
            self.s_MinMax_wavenumber = '50.,100000.'
        #elif self.planet == 'GJ12b':
        #    self.s_MinMax_wavenumber = '1000.,20000.' 
        else:
            self.s_MinMax_wavenumber = '330.,20000.' #'50.,100000.'
        self.s_GridType = 2 # 2 - slit
        self.s_SpectralResponseFxn = 2 # 2 - Triangular Spectral Response Function
        self.s_FWHM = 1.0
        self.s_SampleRes = 1.0 # Sample resolution [cm**-1]
        self.s_TauBinError = 0.25 # Error for Tau Binning
        self.s_pi0BinError = 0.15 # Error for pi0 Binning
        self.s_gError = 0.15 # g Error
        self.s_AlbedoError = 0.02 # Albedo Error
        self.s_OutputFileFormat = 1 # 1 - ascii

        #if self.planet == 'T1c':





    ### Purpose: Make the vpl climate run file based on the climate settings
    #
    ## Attribute Dependencies:
    # molecule_dict - All the molecules of interest
    # casename - for naming convention
    # MMW - the atmospheric MMW found by photochem
    # vplclimate_RunScriptDir - the output location of climate runscripts
    # AtmProfPath - Path to the PT and Mixing ratio profiles
    # planetary_gravity - gravity [in m/s]
    # planetary_radius - planet radius [in km]
    # surface_temp - Surface temperature [K] gets set when calling degrade_PT()
    #
    ## Fxn-specific Inputs:
    # trynum - the iteration number youre on for the specific case, defined by self.num_climate_runs
    ##
    def make_climate_runscript(self, trynum=1):

        # Load default climate settings
        self.set_climate_settings()

        # Start a new runscript to create
        f = open(self.vplclimate_RunScriptDir+'RunVPLClimate_'+self.casename+'.script', 'w')

        f.write(str(self.c_NumberTimesteps)+'			Number of timesteps\n')
        f.write(str(self.c_TimestepLength)+'			Timestep length [s]\n')
        f.write(str(self.c_SubstepIntervals)+'			Substep intervals\n')
        f.write(str(self.c_TimestepOutputIntervals)+'			Timestep output intervals\n')
        f.write(str(self.c_TimeStepMethodIndex)+'			Index of time-stepping method\n')
        f.write(str(self.c_TempChangeTolerance)+'			Temperature change tolerance\n')
        f.write(str(self.c_DoubledRadiationGrid)+'			Doubled radiation grid\n')
        f.write(str(self.c_HRTCalcType)+'			HRT calc type index [global hrt]\n') # Andrew runs as latitude, then 4 separate versions (at each zenith angle)
        f.write(str(self.c_NumberSolarZeniths)+'			No. SZAs [solar zenith angles] used in avg\n') # maybe just start with 1 when doing sweeps (just the 60 deg)
        f.write(str(self.c_IncludeRadiativeHrt)+'			Include radiative hrt\n')
        f.write(str(self.c_IncludeConvectiveHrt)+'			Include convective hrt\n')
        f.write(str(self.c_IncludeConductiveHrt)+'			Include conductive heating rates\n')
        f.write(str(self.c_DayLength)+'			Length of day [s]\n')
        f.write(str(self.c_YearLength)+'			Length of year [days]\n')
        f.write(str(self.c_OrbitalCalcType)+'			Orbital calc type [fixed orbital distance]\n')
        f.write(str(self.c_SemiMajorAxis)+'			Distance from Star [AU]\n')
        f.write(str(self.planetary_gravity)+'			Surface gravity [m/s^2]\n')
        f.write(str(self.planetary_radius)+'			Planet radius [km]\n')
        f.write(str(self.MMW)+'			atm mean mol wgt [g/mol]\n')
        f.write(str(self.c_NumberMajorGases)+'			No. major atm gases\n')
        f.write(str(self.c_MajorAbsorbingGas)+'			***GAS ABSORBER INDEX*** [o2]\n')
        f.write('1.0			Mixing ratio\n') ### This is currently hard coded
        f.write(str(self.c_PressureJacobians)+'			Pressure jacobians? [0 = None; 1 = Radiance; 2 = Flux]\n')
        f.write(str(self.c_TempJacobians)+'			Temperature jacobians? [0 = None; 1 = Radiance; 2 = Flux]\n')
        f.write(str(self.c_Fractional_dtemp)+'			Fractional dtemp\n')
        f.write(str(self.c_SolarTolerance)+'			Solar tolerance [0-1]\n')
        f.write(str(self.c_ThermalTolerance)+'			Thermal tolerance [0-1]\n')
        f.write('3			List directed input\n')
        f.write(self.AtmProfPath+'PT_profile_'+self.casename+'.pt\n')
        f.write('1			Lines to skip\n')
        f.write('1,2			columns of P,T\n')
        f.write('100000.			Scale to Pa\n')
        f.write(str(self.surface_temp)+'			Surface temperature [K]\n') 
        f.write(str(self.c_InternalSurfaceFlux)+'			Internal surface flux [W/m2]\n')
        f.write(str(self.c_ConvectiveType)+'			Convective type [1 = adjustment; 2 = Mixing length scheme; 3 = turbulent; 4 = moist mixing]\n')
        f.write(str(self.c_MixingLengthType)+'			Mixing length type [1 = fixed; 2 = prop to scale height; 3 = Blackadar aymptotic ML]\n')
        f.write(str(self.c_MixingLengthProportionality)+'			Mixing length proportionality\n')
        f.write(str(self.c_MinEddyDiffusivity)+'			Minimum eddy diffusivity [m2/s]\n')
        f.write(str(self.c_SurfaceWindSpeed)+'			Surface Wind Speed [m/s]\n')
        f.write(str(self.c_SurfRoughnessHeight)+'			Surface roughness height [m]\n')
        f.write(str(self.c_NumberCondensibles)+'			No. condensibles\n')
        f.write('1			Generate fluxes using SMART\n')
        f.write(str(len(self.molecule_dict['Gas_names']))+'			******NO. ABSORBING GASES******\n')
        
        # Gases Get written here --------------------

        for m in self.molecule_dict['Gas_names']:
            f.write(str(self.molecule_dict[m])+'			************HITRAN GAS CODE ['+m+']*************\n')
            f.write('0			Mixing ratio jacobians [0 = None; 1 = Radiance; 2 = Flux]\n')

            no_abs_coef = 1
            use_xsec = False
            if m in ['O2', 'H2']: # Then CIA is also available
                no_abs_coef += 1
            if os.path.exists(self.xsec_Path+m.lower()+'xsec.dat'): # If an xsec file is available use that
                use_xsec = True
                no_abs_coef += 1

            f.write(str(no_abs_coef)+'			No. of abs coeff types\n')
            f.write('1			HITRAN Line Absorbers\n') # Always start with the line by line
            f.write('2			Dont compute absorption coefficients\n')
            f.write(self.LBLABC_AbsFilesDir+m+'_'+self.casename+'.abs\n')
            if m in ['O2', 'H2']: # Add CIA
                f.write('2			Collisionally-induced absorption [CIA] files\n')
                f.write(self.xsec_Path+m.lower()+'-'+m.lower()+'_abs.cia\n')
            if use_xsec == True: # Use cross sections
                f.write('3			UV-VIS absorption cross sections\n')
                f.write(self.xsec_Path+m.lower()+'xsec.dat\n')
            
            # Now pass the mixing ratios
            f.write('3			List directed input\n')
            f.write(self.AtmProfPath+'MixingRs_'+self.casename+'.dat\n')
            f.write('1 			Lines to skip\n')
            f.write('1			Use pressure coordinate\n')
            f.write('1,'+str(self.molecule_dict[m+'_RmixCol'])+'			columns of P,rmix\n')
            f.write('1			Mixing type [1 = volume; 2 = mass]\n')
            f.write('100000.,1.0			scale P and rmix\n')

        # End Gases block -------------------------

        f.write(str(self.c_NumberAerosols)+'			******NO. AEROSOLS******\n')
        f.write(str(self.c_SurfaceType)+'			Lambertian surface\n')
        f.write(str(self.c_AlbedoJacobians)+'			No albedo jacobians\n')
        f.write('3			List directed input\n')
        f.write(self.c_SurfaceProfile+'\n')
        f.write(str(self.c_SurfProfile_SkipLines)+'			Lines to skip\n')
        f.write(str(self.c_SurfProfile_wlalbedo_cols)+'			columns of wl, albedo\n')
        f.write(str(self.c_SurfProfile_wlType)+'			Wavelength\n')
        f.write(str(self.c_Convert_SurfProfilewl_microns)+'			Convert to microns\n')
        f.write(str(self.c_ScaleAlbedo)+'			Scale albedo\n')
        f.write(str(self.c_NumberStreams)+'			No. of Streams\n')
        f.write(str(self.c_HRTSources)+'			Sources: 1) Solar, 2) thermal, 3) both\n')
        f.write('3			List directed input\n')
        f.write(self.c_StellarSpectrum+'\n')
        f.write(str(self.c_StellarSpect_SkipLines)+'			Lines to skip\n')
        f.write(str(self.c_SolarFluxUnits)+'			Solar flux units\n')
        f.write(str(self.c_SolarSpectralUnits)+'			Solar spectral units\n')
        f.write(str(self.c_Convert_Stellar_microns)+'			Conversion factor to microns\n')
        f.write(str(self.c_StellarSpect_wnflux_col)+'			Cols of wn and flux\n')
        f.write(str(self.c_OutputsType)+'			Outputs [fluxes]\n')
        f.write(str(self.c_NumberOutputAzimuths)+'			No. of output azimuth angles\n')
        f.write(str(self.c_Azimuth)+'			Azimuth\n')
        f.write(str(self.c_OutputsLevel)+'			Outputs [TOA]\n')
        f.write(str(self.c_OutputsUnits)+'			Output units [radiance]\n')
        f.write(str(self.c_ThermalMinMaxWn)+'			Thermal min, max wn [cm^-1]\n')
        f.write(str(self.c_SolarMinMaxWn)+'			Solar min, max wn [cm^-1]\n')
        f.write(str(self.c_SlitFunction)+'			Triangular slit function\n')
        f.write(str(self.c_SolarHWHM)+'			Solar HWHM\n')
        f.write(str(self.c_SolarRes)+'			Solar resolution\n')
        f.write(str(self.c_ThermalHWHM)+'			Thermal HWHM\n')
        f.write(str(self.c_ThermalRes)+'			Thermal resolution\n')
        f.write(str(self.c_TauError)+'			Tau error\n')
        f.write(str(self.c_SingleScatterError)+'			Single scattering error\n')
        f.write(str(self.c_AsymmetryParamError)+'			Asymmetry parameter error\n')
        f.write(str(self.c_AlbedoError)+'			Albedo error\n')
        f.write(str(self.c_OutputFileType)+'			1)ASCII 3) binary no header output\n')
        if trynum == 1:
            f.write(self.DataOutPath+'vpl_climate_'+self.casename+'\n')
        else:
            f.write(self.DataOutPath+'vpl_climate_'+self.casename+'_tryNum'+str(trynum)+'\n')
        f.write('2			Overwrite existing files\n')

        f.close()

    
    ### Purpose: Make the vpl climate run file for 2 column mode based on the climate settings
    #
    ## Attribute Dependencies:
    # molecule_dict - All the molecules of interest
    # casename - for naming convention
    # MMW - the atmospheric MMW found by photochem
    # vplclimate_RunScriptDir - the output location of climate runscripts
    # AtmProfPath - Path to the PT and Mixing ratio profiles
    # planetary_gravity - gravity [in m/s]
    # planetary_radius - planet radius [in km]
    # surface_temp - Surface temperature [K] gets set when calling degrade_PT()
    #
    ## Fxn-specific Inputs:
    # trynum - the iteration number youre on for the specific case, defined by self.num_climate_runs
    ##
    def make_2column_climate_runscript(self, trynum=1):

        # Load default climate settings
        self.set_climate_settings()

        # Start a new runscript to create
        f = open(self.vplclimate_RunScriptDir+'RunVPLClimate_2column_'+self.casename+'.script', 'w')

        f.write(str(self.c_NumberTimesteps)+'			Number of timesteps\n')
        f.write(str(self.c_TimestepLength)+'			Timestep length [s]\n')
        f.write(str(self.c_SubstepIntervals)+'			Substep intervals\n')
        f.write(str(self.c_TimestepOutputIntervals)+'			Timestep output intervals\n')
        f.write(str(self.c_TimeStepMethodIndex)+'			Index of time-stepping method\n')
        f.write(str(self.c_TempChangeTolerance)+'			Temperature change tolerance\n')
        f.write(str(self.c_DoubledRadiationGrid)+'			Doubled radiation grid\n')
        f.write('4			HRT calc type index [two-column day-night]\n') # Andrew runs as latitude, then 4 separate versions (at each zenith angle)
        f.write(str(self.c_NumberSolarZeniths)+'			No. SZAs [solar zenith angles] used in avg\n') # maybe just start with 1 when doing sweeps (just the 60 deg)
        f.write(str(self.c_IncludeRadiativeHrt)+'			Include radiative hrt\n')
        f.write(str(self.c_IncludeConvectiveHrt)+'			Include convective hrt\n')
        f.write(str(self.c_IncludeConductiveHrt)+'			Include conductive heating rates\n')
        f.write(str(self.planetary_mass)+'			Mass of planet [kg]\n')
        f.write(str(self.c_DayLength)+'			Length of day [s]\n')
        f.write(str(self.c_YearLength)+'			Length of year [days]\n')
        f.write(str(self.c_OrbitalCalcType)+'			Orbital calc type [fixed orbital distance]\n')
        f.write(str(self.c_SemiMajorAxis)+'			Distance from Star [AU]\n')
        f.write(str(self.planetary_gravity)+'			Surface gravity [m/s^2]\n')
        f.write(str(self.planetary_radius)+'			Planet radius [km]\n')
        ### Dayside:
        f.write(str(self.MMW)+'			atm mean mol wgt [g/mol]\n')
        f.write(str(self.c_NumberMajorGases)+'			No. major atm gases\n')
        f.write(str(self.c_MajorAbsorbingGas)+'			***GAS ABSORBER INDEX*** [o2]\n')
        f.write('1.0			Mixing ratio\n') ### This is currently hard coded
        f.write(str(self.c_PressureJacobians)+'			Pressure jacobians? [0 = None; 1 = Radiance; 2 = Flux]\n')
        f.write(str(self.c_TempJacobians)+'			Temperature jacobians? [0 = None; 1 = Radiance; 2 = Flux]\n')
        f.write(str(self.c_Fractional_dtemp)+'			Fractional dtemp\n')
        f.write(str(self.c_SolarTolerance)+'			Solar tolerance [0-1]\n')
        f.write(str(self.c_ThermalTolerance)+'			Thermal tolerance [0-1]\n')
        f.write('3			List directed input\n')
        f.write(self.AtmProfPath+'PT_profile_dayside_'+self.casename+'.pt\n')
        f.write('1			Lines to skip\n')
        f.write('1,2			columns of P,T\n')
        f.write('100000.			Scale to Pa\n')
        f.write(str(self.surface_temp_dayside)+'			Surface temperature [K]\n') 
        ###
        ### Nightside:
        f.write(str(self.MMW)+'			atm mean mol wgt [g/mol]\n')
        f.write(str(self.c_NumberMajorGases)+'			No. major atm gases\n')
        f.write(str(self.c_MajorAbsorbingGas)+'			***GAS ABSORBER INDEX*** [o2]\n')
        f.write('1.0			Mixing ratio\n') ### This is currently hard coded
        f.write(str(self.c_PressureJacobians)+'			Pressure jacobians? [0 = None; 1 = Radiance; 2 = Flux]\n')
        f.write(str(self.c_TempJacobians)+'			Temperature jacobians? [0 = None; 1 = Radiance; 2 = Flux]\n')
        f.write(str(self.c_Fractional_dtemp)+'			Fractional dtemp\n')
        f.write(str(self.c_SolarTolerance)+'			Solar tolerance [0-1]\n')
        f.write(str(self.c_ThermalTolerance)+'			Thermal tolerance [0-1]\n')
        f.write('3			List directed input\n')
        f.write(self.AtmProfPath+'PT_profile_nightside_'+self.casename+'.pt\n')
        f.write('1			Lines to skip\n')
        f.write('1,2			columns of P,T\n')
        f.write('100000.			Scale to Pa\n')
        f.write(str(self.surface_temp_nightside)+'			Surface temperature [K]\n') 
        ###
        f.write(str(self.c_InternalSurfaceFlux)+'			Internal surface flux [W/m2]\n')
        f.write('0.0         Day-night surface transport heat flux [W/m2]\n')
        f.write(str(self.c_ConvectiveType)+'			Convective type [1 = adjustment; 2 = Mixing length scheme; 3 = turbulent; 4 = moist mixing]\n')
        f.write(str(self.c_MixingLengthType)+'			Mixing length type [1 = fixed; 2 = prop to scale height; 3 = Blackadar aymptotic ML]\n')
        f.write(str(self.c_MixingLengthProportionality)+'			Mixing length proportionality\n')
        f.write(str(self.c_MinEddyDiffusivity)+'			Minimum eddy diffusivity [m2/s]\n')
        f.write(str(self.c_SurfaceWindSpeed)+'			Surface Wind Speed [m/s]\n')
        f.write(str(self.c_SurfRoughnessHeight)+'			Surface roughness height [m]\n')
        f.write(str(self.c_NumberCondensibles)+'			No. condensibles\n')
        ### advection schtuff
        f.write('5			Advection type [1 = const diffusion]\n')
        f.write('1.0			Advection constant / diffusion parameter\n')
        ###
        f.write('1			Generate fluxes using SMART\n')

        #### Dayside Gases here

        f.write(str(len(self.molecule_dict['Gas_names']))+'			******NO. ABSORBING GASES******\n')
        
        # Gases Get written here --------------------

        for m in self.molecule_dict['Gas_names']:
            f.write(str(self.molecule_dict[m])+'			************HITRAN GAS CODE ['+m+']*************\n')
            f.write('0			Mixing ratio jacobians [0 = None; 1 = Radiance; 2 = Flux]\n')

            no_abs_coef = 1
            use_xsec = False
            if m in ['O2', 'H2']: # Then CIA is also available
                no_abs_coef += 1
            if os.path.exists(self.xsec_Path+m.lower()+'xsec.dat'): # If an xsec file is available use that
                use_xsec = True
                no_abs_coef += 1

            f.write(str(no_abs_coef)+'			No. of abs coeff types\n')
            f.write('1			HITRAN Line Absorbers\n') # Always start with the line by line
            f.write('2			Dont compute absorption coefficients\n')
            f.write(self.LBLABC_AbsFilesDir+m+'_'+self.casename+'.abs\n')
            if m in ['O2', 'H2']: # Add CIA
                f.write('2			Collisionally-induced absorption [CIA] files\n')
                f.write(self.xsec_Path+m.lower()+'-'+m.lower()+'_abs.cia\n')
            if use_xsec == True: # Use cross sections
                f.write('3			UV-VIS absorption cross sections\n')
                f.write(self.xsec_Path+m.lower()+'xsec.dat\n')
            
            # Now pass the mixing ratios
            f.write('3			List directed input\n')
            f.write(self.AtmProfPath+'MixingRs_'+self.casename+'.dat\n')
            f.write('1 			Lines to skip\n')
            f.write('1			Use pressure coordinate\n')
            f.write('1,'+str(self.molecule_dict[m+'_RmixCol'])+'			columns of P,rmix\n')
            f.write('1			Mixing type [1 = volume; 2 = mass]\n')
            f.write('100000.,1.0			scale P and rmix\n')

        # End Gases block -------------------------

        f.write(str(self.c_NumberAerosols)+'			******NO. AEROSOLS******\n')
        f.write(str(self.c_SurfaceType)+'			Lambertian surface\n')
        f.write(str(self.c_AlbedoJacobians)+'			No albedo jacobians\n')
        f.write('3			List directed input\n')
        f.write(self.c_SurfaceProfile+'\n')
        f.write(str(self.c_SurfProfile_SkipLines)+'			Lines to skip\n')
        f.write(str(self.c_SurfProfile_wlalbedo_cols)+'			columns of wl, albedo\n')
        f.write(str(self.c_SurfProfile_wlType)+'			Wavelength\n')
        f.write(str(self.c_Convert_SurfProfilewl_microns)+'			Convert to microns\n')
        f.write(str(self.c_ScaleAlbedo)+'			Scale albedo\n')

        ######

        ###### Nightside Gases here:
        
        f.write(str(len(self.molecule_dict['Gas_names']))+'			******NO. ABSORBING GASES******\n')
        
        # Gases Get written here --------------------

        for m in self.molecule_dict['Gas_names']:
            f.write(str(self.molecule_dict[m])+'			************HITRAN GAS CODE ['+m+']*************\n')
            f.write('0			Mixing ratio jacobians [0 = None; 1 = Radiance; 2 = Flux]\n')

            no_abs_coef = 1
            use_xsec = False
            if m in ['O2', 'H2']: # Then CIA is also available
                no_abs_coef += 1
            if os.path.exists(self.xsec_Path+m.lower()+'xsec.dat'): # If an xsec file is available use that
                use_xsec = True
                no_abs_coef += 1

            f.write(str(no_abs_coef)+'			No. of abs coeff types\n')
            f.write('1			HITRAN Line Absorbers\n') # Always start with the line by line
            f.write('2			Dont compute absorption coefficients\n')
            f.write(self.LBLABC_AbsFilesDir+m+'_'+self.casename+'.abs\n')
            if m in ['O2', 'H2']: # Add CIA
                f.write('2			Collisionally-induced absorption [CIA] files\n')
                f.write(self.xsec_Path+m.lower()+'-'+m.lower()+'_abs.cia\n')
            if use_xsec == True: # Use cross sections
                f.write('3			UV-VIS absorption cross sections\n')
                f.write(self.xsec_Path+m.lower()+'xsec.dat\n')
            
            # Now pass the mixing ratios
            f.write('3			List directed input\n')
            f.write(self.AtmProfPath+'MixingRs_'+self.casename+'.dat\n')
            f.write('1 			Lines to skip\n')
            f.write('1			Use pressure coordinate\n')
            f.write('1,'+str(self.molecule_dict[m+'_RmixCol'])+'			columns of P,rmix\n')
            f.write('1			Mixing type [1 = volume; 2 = mass]\n')
            f.write('100000.,1.0			scale P and rmix\n')

        # End Gases block -------------------------

        f.write(str(self.c_NumberAerosols)+'			******NO. AEROSOLS******\n')
        f.write(str(self.c_SurfaceType)+'			Lambertian surface\n')
        f.write(str(self.c_AlbedoJacobians)+'			No albedo jacobians\n')
        f.write('3			List directed input\n')
        f.write(self.c_SurfaceProfile+'\n')
        f.write(str(self.c_SurfProfile_SkipLines)+'			Lines to skip\n')
        f.write(str(self.c_SurfProfile_wlalbedo_cols)+'			columns of wl, albedo\n')
        f.write(str(self.c_SurfProfile_wlType)+'			Wavelength\n')
        f.write(str(self.c_Convert_SurfProfilewl_microns)+'			Convert to microns\n')
        f.write(str(self.c_ScaleAlbedo)+'			Scale albedo\n')

        ############

        f.write(str(self.c_NumberStreams)+'			No. of Streams\n')
        f.write(str(self.c_HRTSources)+'			Sources: 1) Solar, 2) thermal, 3) both\n')
        f.write('3			List directed input\n')
        f.write(self.c_StellarSpectrum+'\n')
        f.write(str(self.c_StellarSpect_SkipLines)+'			Lines to skip\n')
        f.write(str(self.c_SolarFluxUnits)+'			Solar flux units\n')
        f.write(str(self.c_SolarSpectralUnits)+'			Solar spectral units\n')
        f.write(str(self.c_Convert_Stellar_microns)+'			Conversion factor to microns\n')
        f.write(str(self.c_StellarSpect_wnflux_col)+'			Cols of wn and flux\n')
        f.write(str(self.c_OutputsType)+'			Outputs [fluxes]\n')
        f.write(str(self.c_NumberOutputAzimuths)+'			No. of output azimuth angles\n')
        f.write(str(self.c_Azimuth)+'			Azimuth\n')
        f.write(str(self.c_OutputsLevel)+'			Outputs [TOA]\n')
        f.write(str(self.c_OutputsUnits)+'			Output units [radiance]\n')
        f.write(str(self.c_ThermalMinMaxWn)+'			Thermal min, max wn [cm^-1]\n')
        f.write(str(self.c_SolarMinMaxWn)+'			Solar min, max wn [cm^-1]\n')
        f.write(str(self.c_SlitFunction)+'			Triangular slit function\n')
        f.write(str(self.c_SolarHWHM)+'			Solar HWHM\n')
        f.write(str(self.c_SolarRes)+'			Solar resolution\n')
        f.write(str(self.c_ThermalHWHM)+'			Thermal HWHM\n')
        f.write(str(self.c_ThermalRes)+'			Thermal resolution\n')
        f.write(str(self.c_TauError)+'			Tau error\n')
        f.write(str(self.c_SingleScatterError)+'			Single scattering error\n')
        f.write(str(self.c_AsymmetryParamError)+'			Asymmetry parameter error\n')
        f.write(str(self.c_AlbedoError)+'			Albedo error\n')
        f.write(str(self.c_OutputFileType)+'			1)ASCII 3) binary no header output\n')
        if trynum == 1:
            f.write(self.DataOutPath+'vpl_2column_climate_'+self.casename+'\n')
        else:
            f.write(self.DataOutPath+'vpl_2column_climate_'+self.casename+'_tryNum'+str(trynum)+'\n')
        f.write('2			Overwrite existing files\n')

        f.close()

    ### Purpose: Make the SMART run file based on the SMART settings
    #
    ## Attribute Dependencies:
    # molecule_dict - All the molecules of interest
    # casename - for naming convention
    # MMW - the atmospheric MMW found by photochem
    # SMART_RunScriptDir - the output location of SMART runscripts
    # AtmProfPath - Path to the PT and Mixing ratio profiles
    # planetary_gravity - gravity [in m/s]
    # planetary_radius - planet radius [in km]
    # surface_temp - Surface temperature [K] gets set when calling degrade_PT()
    #
    ## Fxn-specific Inputs:
    # None
    ##
    def make_smart_runscript(self, whichcol=None):

        # Load default climate settings
        self.set_climate_settings()
        self.set_smart_settings()

        # Start a new runscript to create
        if whichcol == None:
            f = open(self.SMART_RunScriptDir+'RunSMART_'+self.casename+'.run', 'w')
        elif whichcol == 'dayside':
            f = open(self.SMART_RunScriptDir+'RunSMART_dayside_'+self.casename+'.run', 'w')
        elif whichcol == 'nightside':
            f = open(self.SMART_RunScriptDir+'RunSMART_nightside_'+self.casename+'.run', 'w')

        f.write(str(self.s_PressJacobians)+'			Pressure Jacobians (0-None)\n')
        f.write(str(self.s_TempJacobians)+'			Temperature Jacobians (0-None)\n')
        
        # PT profile should remain the same
        f.write('3			Formatted Atmospheric Structure File\n')
        if whichcol == None:
            f.write(self.AtmProfPath+'PT_profile_'+self.casename+'.pt\n')
        elif whichcol == 'dayside':
            f.write(self.AtmProfPath+'PT_profile_dayside_'+self.casename+'.pt\n')
        elif whichcol == 'nightside':
            f.write(self.AtmProfPath+'PT_profile_nightside_'+self.casename+'.pt\n')
        f.write('1			Lines to skip\n')
        f.write('1,2			columns of P,T\n')
        f.write('100000.			Scale to Pa\n')
        if whichcol == None:
            f.write(str(self.surface_temp)+'			Surface temperature [K]\n')
        elif whichcol == 'dayside':
            f.write(str(self.surface_temp_dayside)+'			Surface temperature [K]\n')
        elif whichcol == 'nightside':
            f.write(str(self.surface_temp_nightside)+'			Surface temperature [K]\n')

        # Gases Get written here --------------------

        f.write(str(len(self.molecule_dict['Gas_names']))+'			Number of Gas Absorbers\n')
        for m in self.molecule_dict['Gas_names']:
            f.write(str(self.molecule_dict[m])+'			Hitran gas code for '+m+'\n')
            f.write('0			Absorber jacobians [0 = None; 1 = Radiance; 2 = Flux]\n')

            no_abs_coef = 1
            use_xsec = False
            if m in ['O2', 'H2']: # Then CIA is also available
                no_abs_coef += 1
            if os.path.exists(self.xsec_Path+m.lower()+'xsec.dat'): # If an xsec file is available use that
                use_xsec = True
                no_abs_coef += 1

            f.write(str(no_abs_coef)+'			No. of abs coeff types\n')
            f.write('1			HITRAN Line Absorbers\n') # Always start with the line by line
            if whichcol == None:
                f.write(self.LBLABC_AbsFilesDir+m+'_'+self.casename+'.abs\n')
            elif whichcol == 'dayside':
                f.write(self.LBLABC_AbsFilesDir+m+'_d_'+self.casename+'.abs\n')
            elif whichcol == 'nightside':
                f.write(self.LBLABC_AbsFilesDir+m+'_n_'+self.casename+'.abs\n')
            if m in ['O2', 'H2']: # Add CIA
                f.write('2			Collisionally-induced absorption [CIA] files\n')
                f.write(self.xsec_Path+m.lower()+'-'+m.lower()+'_abs.cia\n')
            if use_xsec == True: # Use cross sections
                f.write('3			UV-VIS absorption cross sections\n')
                f.write(self.xsec_Path+m.lower()+'xsec.dat\n')
            
            # Now pass the mixing ratios
            f.write('3			List directed input\n')
            f.write(self.AtmProfPath+'MixingRs_'+self.casename+'.dat\n')
            f.write('1 			Lines to skip\n')
            f.write('1			Use pressure coordinate\n')
            f.write('1,'+str(self.molecule_dict[m+'_RmixCol'])+'			columns of P,rmix\n')
            f.write('1			Mixing type [1 = volume; 2 = mass]\n')
            f.write('100000.,1.0			scale P and rmix\n')

        # End Gases block -------------------------

        f.write(str(self.s_NumberAerosols)+'			Number of Aerosols\n')
        f.write(str(self.s_SurfaceType)+'           Lambert Reflection\n')
        f.write(str(self.s_AlbedoJacobians)+'           Albedo Jacobians (0 - None)\n')
        f.write('3           List Directed Surface Albedo Profile\n')
        f.write(str(self.s_SurfaceProfile)+'\n')
        f.write(str(self.s_SurfProfile_SkipLines)+'          Lines to Skip\n')
        f.write(str(self.s_SurfProfile_wlalbedo_cols)+'         Columns of wl and albedo\n')
        f.write(str(self.s_SurfProfile_wlType)+'			Wavelength Grid Type\n')
        f.write(str(self.s_Convert_SurfProfilewl_microns)+'			Convert to microns\n')
        f.write(str(self.s_ScaleAlbedo)+'			Scale albedo\n')

        # Physical Parameters
        f.write(str(self.s_SemiMajorAxis)+'   	Distance from Star (AU)\n')
        f.write(str(self.planetary_gravity)+'			Surface Gravity [m/s*s]\n')
        f.write(str(self.planetary_radius)+'	    Planetary Radius [km]\n')
        f.write(str(self.MMW)+'		Mean Molecular Weight\n')

        # Rayleigh scattering params
        f.write(str(self.s_NumRayleigh)+'			Number of Rayleigh Scatterers\n')
        f.write(str(self.s_RayleighIndex)+'			Rayleigh Scatter Index\n')
        f.write(str(self.s_vmrRayleigh)+'			Volume Mixing Ratio\n')

        # Heating Sources / streams params
        f.write(str(self.s_NumberStreams)+'			Number of Streams\n')
        f.write(str(self.s_HRTSources)+'			HRT Source\n')
        f.write('3			File Format - List Directed\n')
        f.write(str(self.s_StellarSpectrum)+'\n')
        f.write(str(self.s_StellarSpect_SkipLines)+'			Lines to Skip\n')
        f.write(str(self.s_SolarFluxUnits)+'			Units solar flux\n')
        f.write(str(self.s_SolarSpectralUnits)+'			Solar spectral units\n')
        if self.s_SolarSpectralUnits != 2:
            f.write(str(self.s_Convert_Stellar_microns)+'			Micron Conversion Factor\n')
        f.write(str(self.s_StellarSpect_wnflux_col)+' 		Columns of wn and Flux\n')
        f.write(str(self.s_NumZenith)+'			Number of Solar Zenith Angles\n')
        f.write(str(self.s_ZenithAzimuth_Angles)+'			Zenith and Azimuth Angles\n')
        f.write(str(self.s_ConvergenceCriteria)+'			Convergence Criterion\n')

        # Output info
        f.write(str(self.s_OutputFormat)+'			Output Format\n')
        f.write(str(self.s_TranistType)+'           Ray Tracing Transit\n')
        f.write(str(self.s_StellarRadius)+'      Stellar Radius (Rsun)\n')
        f.write(str(self.s_ImpactParam)+'          Impact Parameter (Rstar)\n')
        f.write('1           Include Refraction\n')
        f.write(str(self.s_LimbDarkening)+'           Include Limb Darkening?\n')
        f.write(str(self.s_NumAzimuths)+'			Number of Azimuth Angles\n')
        f.write(str(self.s_AzimuthAngles)+'			Azimuth Angles\n')
        f.write(str(self.s_OutputLevels)+'			Output levels (1 - TOA only)\n')
        f.write(str(self.s_OutputUnits)+'			Output Radiance Units [2 - W/m*m/sr/micron]\n')
        f.write(str(self.s_MinMax_wavenumber)+'	Min and Max Wavenumber\n')
        f.write(str(self.s_GridType)+'			Type of Grid (2-slit)\n')
        f.write(str(self.s_SpectralResponseFxn)+'			Triangular Spectral Response Function\n')
        f.write(str(self.s_FWHM)+'			FWHM\n')
        f.write(str(self.s_SampleRes)+'			Sampling Resolution [per cm]\n')
        f.write(str(self.s_TauBinError)+'		Error for tau Binning\n')
        f.write(str(self.s_pi0BinError)+'		Error for pi0 Binning\n')
        f.write(str(self.s_gError)+'		g Error\n')
        f.write(str(self.s_AlbedoError)+'		Albedo Error\n')
        f.write(str(self.s_OutputFileFormat)+'			Output Format (ascii)\n')
        if whichcol == None:
            f.write(self.OutPath+self.casename+'_SMART\n')
        elif whichcol == 'dayside':
            f.write(self.OutPath+self.casename+'_dayside_SMART\n')
        elif whichcol == 'nightside':
            f.write(self.OutPath+self.casename+'_nightside_SMART\n')
        f.write('2			Overwrite\n')

        f.close()


    ### Purpose: Read in the out.dist file to create a python dict with its values (e.g., mixing ratios, T, EDD, Ndens, etc)
    ####  used for finding new surface pressures
    #
    ## Attribute Dependencies:
    # photochem_InputsDir - The directory of inputs using the species.dat file
    # photochemDir - To retrieve the out.dist file
    # nlevel_fine - The NZ of photochem
    #
    ## Fxn-specific Inputs:
    # None
    ##
    def ingest_outdist(self):

        # Read in species.dat file
        species = self.photochem_InputsDir+'species.dat'
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
        blockend = self.nlevel_fine

        NQblocks = np.ceil(NQ/10)
        species_ind = 0

        # Loop through out.dist file NQ blocks and append species mixing ratios to dictionary
        for i in range(int(NQblocks)):
            curr_nq_block = ascii.read(self.photochemDir+'OUTPUT/out.dist', data_start=blockstart, data_end=blockend, header_start=None)
            blockstart = blockend
            blockend = blockend+self.nlevel_fine

            for i in curr_nq_block.columns:
                d[new_gases[species_ind]] = list(curr_nq_block[i])
                species_ind += 1

        # Add Temp, Edd, and number density to output dictionary 
        T_edd_block = ascii.read(self.photochemDir+'OUTPUT/out.dist', data_start=blockstart, data_end=blockend)
        blockstart = blockend

        d['Temp'] = list(T_edd_block['col1'])
        d['Edd'] = list(T_edd_block['col2'])
        d['NDens'] = list(T_edd_block['col3'])

        return d
    

    ### Purpose: Create a new in.dist file for an updated surface pressure / number density grid
    #
    ## Attribute Dependencies:
    # photochem_InputsDir - The directory of inputs using the species.dat and in.dist files
    # photochemDir - To retrieve the out.dist file
    # nlevel_fine - The NZ of photochem
    #
    ## Fxn-Specific Inputs:
    # newndens_perspecies - The new/updated number densities per species (as a python dictionary)... 
    #     ...created within self.change_atmospheric_pressure() using the New Pressure helper fxn get_true_number_densities()
    # newtotndens - The new/updated total number density created within self.change_atmospheric_pressure()...
    #     ...using the New Pressure helper fxn new_total_Ndens()
    # newmixings - The new/updated mixing ratios per species (as a python dictionary)...
    #     ...created within self.change_atmospheric_pressure() using the New Pressure helper fxn new_mixing_rats()
    # oldoutdict - The ingested out.dist file as a python dictionary (from self.ingest_outdist() called by self.change_atmospheric_pressure())
    ##
    def new_indist_new_pressure(self, newndens_perspecies, newtotndens, newmixings, oldoutdict):

        # Define needed files
        oldout = self.photochemDir+'OUTPUT/out.dist'

        fnew = open(self.photochem_InputsDir+'NEWpressure_in.dist', 'w')

        # Read in species.dat file
        species = self.photochem_InputsDir+'species.dat'
        nsp = open(species, 'r')

        done = False
        new_gases = []
        for l in nsp.readlines():
            if l.split()[0][0] == '*':
                if done == True:
                    break
            else:
                new_gases.append(l.split()[0])
                done = True

        NQ = len(new_gases)
        NQblocks = np.ceil(NQ/10)

        for i in range(int(NQblocks)):
            gases_to_add = list(newmixings.keys())[i*10:i*10+10]
            write_table = Table()
            counter = 1
            for gas in gases_to_add:
                write_table.add_column(newmixings[gas], name='col'+str(counter))
                counter += 1

            for line in range(len(write_table)):
                fnew.write('   ')
                for col in write_table.columns:
                    if write_table[col][line] < 9e-99:
                        val = 9e-99
                    else:
                        val = write_table[col][line]
                    fnew.write("{:.8E}".format(val)+'   ')
                fnew.write('\n')

        write_table = Table()
        write_table.add_column(oldoutdict['Temp'], name='col1')
        write_table.add_column(oldoutdict['Edd'], name='col2')
        write_table.add_column(newtotndens, name='col3')
        write_table.add_column(newmixings['O3'], name='col4')
        write_table.add_column(newndens_perspecies['CO2'], name='col5')

        for line in range(len(write_table)):
            fnew.write('   ')
            for col in write_table.columns:
                if write_table[col][line] < 9e-99:
                    val = 9e-99
                else:
                    val = write_table[col][line]
                fnew.write("{:.8E}".format(val)+'   ')
            fnew.write('\n')

        lastbloc_start = (self.nlevel_fine*NQblocks)+self.nlevel_fine

        olddist_txt = open(oldout, 'r')
        alllines = olddist_txt.readlines()
        olddist_txt.close()
        for line in range(len(alllines)):
            if line >= lastbloc_start:
                fnew.write(alllines[line])

        fnew.close()

        subprocess.run('rm '+self.photochem_InputsDir+'in.dist', shell=True)
        subprocess.run('mv '+self.photochem_InputsDir+'NEWpressure_in.dist '+self.photochem_InputsDir+'in.dist', shell=True)

    ### Purpose: If the atmospheric pressure needs to be changed/adjusted, this function finds the new pressure... 
    ###     ...and updates the in.dist and PLANET.dat files
    # 
    ## Attribute Dependencies:
    # photochem_InputsDir - The directory of photochem inputs
    # photochemDir - Photochem directory (for OUTPUTS/)
    # nlevel_fine - The NZ of photochem 
    # planetary_gravity - Planet gravity (to calculate pressure)
    # MMW - mean molecular weight
    #
    ## Fxn-specific Inputs:
    # None
    #  
    ##
    def change_atmospheric_pressure(self, after_sgbsl_err=False):

        # Load the current out.dist values
        outdistdic = self.ingest_outdist()

        # sum of all mixing ratios to find new total number density
        new_total_VMR = sum_mixing_ratios(outdistdic, self.nlevel_fine, self.n2mixingrat)

        # New number densities for each species
        new_Ndens_species = get_true_number_densities(outdistdic)

        # New total number density and the change in number density on each level
        new_Ndens_tot, change_in_Ndens = new_total_Ndens(new_total_VMR, outdistdic['NDens'])
        
        # Update in.dist file for new number densities
        new_VMR_species = new_mixing_rats(new_Ndens_species, new_Ndens_tot)
        if after_sgbsl_err == False: # only update if it wasn't following an error
            self.new_indist_new_pressure(new_Ndens_species, new_Ndens_tot, new_VMR_species, outdistdic)
        else:
            if self.fixsgbsl == False:
                self.new_indist_new_pressure(new_Ndens_species, new_Ndens_tot, new_VMR_species, outdistdic)

        # Find new surface pressure
        loaded_ptz_out = ascii.read(self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist')
        colmass = find_tot_column_mass_dens(new_Ndens_tot, loaded_ptz_out['ALT'], self.MMW)
        new_surfP = (colmass*(self.planetary_gravity*(u.m*u.s**-2).to(u.cm*u.s**-2))*100)*u.Pa.to(u.bar) # 100 converts per cm to per m

        # Check if the pressure converged or if photochem needs to be rerun with new pressure
        maxchange = np.abs(self.updated_atm_pressure - new_surfP)/self.updated_atm_pressure
        if maxchange <= self.NewPressure_Psurf_tolerance:
            pressure_converged = True
            self.updated_atm_pressure = new_surfP 
        else:
            pressure_converged = False
            if after_sgbsl_err == True:
                self.updated_atm_pressure = 2*self.updated_atm_pressure
            else:
                self.updated_atm_pressure = (self.updated_atm_pressure + new_surfP)/2

        # Change surface pressure in PLANET.dat
        planetdat_new = open(self.photochem_InputsDir+'New_PLANET.dat', 'w')
        planetdat_old = open(self.photochem_InputsDir+'PLANET.dat', 'r')

        planetdat_lines = planetdat_old.readlines()
        for l in planetdat_lines:
            hold = l.split()
            if 'surface' in hold and 'pressure' in hold:
                planetdat_new.write("{:.2e}".format(self.updated_atm_pressure))
                planetdat_new.write(' = P0, surface pressure [bar] \n')
            else:
                if self.num_climate_runs == 0 or self.force_dz_adjust == True: # Only change if we don't have dzg already
                    if 'DZGRID' in hold:
                        if self.updated_atm_pressure >= 1 and self.updated_atm_pressure < 10:
                            self.dzg = 0.75e5
                            planetdat_new.write("{:.2E}".format(self.dzg)+' = DZGRID [cm] \n')
                        elif self.updated_atm_pressure > 10:
                            # Only run after photochem has run at least 1x
                            ptzcurr = ascii.read(self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist')
                            presscurr = ptzcurr['PRESS']
                            if presscurr[199] > 1e-4:
                                self.dzg = self.dzg*2.5
                                planetdat_new.write("{:.2E}".format(self.dzg)+' = DZGRID [cm] \n')
                            elif presscurr[199] > 1e-5:
                                self.dzg = self.dzg*2
                                planetdat_new.write("{:.2E}".format(self.dzg)+' = DZGRID [cm] \n')
                            elif presscurr[199] < 1e-10:
                                self.dzg = self.dzg*0.5
                                planetdat_new.write("{:.2E}".format(self.dzg)+' = DZGRID [cm] \n')
                            else:
                                self.dzg = self.dzg
                                planetdat_new.write("{:.2E}".format(self.dzg)+' = DZGRID [cm] \n')
                        else:
                            self.dzg = 5e4
                            planetdat_new.write("{:.2E}".format(self.dzg)+' = DZGRID [cm] \n')
                    else:
                        planetdat_new.write(l)
                else:
                    planetdat_new.write(l)

        planetdat_new.close()
        planetdat_old.close()

        subprocess.run('rm '+self.photochem_InputsDir+'PLANET.dat', shell=True)
        subprocess.run('mv '+self.photochem_InputsDir+'New_PLANET.dat '+self.photochem_InputsDir+'PLANET.dat', shell=True)

        return pressure_converged, maxchange, new_surfP


    ### Purpose: Change the T/EDD profiles from climate in photochems in.dist to rerun photochem after a climate run
    ##
    ## Attribute Dependencies:
    # photochem_InputsDir - to get NZ and NQ from the parameters.inc file 
    #
    ## Fxn-specific Inputs:
    # oldptz - path to the PTZ_out file from the last photochem run (use PHOTOCHEM/OUTPUT/ dir)
    # climateprof - the last output profile from vpl climate, found with self.get_final_climate_output_temp_profile()
    ##
    def update_indist_T_EDD(self, oldptz, climateprof):

        # WARNING: this may be hard coded if folks don't write their parameters.inc file the same way as templates etc
        # I don't think it will be a problem but just noting in case
        # This retrieves the NQ and NZ from parameters.inc
        NQFile = open(self.photochem_InputsDir+'parameters.inc', 'r')
        NQFileLines = NQFile.readlines()
        for l in NQFileLines:
            if len(l.split('NQ=')) > 1:
                NQ = l.split('NQ=')[1]
                NQ = int(NQ.split(',')[0])
                NZ = l.split('NZ=')[1]
                NZ = int(NZ.split(',')[0])
                break

        # Find number of blocks of mixing ratios until T/EDD columns
        NQblocks = np.ceil(NQ/10)

        # Read old in.dist file to change
        #distfi = open(olddist, 'r')

        # Read ptz file to interpolate pressure grid
        ptzold = ascii.read(oldptz)
        photochem_pressure = ptzold['PRESS']

        ### Interpolation

        # Put in surface to TOA order following photochem
        climate_press = np.array(climateprof['P[Pa]'])
        climate_T = np.array(climateprof['T[K]'])
        climate_edd = np.array(climateprof['Km[m2/s]'])
        #climate_h2o = np.array(climateprof['rmix[kg/kg]'])

        # convert to bar
        climate_press = climate_press * 1e-5

        # convert to vol mixing ratio
        #climate_h2o = climate_h2o*(self.MMW/18.01)

        # convert EDD to whatever andrew says
        climate_edd = climate_edd*1e4

        # interpolation to put in in.dist
        new_T = np.interp(photochem_pressure[::-1], climate_press, climate_T)[::-1]
        new_edd = np.interp(photochem_pressure[::-1], climate_press, climate_edd)[::-1]
        #new_h2o = np.interp(photochem_pressure[::-1], climate_press, climate_h2o)[::-1]

        ### Now to create new in.dist

        fnew = open(self.photochem_InputsDir+'NEWin.dist', 'w')

        # Keep track of the lines that are starting and ending the current block in the in.dist file
        blockstart = 0
        blockend = NZ

        # this should handle all the mixing ratio blocks
        for i in range(int(NQblocks)):
            curr_nq_block = ascii.read(self.photochem_InputsDir+'in.dist', data_start=blockstart, data_end=blockend)
            blockstart = blockend
            blockend = blockend+NZ

            # Actually decided not to change H2O
            #if i == block_w_h2o:
            #    curr_nq_block['col'+str(H2OCol_inblock)] = new_h2o

            for line in range(len(curr_nq_block)):
                fnew.write('   ')
                for col in curr_nq_block.columns:
                    fnew.write("{:.8E}".format(curr_nq_block[col][line])+'   ')
                fnew.write('\n')

        # now the T/EDD/DEN/O3/CO2 block
        T_edd_block = ascii.read(self.photochem_InputsDir+'in.dist', data_start=blockstart, data_end=blockend)
        blockstart = blockend
        T_edd_block['col1'] = new_T
        T_edd_block['col2'] = new_edd # Need to convert from m**2/ to cm**2/s; climate model will be between 1/2 and 1000, photochem will be 1e4-1e6 ish

        for line in range(len(T_edd_block)):
            fnew.write('   ')
            for col in T_edd_block.columns:
                fnew.write("{:.8E}".format(T_edd_block[col][line])+'   ')
            fnew.write('\n')

        # now the rest of in.dist
        olddist_txt = open(self.photochem_InputsDir+'in.dist', 'r')
        alllines = olddist_txt.readlines()
        olddist_txt.close()
        for line in range(len(alllines)):
            if line >= blockstart:
                fnew.write(alllines[line])

        fnew.close()

        subprocess.run('rm '+self.photochem_InputsDir+'in.dist', shell=True)
        subprocess.run('mv '+self.photochem_InputsDir+'NEWin.dist '+self.photochem_InputsDir+'in.dist', shell=True)

    ### 

    ### Purpose: Update the inert N2 mixing ratio in species.dat
    ##
    ## Attribute Dependencies:
    # self.photochem_InputsDir - to find the species.dat file to update
    #
    ## Fxn-specific Inputs:
    # n2mixingrat - new mixing ratio to use
    ##
    def update_N2_mixingrat_speciesdat(self):

        species_new = open(self.photochem_InputsDir+'New_species.dat', 'w')
        species_old = open(self.photochem_InputsDir+'species.dat', 'r')

        lines = species_old.readlines()
        for l in lines:
            if len(l.split()) > 1:
                if l.split()[0] == 'N2':
                    species_new.write('N2         IN  0 0 0 0 2 0    '+"{:.3E}".format(self.n2mixingrat)+'\n')
                else:
                    species_new.write(l)
            else:
                species_new.write(l)

        species_new.close()
        species_old.close()
        
        # Delete old species and rename fixed version to be species.dat
        subprocess.run('rm '+self.photochem_InputsDir+'species.dat', shell=True)
        subprocess.run('mv '+self.photochem_InputsDir+'New_species.dat '+self.photochem_InputsDir+'species.dat', shell=True)


    ### Purpose: Update the planet.dat total pressure if N2 was meant to be much greater than the total 
    ##
    ## Attribute Dependencies:
    # self.photochem_InputsDir - to find the planet.dat file to update
    #
    ## Fxn-specific Inputs:
    # 
    ##
    def update_N2_totalpressure_planetdat(self):

        # Updating planet dat:
        planetdat_new = open(self.photochem_InputsDir+'New_PLANET.dat', 'w')
        planetdat_old = open(self.photochem_InputsDir+'PLANET.dat', 'r')

        planetdat_lines = planetdat_old.readlines()
        for l in planetdat_lines:
            hold = l.split()
            if 'surface' in hold and 'pressure' in hold:
                planetdat_new.write("{:.2e}".format(self.updated_atm_pressure))
                planetdat_new.write(' = P0, surface pressure [bar] \n')
            else:
                if self.num_climate_runs == 0: # Only change if we don't have dzg already
                    if 'DZGRID' in hold:
                        if self.updated_atm_pressure >= 1:
                            self.dzg = 0.75e5
                            planetdat_new.write("{:.2E}".format(self.dzg)+' = DZGRID [cm] \n')
                        else:
                            self.dzg = 5e4
                            planetdat_new.write("{:.2E}".format(self.dzg)+' = DZGRID [cm] \n')
                    else:
                        planetdat_new.write(l)
                else:
                    planetdat_new.write(l)

        planetdat_new.close()
        planetdat_old.close()

        subprocess.run('rm '+self.photochem_InputsDir+'PLANET.dat', shell=True)
        subprocess.run('mv '+self.photochem_InputsDir+'New_PLANET.dat '+self.photochem_InputsDir+'PLANET.dat', shell=True)

    ### AUTOMATIC PIPELINE, calls all above functions sequentially defined by flow chart, updates necessary object parameters
    ##
    ## Dependent on all Attributes, no fxn-specific inputs
    ##
    def run_automatic(self):
        
        # Prepare your backup directory for photochem data
        if self.BackupPhotochemRuns == True:
            if not os.path.exists(self.photochemBackupDir):
                os.mkdir(self.photochemBackupDir)
            else:
                subprocess.run('rm -rf '+self.photochemBackupDir+'*', shell=True)
        # Prepare your directory for storing lblabc .abs files
        if not os.path.exists(self.LBLABC_AbsFilesDir):
            os.mkdir(self.LBLABC_AbsFilesDir)
        #else:
        #    subprocess.run('rm -rf '+self.LBLABC_AbsFilesDir+'*', shell=True)
        # Prepare directory for storing lblabc run script files
        if not os.path.exists(self.lblabc_RunScriptDir):
            os.mkdir(self.lblabc_RunScriptDir)
        # Prepare directory for storing vpl climate run script files
        if not os.path.exists(self.vplclimate_RunScriptDir):
            os.mkdir(self.vplclimate_RunScriptDir)
        # Prepare directory for storing new photochem inputs
        #if not os.path.exists(self.photochem_InputsDir):
        #    os.mkdir(self.photochem_InputsDir)
        self.setup_intial_photochem_dir()
        # Make sure model and data output dirs exist
        if not os.path.exists(self.OutPath):
            os.mkdir(self.OutPath)
        if not os.path.exists(self.DataOutPath):
            os.mkdir(self.DataOutPath)
        if not os.path.exists(self.SMART_RunScriptDir):
            os.mkdir(self.SMART_RunScriptDir)
        # Prepare the Hyak environment
        #self.prepare_hyak_env()

        if self.verbose == True:
            if not os.path.exists(self.OutPath+self.casename+'_SavingInfoOut.txt'):
                ftestingoutput = open(self.OutPath+self.casename+'_SavingInfoOut.txt', 'w')
            else:
                ftestingoutput = open(self.OutPath+self.casename+'_SavingInfoOut.txt', 'a')


        ##########################
        # Find the atmospheric pressure used 
        self.updated_atm_pressure = -1
        self.dzg = 0.75e5 # Just a default in case

        # Get original pressure and update DZGrid
        if self.adjust_atmospheric_pressure == True:
            planetdat_new = open(self.photochem_InputsDir+'New_PLANET.dat', 'w')
        planetdat_old = open(self.photochem_InputsDir+'PLANET.dat', 'r')

        planetdat_lines = planetdat_old.readlines()
        for l in planetdat_lines:
            hold = l.split()
            if 'surface' in hold and 'pressure' in hold:
                self.updated_atm_pressure = float(hold[0])
                if self.adjust_atmospheric_pressure == True:
                    planetdat_new.write(l)
            else:
                if self.adjust_atmospheric_pressure == True:
                    if 'DZGRID' in hold:
                        if self.updated_atm_pressure >= 1 and self.updated_atm_pressure < 10:
                            self.dzg = 0.75e5
                            planetdat_new.write("{:.2E}".format(self.dzg)+' = DZGRID [cm] \n')
                        elif self.updated_atm_pressure >= 10:
                            self.dzg = 2e5
                            planetdat_new.write("{:.2E}".format(self.dzg)+' = DZGRID [cm] \n')
                        else:
                            self.dzg = 5e4
                            planetdat_new.write("{:.2E}".format(self.dzg)+' = DZGRID [cm] \n')
                    else:
                        planetdat_new.write(l)

        if self.adjust_atmospheric_pressure == True:
            planetdat_new.close()
        planetdat_old.close()

        if self.adjust_atmospheric_pressure == True:
            subprocess.run('rm '+self.photochem_InputsDir+'PLANET.dat', shell=True)
            subprocess.run('mv '+self.photochem_InputsDir+'New_PLANET.dat '+self.photochem_InputsDir+'PLANET.dat', shell=True)
        ##########################

        ##########################
        # If N2 is supposed to be adjusted, find the mixing ratio and set in the species file
        if self.adjust_N2_amount == True:
            self.n2mixingrat = self.N2_fixed_pressure / self.updated_atm_pressure

            # If you need more N2 than pressure available, update the species and PLANET.dat
            if self.n2mixingrat >= 1:
                self.updated_atm_pressure = self.N2_fixed_pressure/0.99
                self.n2mixingrat = 0.99

                self.update_N2_totalpressure_planetdat()

            # Update species file for N2:
            self.update_N2_mixingrat_speciesdat()
        
        else:
            species_old = open(self.photochem_InputsDir+'species.dat', 'r')

            lines = species_old.readlines()
            for l in lines:
                if len(l.split()) > 1:
                    if l.split()[0] == 'N2':
                        self.n2mixingrat = float(l.split()[8])
                        break

        ##########################

        if self.verbose == True:
            ftestingoutput.write('Starting pressure: '+str(self.updated_atm_pressure)+' bars\n')

        climate_subtries = 0

        # Start loop to find global convergence with photochem + lblabc + vpl climate
        while self.global_convergence == False:

            ### Run photochem section start ---------------------

            # Need to adjust DZGrid
            if self.num_climate_runs > 0: 

                if self.force_dz_adjust == False:
                    planetdat_new = open(self.photochem_InputsDir+'New_PLANET.dat', 'w')
                    planetdat_old = open(self.photochem_InputsDir+'PLANET.dat', 'r')

                    planetdat_lines = planetdat_old.readlines()
                    for l in planetdat_lines:
                        hold = l.split()
                        if 'DZGRID' in hold:
                            planetdat_new.write("{:.2E}".format(self.dzg)+' = DZGRID [cm] \n')
                        else:
                            planetdat_new.write(l)

                    planetdat_new.close()
                    planetdat_old.close()

                    subprocess.run('rm '+self.photochem_InputsDir+'PLANET.dat', shell=True)
                    subprocess.run('mv '+self.photochem_InputsDir+'New_PLANET.dat '+self.photochem_InputsDir+'PLANET.dat', shell=True)
                
                else:
                    planetdat_new = open(self.photochem_InputsDir+'New_PLANET.dat', 'w')
                    planetdat_old = open(self.photochem_InputsDir+'PLANET.dat', 'r')

                    planetdat_lines = planetdat_old.readlines()
                    for l in planetdat_lines:
                        hold = l.split()
                        if 'DZGRID' in hold:
                            planetdat_new.write("{:.2E}".format(self.dzg)+' = DZGRID [cm] \n')
                        else:
                            planetdat_new.write(l)

                    planetdat_new.close()
                    planetdat_old.close()

                    subprocess.run('rm '+self.photochem_InputsDir+'PLANET.dat', shell=True)
                    subprocess.run('mv '+self.photochem_InputsDir+'New_PLANET.dat '+self.photochem_InputsDir+'PLANET.dat', shell=True)

            self.num_photochem_runs += 1
            if self.verbose == True:
                #print('----> Beginning photochem run Try number '+str(self.num_photochem_runs))
                ftestingoutput.write('----> Beginning photochem run Try number '+str(self.num_photochem_runs)+'\n')
            self.run_photochem_1instance(CleanMake=True, InputCopy=self.photochem_InputsDir, trynum=self.num_photochem_runs)
            
            ### If this is the first run, Retrieve the surf gravity and radius of the planet from PLANET.dat
            if self.num_photochem_runs == 1:
                planet = open(self.photochemDir+'INPUTFILES/PLANET.dat', 'r')
                planet_lines = planet.readlines()
                planet.close()
                grav = None
                rad = None
                for i in planet_lines:
                    if len(i.split('= G')) > 1:
                        grav = float(i.split()[0])*(u.cm*u.s**-2).to(u.m*u.s**-2) #*1e-2 # Get the gravity from the first line of PLANET.dat and convert to m/s**2 (should be originally cm/s**2)
                    elif len(i.split('= R0')) > 1:
                        rad = float(i.split()[0])*u.cm.to(u.km) #*1e-5 # Get radius from PLANET.dat and convert from cm to km
                    elif grav != None and rad != None:
                        break
                    
                # Set object values
                self.planetary_gravity = grav
                self.planetary_radius = rad

            photochem_subtries = 1
            # Currently, the trynum will only refer to the converged case (should we save every single try no matter what?)
            local_photochem_conv, grosserr, l2err, finaltime, nsteps_photo, sgbslerror = self.check_photochem_conv(trynum=self.num_photochem_runs)

            if sgbslerror == True and local_photochem_conv == False:
                if self.verbose == True:
                    #print('SGBSL Error occured in photochem run')
                    ftestingoutput.write('\n SGBSL Error occured in photochem run\n')
                
                if self.fixsgbsl == True:
                    # If SGBSL error occured and youre allowed to update pressure, try updating pressure and running to convergence
                    if self.adjust_atmospheric_pressure == True:
                        # Need to set MMW early 
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

                        pressure_converged, maxchange, holdnewsurfp = self.change_atmospheric_pressure(after_sgbsl_err=True)

                        # If N2 needs to be adjusted:
                        if self.adjust_N2_amount == True:
                            self.n2mixingrat = self.N2_fixed_pressure / self.updated_atm_pressure

                            # If you need more N2 than pressure available, update the species and PLANET.dat
                            if self.n2mixingrat >= 1:
                                self.updated_atm_pressure = self.N2_fixed_pressure/0.99
                                self.n2mixingrat = 0.99

                                self.update_N2_totalpressure_planetdat()

                            # Update species file for N2:
                            self.update_N2_mixingrat_speciesdat()

                        if self.verbose == True:
                            #print('Attempting to adjust pressure to fix SGBSL error')
                            #print('New Pressure would be: '+str(holdnewsurfp)+' bars, using pressure of: '+str(self.updated_atm_pressure)+' bars')
                            ftestingoutput.write('Attempting to adjust pressure to fix SGBSL error\n')
                            ftestingoutput.write('New Pressure would be: '+str(holdnewsurfp)+' bars, using pressure of: '+str(self.updated_atm_pressure)+' bars\n\n')


            # If photochem did not converge, try try again
            while local_photochem_conv == False:
                
                if self.num_photochem_runs == 1:
                    subprocess.run('cp '+self.OutPath+'photochem_run_output_'+self.casename+'.run '+self.OutPath+'photochem_run_output_'+self.casename+'_subtry'+str(photochem_subtries)+'.run', shell=True)
                    subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.OutPath+'PTZMix_subtry'+str(photochem_subtries)+'.dist', shell=True)
                    subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.OutPath+'outdist_subtry'+str(photochem_subtries)+'.dist', shell=True)
                    subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.OutPath+'outout_subtry'+str(photochem_subtries)+'.out', shell=True)
                else:
                    subprocess.run('cp '+self.OutPath+'photochem_run_output_'+self.casename+'_Try'+str(self.num_photochem_runs)+'.run '+self.OutPath+'photochem_run_output_'+self.casename+'_Try'+str(self.num_photochem_runs)+'_subtry'+str(photochem_subtries)+'.run', shell=True)
                    subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.OutPath+'PTZMix_Try'+str(self.num_photochem_runs)+'_subtry'+str(photochem_subtries)+'.dist', shell=True)
                    subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.OutPath+'outdist_Try'+str(self.num_photochem_runs)+'_subtry'+str(photochem_subtries)+'.dist', shell=True)                
                    subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.OutPath+'outout_Try'+str(self.num_photochem_runs)+'_subtry'+str(photochem_subtries)+'.out', shell=True)

                if photochem_subtries > self.max_iterations_master:
                    break
                    
                # If you had an SGBSL error that is being fixed by pressure change, need a clean make of photochem
                if self.fixsgbsl == True and sgbslerror == True and self.adjust_atmospheric_pressure == True:
                    self.run_photochem_1instance(CleanMake=True, InputCopy=self.photochem_InputsDir, trynum=self.num_photochem_runs)
                
                else:
                    subprocess.run('rm -rf '+self.photochemDir+'in.dist', shell=True)
                    subprocess.run('rm -rf '+self.photochemDir+'PTZ_mixingratios_in.dist', shell=True)
                    subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.photochemDir+'in.dist', shell=True)
                    self.run_photochem_1instance(CleanMake=False, InputCopy=False, trynum=self.num_photochem_runs)

                photochem_subtries += 1
                local_photochem_conv, grosserr, l2err, finaltime, nsteps_photo, sgbslerror = self.check_photochem_conv(trynum=self.num_photochem_runs, subtries=photochem_subtries, prevtime=finaltime)

                if self.verbose == True:
                    ftestingoutput.write('Normalized Gross error: '+str(grosserr)+'\n')
                    ftestingoutput.write('L2 Error: '+str(l2err)+'\n')
                    ftestingoutput.write('Time of final timestep: '+str(finaltime)+'\n')
                    if sgbslerror == True:
                        ftestingoutput.write('SGBSL Error occured\n')
                    if local_photochem_conv == False:
                        ftestingoutput.write('Photochem subtry '+str(photochem_subtries)+' NOT converged\n\n')

                # If local convergence cant be found, might need to try adjusting the pressure early (could just be << 100% and running into issues)
                if photochem_subtries in [50, 100, 150] and self.adjust_atmospheric_pressure == True:
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
                    pressure_converged, maxchange, holdnewsurfp = self.change_atmospheric_pressure()

                    if self.verbose == True:
                        ftestingoutput.write('Subtries reached '+str(photochem_subtries)+', attempting to adjust pressure \n')
                        ftestingoutput.write('New Pressure: '+str(holdnewsurfp)+', using '+str(self.updated_atm_pressure)+' Bars \n\n')
                        ftestingoutput.close()
                        ftestingoutput = open(self.OutPath+self.casename+'_SavingInfoOut.txt', 'a')


                if self.fixsgbsl == True and sgbslerror == True and self.adjust_atmospheric_pressure == True:
                    # Need to reset MMW early 
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
                    pressure_converged, maxchange, holdnewsurfp = self.change_atmospheric_pressure(after_sgbsl_err=True)

                    # If N2 needs to be adjusted:
                    if self.adjust_N2_amount == True:
                        self.n2mixingrat = self.N2_fixed_pressure / self.updated_atm_pressure

                        # If you need more N2 than pressure available, update the species and PLANET.dat
                        if self.n2mixingrat >= 1:
                            self.updated_atm_pressure = self.N2_fixed_pressure/0.99
                            self.n2mixingrat = 0.99

                            self.update_N2_totalpressure_planetdat()

                        # Update species file for N2:
                        self.update_N2_mixingrat_speciesdat()

                    if self.verbose == True:
                        #print('Attempting to adjust pressure to fix SGBSL error')
                        #print('New Pressure: '+str(holdnewsurfp)+' bars, using pressure of: '+str(self.updated_atm_pressure)+' bars')
                        ftestingoutput.write('Attempting to adjust pressure to fix SGBSL error\n')
                        ftestingoutput.write('New Pressure would be: '+str(holdnewsurfp)+' bars, using pressure of: '+str(self.updated_atm_pressure)+' bars\n\n')



            # Copy the out.dist to be the new in.dist in the photochem inputs directory
            subprocess.run('rm -rf '+self.photochem_InputsDir+'in.dist', shell=True)
            subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.photochem_InputsDir+'in.dist', shell=True)

            # If convergence failed, raise io error or break run
            if local_photochem_conv == False and self.suppress_IOerrors == False:
                subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.DataOutPath+'FINAL_out_FAILED.dist', shell=True)
                subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.DataOutPath+'FINAL_out_FAILED.out', shell=True)
                subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.DataOutPath+'FINAL_PTZ_mixingratios_out_FAILED.dist', shell=True)
                raise IOError('Photochem ran >'+str(self.max_iterations_master)+' times with no convergence. On photochem run number '+str(self.num_photochem_runs))
            elif local_photochem_conv == False and self.suppress_IOerrors == True:
                subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.DataOutPath+'FINAL_out_FAILED.dist', shell=True)
                subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.DataOutPath+'FINAL_out_FAILED.out', shell=True)
                subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.DataOutPath+'FINAL_PTZ_mixingratios_out_FAILED.dist', shell=True)
                if self.verbose == True:
                    ftestingoutput.write('Max Iterations reached ('+str(self.max_iterations_master)+'), couldnt find new pressure, ending run\n')
                    #print('Max Iterations reached ('+str(self.max_iterations_master)+'), couldnt find new pressure, ending run')
                break

            if self.verbose == True:
                #print('Photochem local convergence found with '+str(photochem_subtries)+' subtries')
                ftestingoutput.write('Photochem local convergence found with '+str(photochem_subtries)+' subtries\n')

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

            # If the atmospheric pressure should be checked / updated, run that
            photochem_newPsurf_subtries = 0
            if self.adjust_atmospheric_pressure == True:
                pressure_converged, maxchange, holdnewsurfp = self.change_atmospheric_pressure()

                # If N2 needs to be adjusted:
                if self.adjust_N2_amount == True and pressure_converged == False:
                    self.n2mixingrat = self.N2_fixed_pressure / self.updated_atm_pressure
                    n2converged = True

                    # If you need more N2 than pressure available, update the species and PLANET.dat
                    if self.n2mixingrat >= 1:
                        self.updated_atm_pressure = self.N2_fixed_pressure/0.99
                        self.n2mixingrat = 0.99
                        n2converged = False

                        self.update_N2_totalpressure_planetdat()

                    # Update species file for N2:
                    self.update_N2_mixingrat_speciesdat()

                if self.verbose == True:
                    #print('New Pressure found: '+"{:.4e}".format(holdnewsurfp)+' bars')
                    ftestingoutput.write('New Pressure found: '+"{:.4e}".format(holdnewsurfp)+' bars\n')
                    if pressure_converged == True:
                        #print('Pressure converged, no need to rerun photochem')
                        ftestingoutput.write('Pressure converged, no need to rerun photochem\n')
                        if self.adjust_N2_amount == True:
                            if n2converged == False:
                                ftestingoutput.write('BUT N2 pressure larger than total, rerunning\n')
                                pressure_converged = False
                
                if pressure_converged == False:

                    if self.verbose == True:
                        #print('Pressure NOT converged, rerunning photochem using '+str(self.updated_atm_pressure)+' bars')
                        ftestingoutput.write('Pressure NOT converged, rerunning photochem using '+str(self.updated_atm_pressure)+' bars \n\n')

                    while pressure_converged == False:

                        if photochem_newPsurf_subtries > self.max_iterations_master:
                            break

                        self.run_photochem_1instance(CleanMake=True, InputCopy=self.photochem_InputsDir, trynum=self.num_photochem_runs)

                        photochem_newPsurf_subtries += 1
                        # Currently, the trynum will only refer to the converged case (should we save every single try no matter what?)
                        local_photochem_conv, grosserr, l2err, finaltime, nsteps_photo, sgbslerror = self.check_photochem_conv(trynum=self.num_photochem_runs)
                        photochem_newPsurf_inner_subtries = 1

                        # If SGBSL Error occured, just set local conv to True so it skips the next bit and tries to find a different pressure 
                        if sgbslerror == True:
                            if self.verbose == True:
                                #print('SGBSL Error Occured, retrying for pressure')
                                ftestingoutput.write('SGBSL Error Occured, retrying for pressure\n')
                            if self.fixsgbsl == True:
                                local_photochem_conv = True

                        # If photochem did not converge, try try again
                        while local_photochem_conv == False:
                            if photochem_newPsurf_inner_subtries > self.max_iterations_master:
                                break

                            if self.num_photochem_runs == 1:
                                subprocess.run('cp '+self.OutPath+'photochem_run_output_'+self.casename+'.run '+self.OutPath+'photochem_run_output_'+self.casename+'_Psurfsubtry'+str(photochem_newPsurf_subtries)+'_Innertry_'+str(photochem_newPsurf_inner_subtries)+'.run', shell=True)
                                subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.OutPath+'outout_Psurfsubtry'+str(photochem_newPsurf_subtries)+'_Innertry_'+str(photochem_newPsurf_inner_subtries)+'.out', shell=True)
                            else:
                                subprocess.run('cp '+self.OutPath+'photochem_run_output_'+self.casename+'_Try'+str(self.num_photochem_runs)+'.run '+self.OutPath+'photochem_run_output_'+self.casename+'_Try'+str(self.num_photochem_runs)+'_Psurfsubtry'+str(photochem_newPsurf_subtries)+'_Innertry_'+str(photochem_newPsurf_inner_subtries)+'.run', shell=True)
                                subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.OutPath+'outout_Try'+str(self.num_photochem_runs)+'_Psurfsubtry'+str(photochem_newPsurf_subtries)+'_Innertry_'+str(photochem_newPsurf_inner_subtries)+'.out', shell=True)

                            subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.OutPath+'PTZMix_Psurfsubtry'+str(photochem_newPsurf_subtries)+'_Innertry_'+str(photochem_newPsurf_inner_subtries)+'.dist', shell=True)
                            subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.OutPath+'outdist_Psurfsubtry'+str(photochem_newPsurf_subtries)+'_Innertry_'+str(photochem_newPsurf_inner_subtries)+'.dist', shell=True) 

                            subprocess.run('rm -rf '+self.photochemDir+'in.dist', shell=True)
                            subprocess.run('rm -rf '+self.photochemDir+'PTZ_mixingratios_in.dist', shell=True)
                            subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.photochemDir+'in.dist', shell=True)
                            self.run_photochem_1instance(CleanMake=False, InputCopy=False, trynum=self.num_photochem_runs)
                            photochem_newPsurf_inner_subtries += 1
                            local_photochem_conv, grosserr, l2err, finaltime, nsteps_photo, sgbslerror = self.check_photochem_conv(trynum=self.num_photochem_runs, subtries=photochem_newPsurf_inner_subtries, prevtime=finaltime)

                            # If SGBSL error occured, break and pick a new pressure
                            if sgbslerror == True:
                                if self.verbose == True:
                                    ftestingoutput.write('SGBSL error occured in inner convergence loop on inner subtry '+str(photochem_newPsurf_inner_subtries)+'\n')
                                    ftestingoutput.write('Breaking loop to pick new pressure\n')
                                if self.fixsgbsl == True:
                                    break

                        if self.num_photochem_runs == 1:
                            subprocess.run('cp '+self.OutPath+'photochem_run_output_'+self.casename+'.run '+self.OutPath+'photochem_run_output_'+self.casename+'_Psurfsubtry'+str(photochem_newPsurf_subtries)+'_Innertry_'+str(photochem_newPsurf_inner_subtries)+'.run', shell=True)
                            subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.OutPath+'outout_Psurfsubtry'+str(photochem_newPsurf_subtries)+'_Innertry_'+str(photochem_newPsurf_inner_subtries)+'.out', shell=True)
                        else:
                            subprocess.run('cp '+self.OutPath+'photochem_run_output_'+self.casename+'_Try'+str(self.num_photochem_runs)+'.run '+self.OutPath+'photochem_run_output_'+self.casename+'_Try'+str(self.num_photochem_runs)+'_Psurfsubtry'+str(photochem_newPsurf_subtries)+'_Innertry_'+str(photochem_newPsurf_inner_subtries)+'.run', shell=True)
                            subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.OutPath+'outout_Try'+str(self.num_photochem_runs)+'_Psurfsubtry'+str(photochem_newPsurf_subtries)+'_Innertry_'+str(photochem_newPsurf_inner_subtries)+'.out', shell=True)
                        
                        subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.OutPath+'PTZMix_Psurfsubtry'+str(photochem_newPsurf_subtries)+'_Innertry_'+str(photochem_newPsurf_inner_subtries)+'.dist', shell=True)
                        subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.OutPath+'outdist_Psurfsubtry'+str(photochem_newPsurf_subtries)+'_Innertry_'+str(photochem_newPsurf_inner_subtries)+'.dist', shell=True) 


                        if photochem_newPsurf_inner_subtries > self.max_iterations_master:
                                break
                        
                        # Need to reset MMW 
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

                        pressure_converged, maxchange, holdnewsurfp = self.change_atmospheric_pressure(after_sgbsl_err=sgbslerror)

                        # If N2 needs to be adjusted:
                        if self.adjust_N2_amount == True and pressure_converged == False:
                            self.n2mixingrat = self.N2_fixed_pressure / self.updated_atm_pressure

                            # If you need more N2 than pressure available, update the species and PLANET.dat
                            if self.n2mixingrat >= 1:
                                self.updated_atm_pressure = self.N2_fixed_pressure/0.99
                                self.n2mixingrat = 0.99
                                ftestingoutput.write('N2 too large, updated pressure should be:'+str(self.updated_atm_pressure)+'\n')
                                
                                self.update_N2_totalpressure_planetdat()
                                ftestingoutput.close()
                                ftestingoutput = open(self.OutPath+self.casename+'_SavingInfoOut.txt', 'a')


                            # Update species file for N2:
                            self.update_N2_mixingrat_speciesdat()

                        if pressure_converged == True and sgbslerror == True:
                            pressure_converged = False

                        if pressure_converged == False and self.verbose == True:
                            ftestingoutput.write('Normalized Gross error: '+str(grosserr)+'\n')
                            ftestingoutput.write('L2 Error: '+str(l2err)+'\n')
                            ftestingoutput.write('Time of final timestep: '+str(finaltime)+'\n')
                            ftestingoutput.write('Max change of a number density layer: '+str(maxchange)+'\n')
                            ftestingoutput.write('New Pressure: '+str(holdnewsurfp)+', using '+str(self.updated_atm_pressure)+' Bars \n')
                            #print('Surf Pressure Subtry '+str(photochem_newPsurf_subtries)+' NOT converged')
                            ftestingoutput.write('Surf Pressure Subtry '+str(photochem_newPsurf_subtries)+' NOT converged\n\n')

                    if photochem_newPsurf_inner_subtries > self.max_iterations_master and self.suppress_IOerrors == False:
                        subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.DataOutPath+'FINAL_out_FAILED.dist', shell=True)
                        subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.DataOutPath+'FINAL_out_FAILED.out', shell=True)
                        subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.DataOutPath+'FINAL_PTZ_mixingratios_out_FAILED.dist', shell=True)
                        raise IOError('Photochem attempted to use new pressure >'+str(self.max_iterations_master)+' times with no inner convergence. On photochem run number '+str(self.num_photochem_runs))
                    elif photochem_newPsurf_inner_subtries > self.max_iterations_master and self.suppress_IOerrors == True:
                        subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.DataOutPath+'FINAL_out_FAILED.dist', shell=True)
                        subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.DataOutPath+'FINAL_out_FAILED.out', shell=True)
                        subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.DataOutPath+'FINAL_PTZ_mixingratios_out_FAILED.dist', shell=True)
                        if self.verbose == True:
                            #print('Max iterations reached ('+str(self.max_iterations_master)+'), using new pressure without finding inner convergence, ending run')
                            ftestingoutput.write('Max iterations reached ('+str(self.max_iterations_master)+'), using new pressure without finding inner convergence, ending run\n')
                        break

                    if pressure_converged == False and self.suppress_IOerrors == False:
                        subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.DataOutPath+'FINAL_out_FAILED.dist', shell=True)
                        subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.DataOutPath+'FINAL_out_FAILED.out', shell=True)
                        subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.DataOutPath+'FINAL_PTZ_mixingratios_out_FAILED.dist', shell=True)
                        raise IOError('Photochem attempted to find new pressure >'+str(self.max_iterations_master)+' times with no pressure convergence. On photochem run number '+str(self.num_photochem_runs))
                    elif pressure_converged == False and self.suppress_IOerrors == True:
                        subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.DataOutPath+'FINAL_out_FAILED.dist', shell=True)
                        subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.DataOutPath+'FINAL_out_FAILED.out', shell=True)
                        subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.DataOutPath+'FINAL_PTZ_mixingratios_out_FAILED.dist', shell=True)
                        if self.verbose == True:
                            #print('Max iterations reached ('+str(self.max_iterations_master)+'), couldnt find new pressure, ending run')
                            ftestingoutput.write('Max iterations reached ('+str(self.max_iterations_master)+'), couldnt find new pressure, ending run\n')
                        break

                    if self.verbose == True:
                        #print('Pressure converged after '+str(photochem_newPsurf_subtries)+' iterations, with '+str(photochem_newPsurf_inner_subtries)+' number of photochem reruns at this pressure')
                        ftestingoutput.write('Pressure converged after '+str(photochem_newPsurf_subtries)+' iterations, with '+str(photochem_newPsurf_inner_subtries)+' number of photochem reruns at this pressure\n')
                        ftestingoutput.write('Converged pressure: '+str(self.updated_atm_pressure)+' bars\n')

            if self.verbose == True:
                ftestingoutput.close()
                ftestingoutput = open(self.OutPath+self.casename+'_SavingInfoOut.txt', 'a')
            # Save backup of photochem output if desired
            if self.BackupPhotochemRuns == True:
                self.backup_photochem_run(trynum=self.num_photochem_runs)

            # If you only had 1 subtry, and its not the first run, check for global convergence
            #if photochem_subtries == 1 and self.num_photochem_runs != 1 and photochem_newPsurf_subtries == 0:
            if photochem_subtries <= 3 and self.num_photochem_runs != 1 and photochem_newPsurf_subtries == 0 and climate_subtries == 1:
                #if nsteps_photo < 1000:
                self.global_convergence = True
                if self.verbose == True:
                    #print('Global Convergence achieved')
                    ftestingoutput.write('Global Convergence achieved')
                subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.DataOutPath+'FINAL_out.dist', shell=True)
                subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.DataOutPath+'FINAL_out.out', shell=True)
                subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.DataOutPath+'FINAL_PTZ_mixingratios_out.dist', shell=True)
                break
            else:
                self.global_convergence = False
                if self.num_photochem_runs > self.max_iterations_master or self.num_climate_runs > self.max_iterations_climate:
                    subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.DataOutPath+'FINAL_out_FAILED.dist', shell=True)
                    subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.DataOutPath+'FINAL_out_FAILED.out', shell=True)
                    subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.DataOutPath+'FINAL_PTZ_mixingratios_out_FAILED.dist', shell=True)

                    if self.suppress_IOerrors == False:
                        raise IOError('Photochem+Climate have run together >'+str(self.max_iterations_master)+' without finding global convergence, run failed.')
                    elif self.suppress_IOerrors == True:
                        break

                #self.global_convergence = True

            # If MCMC is only looking for a pressure convergence (for computational efficiency), just find convergence here
            if self.MCMC_pressure_only == True:
                if sgbslerror == False:
                    subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.DataOutPath+'FINAL_out.dist', shell=True)
                    subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.DataOutPath+'FINAL_out.out', shell=True)
                    subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.DataOutPath+'FINAL_PTZ_mixingratios_out.dist', shell=True)
                    self.global_convergence = True
                break

            ### Run photochem section end --------------------------

            ### Create degraded atmospheric profile to prepare for LBLABC and Climate -------------
            if self.adjust_atmospheric_pressure == True:
                if self.updated_atm_pressure < 1e-2:
                    self.degrade_PT(grid_spacing='log')
                    if self.verbose == True:
                        ftestingoutput.write('Log Spacing Used\n')
                        #print('log spacing used')
                else:
                    self.degrade_PT()
            else:
                self.degrade_PT()
            if self.verbose == True:
                #print('Degraded PT profile created from photochem run '+str(self.num_photochem_runs))
                ftestingoutput.write('Degraded PT profile created from photochem run '+str(self.num_photochem_runs)+'\n')
            ### Degraded PT Profile finished ------------

            ### Create new mixing ratios profile file --------------------------
            self.prep_rmix_file(self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist')
            ### Mixing ratios profile created --------------------------

            ### Run LBLABC section start ----------------------------

            # Create the LBLABC run scripts
            # Wont change much after first run, but needs to be done everytime to update MMW
            # Not much computational expense saved by only replacing MMW after each run ...
            # ... So cleaner just to recreate everytime
            self.make_lblabc_runscripts()

            # Now Run LBLABC for all the gases of interest
            for gas in self.molecule_dict['Gas_names']:
                self.run_lblabc_1instance(self.lblabc_RunScriptDir+'RunLBLABC_'+gas+'_'+self.casename+'.script', gas)
                if self.verbose == True:
                    #print('LBLABC run for '+gas+' complete, LBLABC iteration '+str(self.num_lblabc_runs+1))
                    ftestingoutput.write('LBLABC run for '+gas+' complete, LBLABC iteration '+str(self.num_lblabc_runs+1)+'\n')
            self.num_lblabc_runs += 1

            ### Run LBLABC section end -------------------------------------
            
            ### Run VPL Climate section start ------------------------------

            # Add 1 to the climate run counter
            self.num_climate_runs += 1

            # First Generate the climate runscript
            # Things, e.g., MMW need to be updated, so just do this every time
            self.make_climate_runscript(trynum=self.num_climate_runs)

            if self.verbose == True:
                #print('Climate Runscript created, beginning first climate run')
                ftestingoutput.write('Climate Runscript created, beginning first climate run\n')

            # Now run climate 
            ftestingoutput.write('The runscript: '+self.vplclimate_RunScriptDir+'RunVPLClimate_'+self.casename+'.script\n')
            ftestingoutput.write('The executable: '+self.vplclimate_executable+'\n')
            self.run_climate_1instance(self.vplclimate_RunScriptDir+'RunVPLClimate_'+self.casename+'.script', self.vplclimate_executable, trynum=self.num_climate_runs)

            if self.verbose == True:
                #print('First Climate run completed')
                ftestingoutput.write('First Climate run completed\n')

            # Check for local convergence of climate, similar to process followed for a given try on photochem            
            local_climate_convergence, tropheating, avgflux = self.check_vplclimate_conv(trynum=self.num_climate_runs)

            if self.verbose == True:
                if local_climate_convergence == True:
                    #print('Climate convergence found on first try for run number '+str(self.num_climate_runs))
                    ftestingoutput.write('Climate convergence found on first try for run number '+str(self.num_climate_runs)+'\n')
                else:
                    #print('Climate convergence NOT found on first try for run number '+str(self.num_climate_runs)+', beginning rerun sequence')
                    ftestingoutput.write('Climate convergence NOT found on first try for run number '+str(self.num_climate_runs)+', beginning rerun sequence\n')

            climate_subtries = 1
            # Until climate converges, loop through taking new temp profile
            while local_climate_convergence == False:

                if climate_subtries == self.max_iterations_climate:
                    break

                climate_subtries += 1

                # First, get the final profile from the last climate run to use in the restart
                climate_profile = self.get_final_climate_output_temp_profile(trynum=self.num_climate_runs)

                # Update temperature in the PT profile and update the surface temperature
                self.replace_PT_tempcol(climate_profile['T[K]'])

                # Recreate the runscript to update the surface temp
                self.make_climate_runscript(trynum=self.num_climate_runs)

                if self.num_climate_runs == 1:
                    subprocess.run('mv '+self.OutPath+'vpl_climate_output_'+self.casename+'.run '+self.OutPath+'vpl_climate_output_'+self.casename+'_Subtry'+str(climate_subtries)+'.run', shell=True)
                else:
                    subprocess.run('mv '+self.OutPath+'vpl_climate_output_'+self.casename+'_Try'+str(self.num_climate_runs)+'.run '+self.OutPath+'vpl_climate_output_'+self.casename+'_Try'+str(self.num_climate_runs)+'_Subtry'+str(climate_subtries)+'.run', shell=True)

                # Re run Climate 
                if self.verbose == True:
                    #print('Beginning Climate rerun')
                    ftestingoutput.write('Beginning Climate rerun\n')
                self.run_climate_1instance(self.vplclimate_RunScriptDir+'RunVPLClimate_'+self.casename+'.script', self.vplclimate_executable, trynum=self.num_climate_runs)
                if self.verbose == True:
                    #print('Climate subtry number '+str(climate_subtries)+' completed')
                    ftestingoutput.write('Climate subtry number '+str(climate_subtries)+' completed\n')

                # Check convergence
                local_climate_convergence, tropheating, avgflux = self.check_vplclimate_conv(trynum=self.num_climate_runs, subtries=climate_subtries)

                if self.verbose == True:
                    if local_climate_convergence == True:
                        #print('Climate convergence found on subtry number '+str(climate_subtries)+' for run number '+str(self.num_climate_runs))
                        ftestingoutput.write('Climate convergence found on subtry number '+str(climate_subtries)+' for run number '+str(self.num_climate_runs)+'\n')
                    else:
                        #print('Climate convergence NOT found on subtry number '+str(climate_subtries)+' for run number '+str(self.num_climate_runs)+', continuing rerun sequence')
                        ftestingoutput.write('Climate convergence NOT found on subtry number '+str(climate_subtries)+' for run number '+str(self.num_climate_runs)+', continuing rerun sequence\n')

            if local_climate_convergence == False and self.suppress_IOerrors == False:
                subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.DataOutPath+'FINAL_out_FAILED.dist', shell=True)
                subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.DataOutPath+'FINAL_out_FAILED.out', shell=True)
                subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.DataOutPath+'FINAL_PTZ_mixingratios_out_FAILED.dist', shell=True)
                raise IOError('Climate could not converge in >'+str(self.max_iterations_master)+' tries. For climate run number '+str(self.num_climate_runs))
            elif local_climate_convergence == False and self.suppress_IOerrors == True:
                subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.DataOutPath+'FINAL_out_FAILED.dist', shell=True)
                subprocess.run('cp '+self.photochemDir+'OUTPUT/out.out '+self.DataOutPath+'FINAL_out_FAILED.out', shell=True)
                subprocess.run('cp '+self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist '+self.DataOutPath+'FINAL_PTZ_mixingratios_out_FAILED.dist', shell=True)
                break

            ### Run VPL Climate section end ------------------------------

            ### Update in.dist for next photochem run ------------------------------

            # If this is the first change, save the original in.dist just in case
            if self.num_photochem_runs == 1:
                subprocess.run('cp '+self.photochem_InputsDir+'in.dist '+self.photochem_InputsDir+'Original_in.dist', shell=True)
            
            # Now update in.dist
            climate_profile = self.get_final_climate_output_temp_profile(trynum=self.num_climate_runs)
            #ftestingoutput.write('Retrieved final climate output\n')
            self.update_indist_T_EDD(self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist', climate_profile)
            #ftestingoutput.write('Updated indist TEDD\n')

            ### Update in.dist section end ------------------------------

        if self.verbose == True:
            ftestingoutput.write('Global Convergence, '+str(self.global_convergence)+'\n')
            ftestingoutput.write('Continuing 2 col? '+str(self.clim2col_restarting)+'\n')
            ftestingoutput.write('runspec: '+str(self.run_spectra)+', include 2 col: '+str(self.include_2column_climate)+'\n')
            ftestingoutput.close()


        if self.global_convergence == True and (self.run_spectra == True or self.include_2column_climate == True): 

            if self.verbose == True:
                ftestingoutput = open(self.OutPath+self.casename+'_SavingInfoOut.txt', 'a')

            ##### Generate SMART spectra of the final converged atmosphere ------------------------------

            #'''
            if self.clim2col_restarting == False:
                ### Create degraded atmospheric profile to prepare for LBLABC and Climate -------------
                if self.adjust_atmospheric_pressure == True:
                    if self.updated_atm_pressure < 1e-2:
                        self.degrade_PT(grid_spacing='log')
                        if self.verbose == True:
                            ftestingoutput.write('log spacing used\n')
                    else:
                        self.degrade_PT()
                else:
                    self.degrade_PT()
                if self.verbose == True:
                    #print('Degraded PT profile created from photochem run '+str(self.num_photochem_runs))
                    ftestingoutput.write('Degraded PT profile created from photochem run '+str(self.num_photochem_runs)+'\n')
                ### Degraded PT Profile finished ------------
                
                
                ### Create new mixing ratios profile file --------------------------
                self.prep_rmix_file(self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist')
                ### Mixing ratios profile created --------------------------
                
                
                ### Rerun the LBLABC files for the most recent atmosphere ------------------------------
                
            if self.clim2col_restarting == True:

                planet = open(self.photochem_InputsDir+'PLANET.dat', 'r')
                planet_lines = planet.readlines()
                planet.close()
                grav = None
                rad = None
                for i in planet_lines:
                    if len(i.split('= G')) > 1:
                        grav = float(i.split()[0])*(u.cm*u.s**-2).to(u.m*u.s**-2) #*1e-2 # Get the gravity from the first line of PLANET.dat and convert to m/s**2 (should be originally cm/s**2)
                    elif len(i.split('= R0')) > 1:
                        rad = float(i.split()[0])*u.cm.to(u.km) #*1e-5 # Get radius from PLANET.dat and convert from cm to km
                    elif grav != None and rad != None:
                        break
                    
                # Set object values
                self.planetary_gravity = grav
                self.planetary_radius = rad
                
                f = open(self.vplclimate_RunScriptDir+'RunVPLClimate_2column_'+self.casename+'.script')
                lines = f.readlines()
                for l in lines:
                    if len(l.split('atm mean mol wgt [g/mol]')) > 1:
                        self.MMW = float(l.split()[0])
                        break


            self.make_lblabc_runscripts()

            if self.MultiNest_DataFit == True: # To run LBLABC gases in parallel

                lblabc_input = []
                for gas in self.molecule_dict['Gas_names']:
                    lblabc_input.append([self.lblabc_RunScriptDir+'RunLBLABC_'+gas+'_'+self.casename+'.script', gas, 'Avg'])
                
                with Pool() as p:
                    lblruns = p.map(self.run_lblabc_1instance_Parallel, lblabc_input)
            
            else:

                # Now Run LBLABC for all the gases of interest
                for gas in self.molecule_dict['Gas_names']:
                    self.run_lblabc_1instance(self.lblabc_RunScriptDir+'RunLBLABC_'+gas+'_'+self.casename+'.script', gas)
                    if self.verbose == True:
                        #print('LBLABC run for '+gas+' complete, LBLABC iteration '+str(self.num_lblabc_runs+1))
                        ftestingoutput.write('LBLABC run for '+gas+' complete, LBLABC iteration '+str(self.num_lblabc_runs+1)+'\n')
                self.num_lblabc_runs += 1
                ftestingoutput.close()
                ftestingoutput = open(self.OutPath+self.casename+'_SavingInfoOut.txt', 'a')
            
                
            ### Rerun the LBLABC Section Finish ------------------------------


            ### Setup and Run 2 column climate to convergence ------------------------------
            if self.include_2column_climate == True:

                if self.clim2col_restarting == False:
                    # create two pt profiles
                    if self.dayside_starting_PT == None:
                        shutil.copyfile(self.AtmProfPath+'PT_profile_'+self.casename+'.pt', self.AtmProfPath+'PT_profile_dayside_'+self.casename+'.pt')
                    else:
                        shutil.copyfile(self.dayside_starting_PT, self.AtmProfPath+'PT_profile_dayside_'+self.casename+'.pt')

                    if self.nightside_starting_PT == None:
                        shutil.copyfile(self.AtmProfPath+'PT_profile_'+self.casename+'.pt', self.AtmProfPath+'PT_profile_nightside_'+self.casename+'.pt')
                    else:
                        shutil.copyfile(self.nightside_starting_PT, self.AtmProfPath+'PT_profile_nightside_'+self.casename+'.pt')

                    
                    # create two surface temps
                    if self.dayside_starting_PT == None and self.nightside_starting_PT == None:
                        self.surface_temp_dayside = self.surface_temp
                        self.surface_temp_nightside = self.surface_temp
                    else:
                        openhold = ascii.read(self.dayside_starting_PT)
                        self.surface_temp_dayside = openhold['Temp'][len(openhold['Temp'])-1]
                        openhold = ascii.read(self.nightside_starting_PT)
                        self.surface_temp_nightside = openhold['Temp'][len(openhold['Temp'])-1]

                    # Add 1 to the climate run counter
                    self.num_2col_climate_runs += 1

                    # First Generate the climate runscript
                    # Things, e.g., MMW need to be updated, so just do this every time
                    self.make_2column_climate_runscript(trynum=self.num_2col_climate_runs)

                    if self.verbose == True:
                        #print('Climate 2 column Runscript created, beginning first 2 column climate run')
                        ftestingoutput.write('Climate 2 column Runscript created, beginning first 2 column climate run\n')

                    # Now run climate 
                    ftestingoutput.write('The runscript: '+self.vplclimate_RunScriptDir+'RunVPLClimate_2column_'+self.casename+'.script\n')
                    ftestingoutput.write('The executable: '+self.vplclimate_executable+'\n')
                    self.run_climate_1instance(self.vplclimate_RunScriptDir+'RunVPLClimate_2column_'+self.casename+'.script', self.vplclimate_executable, trynum=self.num_2col_climate_runs, twocol=True)

                    if self.verbose == True:
                        #print('First 2 column Climate run completed')
                        ftestingoutput.write('First 2 column Climate run completed\n')

                    climate_subtries = 1

                        # Check for local convergence of climate, similar to process followed for a given try on photochem            
                    local_climate_convergence, tropheating_dayside, tropheating_nightside, avgflux, cnvtype = self.check_2column_vplclimate_conv(trynum=self.num_2col_climate_runs, climsubtry=climate_subtries)

                    if self.verbose == True:
                        if local_climate_convergence == True:
                            #print('2 column Climate convergence found on first try for run number '+str(self.num_2col_climate_runs))
                            ftestingoutput.write('2 column Climate convergence found on first try for run number '+str(self.num_2col_climate_runs)+'\n')
                        else:
                            #print('2 column Climate convergence NOT found on first try for run number '+str(self.num_2col_climate_runs)+', beginning rerun sequence')
                            ftestingoutput.write('2 column Climate convergence NOT found on first try for run number '+str(self.num_2col_climate_runs)+', beginning rerun sequence\n')

                    
                
                else:
                    self.num_2col_climate_runs += 1
                    patt = re.compile(r'Subtry(\d+)')
                    maxsubt = -1

                    for filename in os.listdir(self.OutPath):
                        match = patt.search(filename)
                        if match:
                            subtnum = int(match.group(1))
                            if subtnum > maxsubt:
                                maxsubt = subtnum
                    
                    climate_subtries = maxsubt
                    
                    # Re run Climate 
                    if os.path.exists(self.OutPath+'vpl_2col_climate_output_'+self.casename+'.run'):
                        subprocess.run('rm '+self.OutPath+'vpl_2col_climate_output_'+self.casename+'.run', shell=True)
                    if self.verbose == True:
                        #print('Beginning 2 column Climate rerun')
                        ftestingoutput.write('Beginning 2 column Climate rerun\n')
                        ftestingoutput.close()
                        ftestingoutput = open(self.OutPath+self.casename+'_SavingInfoOut.txt', 'a')

                    self.run_climate_1instance(self.vplclimate_RunScriptDir+'RunVPLClimate_2column_'+self.casename+'.script', self.vplclimate_executable, trynum=self.num_2col_climate_runs, twocol=True)
                    if self.verbose == True:
                        #print('2 column Climate subtry number '+str(climate_subtries)+' completed')
                        ftestingoutput.write('2 column Climate subtry number '+str(climate_subtries)+' completed\n')

                    # Check convergence
                    local_climate_convergence, tropheating_dayside, tropheating_nightside, avgflux, cnvtype = self.check_2column_vplclimate_conv(trynum=self.num_2col_climate_runs, climsubtry=climate_subtries)

                    if self.verbose == True:
                        if local_climate_convergence == True:
                            #print('2 column Climate convergence found on subtry number '+str(climate_subtries)+' for run number '+str(self.num_2col_climate_runs))
                            ftestingoutput.write('2 column Climate convergence found on subtry number '+str(climate_subtries)+' for run number '+str(self.num_2col_climate_runs)+'\n')
                            ftestingoutput.write('2 col cnv type: '+str(cnvtype)+'\n')
                        else:
                            #print('2 column Climate convergence NOT found on subtry number '+str(climate_subtries)+' for run number '+str(self.num_2col_climate_runs)+', continuing rerun sequence')
                            ftestingoutput.write('2 column Climate convergence NOT found on subtry number '+str(climate_subtries)+' for run number '+str(self.num_2col_climate_runs)+', continuing rerun sequence\n')



                # Until climate converges, loop through taking new temp profile
                while local_climate_convergence == False:

                    if climate_subtries == self.max_iterations_climate:
                        break

                    climate_subtries += 1

                    # First, get the final profile from the last climate run to use in the restart
                    climate_profile = self.get_final_2column_climate_output_temp_profile(trynum=self.num_2col_climate_runs)

                    # Update temperature in the PT profile and update the surface temperature
                    self.replace_PT_tempcol(climate_profile['Dayside']['T[K]'], whichcolumn='dayside')
                    self.replace_PT_tempcol(climate_profile['Nightside']['T[K]'], whichcolumn='nightside')

                    # Recreate the runscript to update the surface temp
                    self.make_2column_climate_runscript(trynum=self.num_2col_climate_runs)

                    if self.num_2col_climate_runs == 1:
                        subprocess.run('mv '+self.OutPath+'vpl_2col_climate_output_'+self.casename+'.run '+self.OutPath+'vpl_2col_climate_output_'+self.casename+'_Subtry'+str(climate_subtries)+'.run', shell=True)
                    else:
                        subprocess.run('mv '+self.OutPath+'vpl_2col_climate_output_'+self.casename+'_Try'+str(self.num_2col_climate_runs)+'.run '+self.OutPath+'vpl_2col_climate_output_'+self.casename+'_Try'+str(self.num_2col_climate_runs)+'_Subtry'+str(climate_subtries)+'.run', shell=True)

                    # Re run Climate 
                    if self.verbose == True:
                        #print('Beginning 2 column Climate rerun')
                        ftestingoutput.write('Beginning 2 column Climate rerun\n')
                    self.run_climate_1instance(self.vplclimate_RunScriptDir+'RunVPLClimate_2column_'+self.casename+'.script', self.vplclimate_executable, trynum=self.num_2col_climate_runs, twocol=True)
                    if self.verbose == True:
                        #print('2 column Climate subtry number '+str(climate_subtries)+' completed')
                        ftestingoutput.write('2 column Climate subtry number '+str(climate_subtries)+' completed\n')

                    # Check convergence
                    local_climate_convergence, tropheating_dayside, tropheating_nightside, avgflux, cnvtype = self.check_2column_vplclimate_conv(trynum=self.num_2col_climate_runs, climsubtry=climate_subtries)

                    if self.verbose == True:
                        if local_climate_convergence == True:
                            #print('2 column Climate convergence found on subtry number '+str(climate_subtries)+' for run number '+str(self.num_2col_climate_runs))
                            ftestingoutput.write('2 column Climate convergence found on subtry number '+str(climate_subtries)+' for run number '+str(self.num_2col_climate_runs)+'\n')
                            ftestingoutput.write('2 col cnv type: '+str(cnvtype)+'\n')
                        else:
                            #print('2 column Climate convergence NOT found on subtry number '+str(climate_subtries)+' for run number '+str(self.num_2col_climate_runs)+', continuing rerun sequence')
                            ftestingoutput.write('2 column Climate convergence NOT found on subtry number '+str(climate_subtries)+' for run number '+str(self.num_2col_climate_runs)+', continuing rerun sequence\n')

                climate_profile = self.get_final_2column_climate_output_temp_profile(trynum=self.num_2col_climate_runs)

                # Update temperature in the PT profile and update the surface temperature
                self.replace_PT_tempcol(climate_profile['Dayside']['T[K]'], whichcolumn='dayside')
                self.replace_PT_tempcol(climate_profile['Nightside']['T[K]'], whichcolumn='nightside')

                if local_climate_convergence == False and self.suppress_IOerrors == False:
                    raise IOError('2 column Climate could not converge in >'+str(self.max_iterations_climate)+' tries. For climate run number '+str(self.num_2col_climate_runs))
                elif local_climate_convergence == False and self.suppress_IOerrors == True:
                    self.global_convergence = False

                if self.verbose == True:
                    ftestingoutput.close()
                

            ### Create SMART runscript ------------------------------

            if self.run_spectra == True and self.global_convergence == True and self.MultiNest_DataFit == False:

                if self.verbose == True:
                    ftestingoutput = open(self.OutPath+self.casename+'_SavingInfoOut.txt', 'a')

                self.make_smart_runscript()

                self.run_smart_1instance(self.SMART_RunScriptDir+'RunSMART_'+self.casename+'.run')

                if self.verbose == True:
                    #print('Terminator SMART run completed')
                    ftestingoutput.write('Terminator SMART run completed\n')
                    ftestingoutput.close()


                if self.rerun_smart_for_2col == True:

                    self.make_lblabc_runscripts(whichcol='dayside')
                    self.make_lblabc_runscripts(whichcol='nightside')
                    lblabc_input = []
                    for gas in self.molecule_dict['Gas_names']:
                        lblabc_input.append([self.lblabc_RunScriptDir+'RunLBLABC_d_'+gas+'_'+self.casename+'.script', gas, 'dayside'])
                        lblabc_input.append([self.lblabc_RunScriptDir+'RunLBLABC_n_'+gas+'_'+self.casename+'.script', gas, 'nightside'])
                    
                    #with Pool() as p:
                    #    lblruns = p.map(self.run_lblabc_1instance_Parallel, lblabc_input)
                    
                    for lblinp in lblabc_input:
                        self.run_lblabc_1instance_Parallel(lblinp)
                    
                    smartruninputs = ['dayside', 'nightside'] #'Avg', 

                    #with Pool() as p:
                        #smartrunsparallel = p.map(self.run_multinest_smart_parallel, smartruninputs)

                    for sminp in smartruninputs:
                        self.run_multinest_smart_parallel(sminp)
                    
            elif self.run_spectra == True and self.global_convergence == True and self.MultiNest_DataFit == True:

                shutil.copyfile(self.dayside_starting_PT, self.AtmProfPath+'PT_profile_dayside_'+self.casename+'.pt')
                shutil.copyfile(self.nightside_starting_PT, self.AtmProfPath+'PT_profile_nightside_'+self.casename+'.pt')

                self.make_lblabc_runscripts(whichcol='dayside')
                self.make_lblabc_runscripts(whichcol='nightside')
                lblabc_input = []
                for gas in self.molecule_dict['Gas_names']:
                    lblabc_input.append([self.lblabc_RunScriptDir+'RunLBLABC_d_'+gas+'_'+self.casename+'.script', gas, 'dayside'])
                    lblabc_input.append([self.lblabc_RunScriptDir+'RunLBLABC_n_'+gas+'_'+self.casename+'.script', gas, 'nightside'])
                
                #with Pool() as p:
                #    lblruns = p.map(self.run_lblabc_1instance_Parallel, lblabc_input)
                
                for lblinp in lblabc_input:
                    self.run_lblabc_1instance_Parallel(lblinp)
                
                smartruninputs = ['dayside', 'nightside'] #'Avg', 

                #with Pool() as p:
                    #smartrunsparallel = p.map(self.run_multinest_smart_parallel, smartruninputs)

                for sminp in smartruninputs:
                    self.run_multinest_smart_parallel(sminp)

        return self.global_convergence


