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

''' Delete after testing

# Input Molecule dict example:
# molecules_I_want = {'O2':7, 'O3':3}
# where 7 and 3 are the HITRAN gas codes for O2 and O3, respectively

def createmoldic():
    moldic = {}
    moldic['O2'] = 7
    moldic['O3'] = 3
    moldic['CO'] = 5
    moldic['CO2'] = 2
    moldic['H2O'] = 1
    moldic['HNO3'] = 12
    moldic['N2O'] = 4
    moldic['NO2'] = 10
    moldic['SO2'] = 9

    return moldic
'''

class VPLModelingPipeline:

    # Set Global and initialize atmosphere object:
    def __init__(self, casename, photochemInitial, verbose, find_molecules_of_interest=False, hitran_year='2020') -> None:
        # Set any and all needed paths
        self.photochemDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/' # path to PHOTOCHEM/ dir
        self.atmosDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/' # path to atmos/ dir
        self.lblabcDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/' # path to lblabc/ dir (such that lblabcDir/lblabc is the executable to call)
        self.OutPath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/'+casename+'/' # path for the raw model run outputs (NOT for created data products like dictionaries)
        self.DataOutPath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ModelRunOutputs/'+casename+'/' # path for created data products like dictionaries
        self.AtmProfPath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/AtmProfiles/' # path to put atmospheric profile files (.pt files really)
        self.photochemBackupDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/Save_Photochem_Output/'+casename+'/' # path to save output from each photochem run
        self.LBLABC_AbsFilesDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/LinebyLine_absFiles/'+casename+'/' # path to put the created lbl .abs files in 
        self.lblabc_RunScriptDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/LBLABC/'+casename+'/' # path to put lbl runscripts in
        self.vplclimate_RunScriptDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/RunFiles/VPLClimate/'+casename+'/' # path to put vpl climate runscripts in
        self.photochem_InputsDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/Bodies/'+casename+'/' # The path to create new photochem inputs in
        self.xsec_Path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/xsec/' # The path where cross section files can be found

        # The climate executable:
        self.vplclimate_executable = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimateModel/vpl_climate/vpl_climate' # The VPL Climate executable you want to use WITH FULL PATH

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
        self.nlevel_coarse = 60 # Number of atm layers for coarse grids (i.e., for all models besides photochem)

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
        #self.vplclimateInitial = vplclimateInitial # path to vpl climate template file
        # User needs to create a mostly filled in vpl climate template run script
        # Pipeline will automatically point to the right PT prof/mixing ratio files and iterate the MMW
        self.verbose = verbose # Boolean, whether or not you want print statements (False for computational efficiency)
        # MOLECULES MUST BE ALL CAPITAL LETTERS AS THEY WILL PRINT OUT FROM PHOTOCHEM

        # put photochem initial inputs into the inputs dir if they aren't there already
        if photochemInitial != self.photochem_InputsDir:
            subprocess.run('cp '+photochemInitial+'input_photchem.dat '+self.photochem_InputsDir, shell=True)
            subprocess.run('cp '+photochemInitial+'parameters.inc '+self.photochem_InputsDir, shell=True)
            #subprocess.run('cp '+photochemInitial+'params.dat '+self.photochem_InputsDir, shell=True)
            subprocess.run('cp '+photochemInitial+'PLANET.dat '+self.photochem_InputsDir, shell=True)
            subprocess.run('cp '+photochemInitial+'reactions.rx '+self.photochem_InputsDir, shell=True)
            subprocess.run('cp '+photochemInitial+'species.dat '+self.photochem_InputsDir, shell=True)
            subprocess.run('cp '+photochemInitial+'in.dist '+self.photochem_InputsDir, shell=True)


        # lookup table to connect Hitran gas codes to molecule names
        self.hitran_lookup = pd.read_csv('HitranTable.csv', index_col='Molecule')

        if find_molecules_of_interest == False:
            self.molecule_dict = {} # key-value pairs of molecules of interest (keys, str) and their hitran codes (value, int)
            gas_names = ['O2', 'O3', 'H2O']#, 'CO', 'CO2', 'HNO3', 'N2O', 'NO2', 'SO2']
            self.molecule_dict['Gas_names'] = gas_names
            for m in range(len(gas_names)):
                self.molecule_dict[gas_names[m]] = self.hitran_lookup.loc[gas_names[m]]['HitranNumber']
                self.molecule_dict[gas_names[m]+'_RmixCol'] = m+1




    ######################################################### Support Functions

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
        self.vplclimate_executable = 'something_Beta17?' # The VPL Climate executable you want to use WITH FULL PATH

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

        # Save surface temperature for use in climate 
        self.surface_temp = temp[0]

        new_grid = np.linspace(alt[0], alt[len(alt)-1], self.nlevel_coarse)
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
    def replace_PT_tempcol(self, new_temp):

        # Read in the old PT profile
        oldpt = ascii.read(self.AtmProfPath+'PT_profile_'+self.casename+'.pt')

        # replace the temp
        oldpt['Temp'] = new_temp

        # Update Surface Temperature
        self.surface_temp = new_temp[len(new_temp)-1]

        # Overwrite the PT profile
        ascii.write(oldpt, self.AtmProfPath+'PT_profile_'+self.casename+'.pt', overwrite=True)

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
                    Nstep_photochemrun = int(hold[2])
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
            #if HasItConverged:
            #    print('Photochem run '+self.casename+' Try '+str(trynum)+' has converged!')
            #else:
            #    print('Photochem run '+self.casename+' Try '+str(trynum)+' has NOT converged.')
            print('Normalized Gross error: '+str(NormGrosserr))
            print('L2 Error: '+str(L2err))
            print('Time of final timestep: '+str(FinalTime)+'\n')

        return HasItConverged, NormGrosserr, L2err, FinalTime, Nstep_photochemrun
    # usage should be 'convergence, grosserr, l2err, finaltime, nstepsphoto = check_photochem_conv()

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
        os.mkdir(self.photochemBackupDir+'RunNumber'+str(trynum)+'/')
        subprocess.run('cp '+self.photochemDir+'PTZ_mixingratios_in.dist '+self.photochemBackupDir+'RunNumber'+str(trynum)+'/', shell=True)
        subprocess.run('cp '+self.photochemDir+'OUTPUT/* '+self.photochemBackupDir+'RunNumber'+str(trynum)+'/', shell=True)
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
            f.write('50,100000                               min, max wavenumber\n')
            f.write('200.                                    maximum line width\n')
            f.write('1.e-5                                   minimum column optical depth\n')
            f.write(self.HITRAN_FundamentalFile+'\n')
            f.write(self.LBLABC_AbsFilesDir+i+'_'+self.casename+'.abs\n')
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

        planet = 'T1c' # Would suggest keeping all settings available for each target, this provides a quick way to switch between them

        if planet == 'T1c':
            self.c_NumberTimesteps = 10000
            self.c_TimestepLength = 86400.0 # [s]
            self.c_SubstepIntervals = 20
            self.c_TimestepOutputIntervals = 50
            self.c_TimeStepMethodIndex = 7
            self.c_TempChangeTolerance = 0.01
            self.c_DoubledRadiationGrid = True
            self.c_HRTCalcType = 3 # 3 - Global Hrt
            self.c_NumberSolarZeniths = 4 # number of solar zenith angles used in avg
            self.c_IncludeRadiativeHrt = True
            self.c_IncludeConvectiveHrt = True
            self.c_IncludeConductiveHrt = False
            self.c_DayLength = 209260.80000000002 # [s], 2.42 days
            self.c_YearLength = 1.0 # Length of year in days according to c_DayLength (1 when tidally locked)
            self.c_OrbitalCalcType = 0 # 0 - fixed orbital distance
            self.c_SemiMajorAxis = 0.0158 # [AU]

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

    ''' Delete after testing, or change to be another option for climate script generation

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
    '''

    ### Purpose: Change the T/EDD profiles from climate in photochems in.dist to rerun photochem after a climate run
    ##
    ## Attribute Dependencies:
    # photochem_InputsDir - to get NZ and NQ from the parameters.inc file 
    #
    ## Fxn-specific Inputs:
    # oldptz - path to the PTZ_out file from the last photochem run (use PHOTOCHEM/OUTPUT/ dir)
    # climateprof - the last output profile from vpl climate, found with self.get_final_climate_output_temp_profile()
    ##
    def update_indist(self, oldptz, climateprof):

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
        T_edd_block['col2'] = new_edd

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




    ### AUTOMATIC PIPELINE, calls all above functions sequentially defined by flow chart, updates necessary object parameters
    ##
    ## Dependent on all Attributes, no fxn-specific inputs
    ##
    def run_automatic(self):
        # Prepare your backup directory for photochem data
        if not os.path.exists(self.photochemBackupDir):
            os.mkdir(self.photochemBackupDir)
        else:
            subprocess.run('rm -rf '+self.photochemBackupDir+'*', shell=True)
        # Prepare your directory for storing lblabc .abs files
        if not os.path.exists(self.LBLABC_AbsFilesDir):
            os.mkdir(self.LBLABC_AbsFilesDir)
        else:
            subprocess.run('rm -rf '+self.LBLABC_AbsFilesDir+'*', shell=True)
        # Prepare directory for storing lblabc run script files
        if not os.path.exists(self.lblabc_RunScriptDir):
            os.mkdir(self.lblabc_RunScriptDir)
        # Prepare directory for storing vpl climate run script files
        if not os.path.exists(self.vplclimate_RunScriptDir):
            os.mkdir(self.vplclimate_RunScriptDir)
        # Prepare directory for storing new photochem inputs
        if not os.path.exists(self.photochem_InputsDir):
            os.mkdir(self.photochem_InputsDir)
        # Make sure model and data output dirs exist
        if not os.path.exists(self.OutPath):
            os.mkdir(self.OutPath)
        if not os.path.exists(self.DataOutPath):
            os.mkdir(self.DataOutPath)
        # Prepare the Hyak environment
        #self.prepare_hyak_env()

        # Start loop to find global convergence with photochem + lblabc + vpl climate
        while self.global_convergence == False:

            ### Run photochem section start ---------------------

            self.num_photochem_runs += 1
            if self.verbose == True:
                print('----> Beginning photochem run Try number '+str(self.num_photochem_runs))
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
            local_photochem_conv, grosserr, l2err, finaltime, nsteps_photo = self.check_photochem_conv()

            # If photochem did not converge, try try again
            while local_photochem_conv == False:
                subprocess.run('rm -rf '+self.photochemDir+'in.dist', shell=True)
                subprocess.run('rm -rf '+self.photochemDir+'PTZ_mixingratios_in.dist', shell=True)
                subprocess.run('cp '+self.photochemDir+'OUTPUT/out.dist '+self.photochemDir+'in.dist', shell=True)
                self.run_photochem_1instance(CleanMake=False, InputCopy=False, trynum=self.num_photochem_runs)
                photochem_subtries += 1
                local_photochem_conv, grosserr, l2err, finaltime, nsteps_photo = self.check_photochem_conv()

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
                if nsteps_photo < 1000:
                    self.global_convergence = True
                    if self.verbose == True:
                        print('Global Convergence achieved')
                    break
                else:
                    self.global_convergence = False

            ### Run photochem section end --------------------------

            ### Create degraded atmospheric profile to prepare for LBLABC and Climate -------------
            self.degrade_PT()
            if self.verbose == True:
                print('Degraded PT profile created from photochem run '+str(self.num_photochem_runs))
            ### Degraded PT Profile finished ------------

            ### Create new mixing ratios profile file --------------------------
            self.prep_rmix_file(self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist')
            ### Mixing ratios profile created --------------------------

            ''' Delete after testing, switched to hitran lookup table and creating mixingRs file

            # If LBLABC hasnt run yet, Get the columns of rmix for the molecules of interest, save molecule names though
            if self.num_lblabc_runs == 0:
                mixingRs = ascii.read(self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist')
                self.molecule_dict['Gas_names'] = list(self.molecule_dict.keys())
                for gas in self.molecule_dict['Gas_names']:
                    for i in range(len(mixingRs.keys())):
                        if mixingRs.keys()[i] == gas:
                            self.molecule_dict[gas+'_RmixCol'] = i+1
                            break
            '''

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
                    print('LBLABC run for '+gas+' complete, LBLABC iteration '+str(self.num_lblabc_runs+1))
            self.num_lblabc_runs += 1

            ### Run LBLABC section end -------------------------------------
            
            ### Run VPL Climate section start ------------------------------

            # Add 1 to the climate run counter
            self.num_climate_runs += 1

            # First Generate the climate runscript
            # Things, e.g., MMW need to be updated, so just do this every time
            self.make_climate_runscript(trynum=self.num_climate_runs)

            if self.verbose == True:
                print('Climate Runscript created, beginning first climate run')

            # Now run climate 
            self.run_climate_1instance(self.vplclimate_RunScriptDir+'RunVPLClimate_'+self.casename+'.script', self.vplclimate_executable, trynum=self.num_climate_runs)

            if self.verbose == True:
                print('First Climate run completed')

            # Check for local convergence of climate, similar to process followed for a given try on photochem            
            local_climate_convergence, tropheating, avgflux = self.check_vplclimate_conv(trynum=self.num_climate_runs)

            if self.verbose == True:
                if local_climate_convergence == True:
                    print('Climate convergence found on first try for run number '+str(self.num_climate_runs))
                else:
                    print('Climate convergence NOT found on first try for run number '+str(self.num_climate_runs)+', beginning rerun sequence')

            climate_subtries = 1
            # Until climate converges, loop through taking new temp profile
            while local_climate_convergence == False:

                climate_subtries += 1

                # First, get the final profile from the last climate run to use in the restart
                climate_profile = self.get_final_climate_output_temp_profile(trynum=self.num_climate_runs)

                # Update temperature in the PT profile and update the surface temperature
                self.replace_PT_tempcol(climate_profile['T[K]'])

                # Recreate the runscript to update the surface temp
                self.make_climate_runscript(trynum=self.num_climate_runs)

                # Re run Climate 
                if self.verbose == True:
                    print('Beginning Climate rerun')
                self.run_climate_1instance(self.vplclimate_RunScriptDir+'RunVPLClimate_'+self.casename+'.script', self.vplclimate_executable, trynum=self.num_climate_runs)
                if self.verbose == True:
                    print('Climate subtry number '+str(climate_subtries)+' completed')

                # Check convergence
                local_climate_convergence, tropheating, avgflux = self.check_vplclimate_conv(trynum=self.num_climate_runs)

                if self.verbose == True:
                    if local_climate_convergence == True:
                        print('Climate convergence found on subtry number '+str(climate_subtries)+' for run number '+str(self.num_climate_runs))
                    else:
                        print('Climate convergence NOT found on subtry number '+str(climate_subtries)+' for run number '+str(self.num_climate_runs)+', continuing rerun sequence')

            ### Run VPL Climate section end ------------------------------

            ## Last thing to do: update in.dist for next photochem run
            ### Update in.dist for next photochem run ------------------------------

            # If this is the first change, save the original in.dist just in case
            if self.num_photochem_runs == 1:
                subprocess.run('cp '+self.photochem_InputsDir+'in.dist '+self.photochem_InputsDir+'Original_in.dist', shell=True)
            
            # Now update in.dist
            climate_profile = self.get_final_climate_output_temp_profile(trynum=self.num_climate_runs)
            self.update_indist(self.photochemDir+'OUTPUT/PTZ_mixingratios_out.dist', climate_profile)

            ### Update in.dist section end ------------------------------


