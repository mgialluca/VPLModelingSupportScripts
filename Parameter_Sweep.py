# Units of fluxes are molecules / cm^2 / s

import sys
#sys.path.append("/gscratch/vsm/gialluca/anaconda3/lib/python3.9/site-packages")

from Pipeline import *
import numpy as np
import astropy.units as u
import astropy.constants as const 
import os
import copy
from multiprocessing import Pool
import shutil
import re
import json
import emcee
import h5py
from functools import partial
from numpy import log, exp, pi
import scipy.stats, scipy
import sys
from spectral_utils import *
from planck import *
import pandas as pd
from scipy.optimize import minimize

# To use pymultinest:
# - Load the intel module
# - export the python path back to anaconda
# - Load the ompi module
# - Add the multinest path to LD_LIBRARY_PATH

import pymultinest
import mpi4py
from mpi4py import MPI


# Need to figure out how parameter sweeps are running:
# 1. Run across a grid (every combination? Or strategic points? former is brute force method, could start with that)
# 2. Fit for a particular case with emcee 

# General Notes:
# - Could do a dict for grid; define the molecules to change as keys, allow user to define either the array of points to test at
###  or allow the user to specify a min/max range, log or linear sampling, and sampling resolution for np
# - Have user give input units of flux via astropy units, can do conversion to molec/cm2/s internally 

class Generate_Atmosphere_Parameter_Sweep:

    def __init__(self, sweepname, photochemInitial, restart_run=False, starting_point='Exact', hitran_year='2020', climate2col=True, spectra=True, planet='T1c'):

        self.sweepname = sweepname # Naming convention for directory structure
        self.photochemInitial = photochemInitial # Input files for photochem to copy and change 
        self.hitran_year = hitran_year # hitran year, 2016 or 2020 (default should be latter)
        self.supernode = True # If you are using the supernode on hyak, needs to be True
        if hitran_year == '2020':
            self.lblabc_qtxt_dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/hitranQtips2020/' # For the hitran distribution you want
        elif hitran_year == '2016':
            self.lblabc_qtxt_dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/hitranQtips/'

        self.planet = planet

        if self.planet == 'T1b':
            self.R_p = 1.116*u.Rearth
        elif self.planet == 'T1c':
            self.R_p = 1.097*u.Rearth # Currently radius of Trappist-1c 
        elif self.planet == 'T1d':
            self.R_p = 0.788*u.Rearth
        elif self.planet == 'T1e':
            self.R_p = 0.920*u.Rearth
        elif self.planet == 'T1f':
            self.R_p = 1.045*u.Rearth
        elif self.planet == 'T1g':
            self.R_p = 1.129*u.Rearth
        elif self.planet == 'T1h':
            self.R_p = 0.755*u.Rearth

        self.atmos_Dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/' # Path to dir containing atmos, will be copied for runs 

        self.mcmc_pressure_only = False
        self.multinest_fit_data = False
        self.multinest_input_options_photochem = restart_run # Gives the input options for photochem files for each simulation in a multinest suite
        self.multinest_input_options_climate = 'Cinit' # Gives the input options for climate profiles which may have differing species, indist etc 
        self.multinest_use_12um = False

        #########  Parameters to set if you want to do a grid sweep  #########

        self.outgass_species_gridsweep = ['H2O', 'CO2', 'SO2'] # Species to vary outgassing rates of
        self.outgass_species_sourcetype = ['Flx', 'Flx', 'FixMR'] # Flux at surface (Flx) or Fixed mixing ratio (FixMR)
        self.outgass_species_molarmass = {} # Molar masses in g/mol
        self.outgass_species_molarmass['H2O'] = [18.015]*(u.g/u.mol)# Molar masses in g/mol
        self.outgass_species_molarmass['CO2'] = [44.01]*(u.g/u.mol)
        self.outgass_species_molarmass['SO2'] = [64.066]*(u.g/u.mol)
        self.escape_species_gridsweep = ['O', 'O2', 'O3', 'H2O2', 'CO2', 'CO'] # Species to vary escape rates of
        self.escape_species_losstype = ['Veff', 'Veff', 'Vdep', 'Vdep', 'Veff', 'Vdep'] # Vdep (depositional velocity at surface) or TOA (flux at top of atmosphere)
        self.escape_species_molarmass = {}
        self.escape_species_molarmass['O'] = [15.999]*(u.g/u.mol) 
        self.escape_species_molarmass['O2'] = [31.998]*(u.g/u.mol) 
        self.escape_species_molarmass['O3'] = [47.997]*(u.g/u.mol) 
        self.escape_species_molarmass['H2O2'] = [34.014]*(u.g/u.mol) 
        self.escape_species_molarmass['CO2'] = [44.01]*(u.g/u.mol) 
        self.escape_species_molarmass['CO'] = [28.01]*(u.g/u.mol) 

        self.outgass_sample_type_gridsweep = ['Log', 'UserDef', 'UserDef'] # How to sample outgassed molecules: 'Linear', 'Log', or 'UserDef' 
        self.escape_sample_type_gridsweep = ['UserDef', 'UserDef', 'UserDef', 'UserDef', 'UserDef', 'UserDef'] #['UserDef', 'UserDef', 'UserDef', 'UserDef']
        # Linear - sample every flux on a linear grid (np.linspace) with some defined resolution
        # Log - sample every flux on a log grid (np.logspace) with some defined resolution
        # UserDef - User defined arrays of samples for every flux to vary 

        # Need to set Min / Max ranges for each molecule to vary in the form of a dictionary if using Linear or Log sampling
        self.outgass_species_MinMax_gridsweep = {}
        self.outgass_species_MinMax_gridsweep['H2O'] = [44552887.2545331, 9.47899801e11]
        self.outgass_species_MinMax_gridsweep['CO2'] = []
        self.outgass_species_MinMax_gridsweep['SO2'] = []
        #[1.00028455e+11, 1.00089313e+11]
        #[9.97550516e+10, 9.99980399e+10]
        #[9.96337789e+10, 1.00608103e+11]
        #[8.48538802e+10, 9.91501611e+10] 
        #[1.07177663e+11, 1.99798009e+11]#[2.25884849e+10, 2.72793861e+11]# full bound: [44552887.2545331, 9.47899801e11] 
        #[90000000000.0, 100000000000.0] #[1.65329797e8, 3.00359578e12] # min, max

        self.escape_species_MinMax_gridsweep = {}
        self.escape_species_MinMax_gridsweep['O'] = []
        self.escape_species_MinMax_gridsweep['O2'] = []
        self.escape_species_MinMax_gridsweep['O3'] = []
        self.escape_species_MinMax_gridsweep['H2O2'] = []
        self.escape_species_MinMax_gridsweep['CO2'] = []
        self.escape_species_MinMax_gridsweep['CO'] = []
        

        # Sample resolution if using Linear / Log sampling 
        self.outgass_sample_resolution_gridsweep = [8, 0, 0] # number of samples for each outgassed species
        self.escape_sample_resolution_gridsweep = []

        # Need to pass samples for user defined option
        self.outgass_samples_gridsweep = {}
        self.outgass_samples_gridsweep['CO2'] = [0]
        self.outgass_samples_gridsweep['SO2'] = [0.0001] #100 ppm

        self.escape_samples_gridsweep = {}
        self.escape_samples_gridsweep['O'] = [0.1, 0.5]#[0, 1e27, 1e29]#[0, 1e27, 1e28, 1e29] #[1e28, 1e29, 1e30] #[0, 1e26, 1e27] #[0, 1e23, 5e23, 1e24, 5e24, 1e25, 5e25, 1e26]
        self.escape_samples_gridsweep['O2'] = [0.05, 0.1]#[1e26, 5e26, 1e27]
        self.escape_samples_gridsweep['O3'] = [0.02, 0.4] 
        self.escape_samples_gridsweep['H2O2'] = [0.02]
        self.escape_samples_gridsweep['CO2'] = [0.02]
        self.escape_samples_gridsweep['CO'] = [0.001]
        


        # Units for either Min/Max values, or the user defined samples 
        self.outgass_species_units_gridsweep = 1 / (u.cm**2 * u.s) # molecules / cm2*s (can convert from mass/time with molar mass or mol/time)
        #self.escape_species_units_gridsweep = 1 / u.s # Molecules per second
        self.escape_species_units_gridsweep = [u.cm/u.s, u.cm/u.s, u.cm/u.s, u.cm/u.s, u.cm/u.s, u.cm/u.s] # Molecules per second

        #######################################################################

        #########  Setting output paths  #########

        self.Restart_Run = restart_run
        self.Starting_Point = starting_point # If restart_run is a string, it gives the initial files to use separately ...
        # ... if this is 'Exact', that means the runs inputs will be the exact same ...
        # ... otherwise this will point to a run statistics file to use to determine the closest available input files
        self.master_out = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+self.sweepname+'/'
        self.include_clim2col = climate2col
        self.include_spectra = spectra
        self.adjust_N2 = False # Set to the partial pressure number if you want to change that in the runs
        

        # UNCOMMENT BEFORE RUNNING:
        if not os.path.exists(self.master_out):
            os.makedirs(self.master_out, exist_ok=True)
        elif self.Restart_Run == True:
            subprocess.run('rm -rf '+self.master_out+'/*', shell=True)

        ##########################################


    # Set relevant variables of a model pipeline object that is constant for all runs in the sweep
    # casename - specific for each pipeline initialization
    # pipelineobj - initialized pipeline object
    def set_pipeline_vars(self, casename, pipelineobj):

        # Paths are the main thing to set, because they will be massive amounts of running/files, want to keep each sweep colocated in one master dir

        # Create casename dir
        if not os.path.exists(self.master_out+casename+'/'):
            os.mkdir(self.master_out+casename+'/')
        
        #subprocess.run('cp -r '+self.atmos_Dir+' '+self.master_out+casename+'/atmos/', shell=True)
        shutil.copytree(self.atmos_Dir, self.master_out + casename + '/atmos/')

        pipelineobj.photochemDir = self.master_out+casename+'/atmos/PHOTOCHEM/' # path to PHOTOCHEM/ dir
        pipelineobj.atmosDir = self.master_out+casename+'/atmos/' # path to atmos/ dir
        pipelineobj.lblabcDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/' # path to lblabc/ dir (such that lblabcDir/lblabc is the executable to call)
        pipelineobj.OutPath = self.master_out+casename+'/' # path for the raw model run outputs (NOT for created data products like dictionaries)
        pipelineobj.DataOutPath = self.master_out+casename+'/' # path for created data products like dictionaries
        pipelineobj.AtmProfPath = self.master_out+casename+'/' # path to put atmospheric profile files (.pt files really)
        pipelineobj.BackupPhotochemRuns = False # Make backups of individual photochem runs
        pipelineobj.photochemBackupDir = self.master_out+casename+'/PhotochemBackup/' # path to save output from each photochem run
        pipelineobj.LBLABC_AbsFilesDir = self.master_out+casename+'/ABSFiles/' # path to put the created lbl .abs files in 
        pipelineobj.lblabc_RunScriptDir = self.master_out+casename+'/' # path to put lbl runscripts in
        pipelineobj.vplclimate_RunScriptDir = self.master_out+casename+'/' # path to put vpl climate runscripts in
        pipelineobj.photochem_InputsDir = self.master_out+casename+'/PhotochemInputs/' # The path to create new photochem inputs in
        pipelineobj.xsec_Path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/xsec/' # The path where cross section files can be found
        pipelineobj.SMART_RunScriptDir = self.master_out+casename+'/'

        # Adjust the atmospheric pressure
        pipelineobj.adjust_atmospheric_pressure = True
        pipelineobj.suppress_IOerrors = True
        pipelineobj.MCMC_pressure_only = self.mcmc_pressure_only
        pipelineobj.MultiNest_DataFit = self.multinest_fit_data
        pipelineobj.include_2column_climate = self.include_clim2col
        pipelineobj.rerun_smart_for_2col = self.include_clim2col
        pipelineobj.run_spectra = self.include_spectra
        pipelineobj.c_NumberSolarZeniths = 1

        if type(self.adjust_N2) != bool:
            pipelineobj.adjust_N2_amount = True
            pipelineobj.N2_fixed_pressure = self.adjust_N2
            pipelineobj.NewPressure_Psurf_tolerance = 0.06
        
        if self.mcmc_pressure_only == True:
            pipelineobj.include_2column_climate = False
            pipelineobj.run_spectra = False            

        if self.multinest_fit_data == True:

            pipelineobj.run_spectra = True
            pipelineobj.multinest_climcopy_dbOptions = 'Cinit'

            # THIS will now be done within the pipeline, after photochem has converged to truly select the correct profile
            #
            # copycase = pipelineobj.multinest_climate_copycase.split('/')[1]
            # pipelineobj.dayside_starting_PT = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+pipelineobj.multinest_climate_copycase+'/PT_profile_dayside_'+copycase+'.pt'
            # pipelineobj.nightside_starting_PT = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+pipelineobj.multinest_climate_copycase+'/PT_profile_nightside_'+copycase+'.pt'

            # ## Set the day/night surface temps
            # climoutcopy = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+pipelineobj.multinest_climate_copycase+'/vpl_2col_climate_output_'+copycase+'.run'
            # fi = open(climoutcopy, 'r')
            # lines = fi.readlines()
            # fi.close()

            # # Want to get the last output trop heating rate and avg flux, should be last two lines
            # # so loop in reversed order, break loop after to conserve efficiency

            # nightside_found = False
            # for i in reversed(range(len(lines))):
            #     hold = lines[i].split()
            #     if len(hold) > 2:
            #         if hold[0] == 'surface:':
            #             if nightside_found == False:
            #                 pipelineobj.surface_temp_nightside = float(hold[8])
            #                 nightside_found = True
            #             else:
            #                 pipelineobj.surface_temp_dayside = float(hold[8])
            #                 # After retrieving surface temp for nightside, will have all values, break
            #                 break

            

        # Testing if climate executable needs to be copied
        if self.supernode == True:
            pipelineobj.vplclimate_executable = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimateModel/vpl_climate_supernode'

        # Molecules for the type of atmosphere we're interested in 

        pipelineobj.molecule_dict = {} # key-value pairs of molecules of interest (keys, str) and their hitran codes (value, int)
        gas_names = ['O2', 'H2O', 'O3', 'CO2', 'CO', 'SO2']
        pipelineobj.molecule_dict['Gas_names'] = gas_names
        for m in range(len(gas_names)):
            pipelineobj.molecule_dict[gas_names[m]] = pipelineobj.hitran_lookup.loc[gas_names[m]]['HitranNumber']
            pipelineobj.molecule_dict[gas_names[m]+'_RmixCol'] = m+2


    # Prepare the hyak environment with ifort, python, and the HITRAN you want to use
    ## lblabc_qtxt_dir - the qtxt dir for the hitran distribution you want to use with lblabc
    def prepare_hyak_env(self):
        subprocess.run('module load intel/oneAPI/2021.1.1', shell=True)
        subprocess.run("alias python='/gscratch/vsm/gialluca/anaconda3/bin/python'", shell=True)
        subprocess.run('export LBLABC_QTXT_DIR='+self.lblabc_qtxt_dir, shell=True)


    # Change any flux unit to molecules / cm**2*s 
    ## var - variable to convert
    ## species - the species the variable corresponds to (to use molar mass properly)
    ## fluxtype - 'escape' or 'outgass' 
    def fix_flux_units(self, var, species, fluxtype, loss_type='TOA', source_type='Flx'):

        if (fluxtype == 'escape' and loss_type == 'TOA') or (fluxtype == 'outgass' and source_type == 'Flx'):
            try:
                var = var.to(1 / (u.cm**2 * u.s))
            except:
                try:
                    var = var.to(u.mol / (u.cm**2 * u.s))
                    var = (var*const.N_A).to(1 / (u.cm**2 * u.s))
                except:
                    try:
                        var = var.to(u.mol / u.s)
                        var = var / (4*np.pi*(self.R_p**2))
                        var = var.to(u.mol / (u.cm**2 * u.s))
                        var = (var*const.N_A).to(1 / (u.cm**2 * u.s))
                    except:
                        try:
                            var = var.to(u.kg / u.s)
                            if fluxtype == 'outgass':
                                var = (var / self.outgass_species_molarmass[species]).to(u.mol/u.s)
                            elif fluxtype == 'escape':
                                var = (var / self.escape_species_molarmass[species]).to(u.mol/u.s)
                            var = var / (4*np.pi*(self.R_p**2))
                            var = var.to(u.mol / (u.cm**2 * u.s))
                            var = (var*const.N_A).to(1 / (u.cm**2 * u.s))
                        except:
                            try:
                                var = var.to(1/u.s)
                                var = var / (4*np.pi*(self.R_p**2))
                                var = var.to(1 / (u.cm**2 * u.s))
                            except:
                                raise IOError('Could not properly convert flux units to molecules/cm^2/s for '+species+' ('+fluxtype+')')
            
        elif (fluxtype == 'escape' and loss_type == 'Vdep') or (fluxtype == 'escape' and loss_type == 'Veff'):
            try:
                var = var.to(u.cm/u.s)
            except:
                raise IOError('Could not properly convert depositional or effusion velocity units to cm/s for '+species+' ('+fluxtype+')')
            
        elif fluxtype == 'outgass' and source_type == 'FixMR':
            pass
        
        return var
    

    # Replace input fluxes in the species.dat file
    ## pipelineobj - the initialized pipeline object for the model with the variables correctly set by self.set_pipeline_vars()
    ## fluxes - the input fluxes to use, same as InputFlux for self.run_one_model but without the number to indicate the unique casename
    def replace_fluxes(self, pipelineobj, fluxes):

        # First need to create the inputs directory and copy the initial master files to change
        # Pipeline already has a function that does this:
        pipelineobj.setup_intial_photochem_dir()

        # Set pipeline objects "photochemInitialInput" dir to the inputs dir to avoid this happening again in the auto run
        pipelineobj.photochemInitialInput = pipelineobj.photochem_InputsDir

        # Set fluxes in the species.dat file
        # Fluxes are in order of outgassed species followed by escaping species 
        all_affected_species = self.outgass_species_gridsweep + self.escape_species_gridsweep
        nsp = open(pipelineobj.photochem_InputsDir+'species.dat', 'r')
        nsp_new = open(pipelineobj.photochem_InputsDir+'species_new.dat', 'w')
        inLLbloc = False
        LLblocdone = False
        for l in nsp.readlines():
            if l.split()[0][0] == '*':
                if inLLbloc == True: # If you're done with the LL species block can just write rest of file out
                    LLblocdone = True
                nsp_new.write(l)

            elif LLblocdone == False: # In the LL species block
                currgas = l.split()[0]
                inLLbloc = True
                if currgas in all_affected_species: # If the current gas needs to change outgassing and/or escape flux
                
                    currline = l.split()
                    nsp_new.write(currline[0])
                    add_spaces = 11-len(currline[0])
                    for space in range(add_spaces):
                        nsp_new.write(' ')
                    nsp_new.write(currline[1]+'  ') # will be 'LL' and then 2 spaces

                    # Now writing the 'O H C S N CL' block, each has a space after with 4 spaces after CL to get to LBOUND
                    nsp_new.write(currline[2]+' '+currline[3]+' '+currline[4]+' '+currline[5]+' '+currline[6]+' '+currline[7]+'    ')

                    # If the outgassing flux is changing, we can explicitly set this now at LBOUND
                    if currgas in self.outgass_species_gridsweep:

                        speciesind = self.outgass_species_gridsweep.index(currgas)
                        source_type = self.outgass_species_sourcetype[speciesind]

                        if source_type == 'Flx':

                            # get the index of the new value in fluxes array
                            fluxind = np.where(np.array(all_affected_species) == currgas)[0][0] # outgassing would always come first so you can use ind of 0 

                            # since we're outgassing, LBOUND will be '2' and VDEP0 and FIXEDMR will be '0.' with fixed amount of spaces
                            nsp_new.write('2     0.      0.      ')

                            # Now SGFLUX is set to be the new flux value 
                            newsgval = "{:.4E}".format(fluxes[fluxind])
                            nsp_new.write(newsgval)
                            add_spaces = 1#10-len(newsgval)
                            for space in range(add_spaces):
                                nsp_new.write(' ')
                            
                            # Now DISTH is set to '0.' with fixed spaces
                            nsp_new.write('0.      ')
                        
                        elif source_type == 'FixMR':

                            # get the index of the new value in fluxes array
                            fluxind = np.where(np.array(all_affected_species) == currgas)[0][0] # outgassing would always come first so you can use ind of 0 

                            # since we're outgassing, LBOUND will be '2' and VDEP0 will be '0.' with fixed amount of spaces
                            nsp_new.write('1     0.      ')

                            # Now FIXEDMR is set to be the new fixed mixing ratio value 
                            newsgval = "{:.4E}".format(fluxes[fluxind])
                            nsp_new.write(newsgval)
                            add_spaces = 1#10-len(newsgval)
                            for space in range(add_spaces):
                                nsp_new.write(' ')
                            
                            # Now SGFLUX and DISTH are set to '0.' with fixed spaces
                            nsp_new.write('0.      0.      ')


                    else: 

                        speciesind = self.escape_species_gridsweep.index(currgas)
                        loss_type = self.escape_species_losstype[speciesind]

                        if loss_type == 'Vdep' or loss_type == 'vdep': # Need the LBOUND to be a depositional velocity loss of the species
                            # Get the index of the new value in fluxes array
                            fluxind_hold = np.where(np.array(all_affected_species) == currgas)[0] # escape will always come last so need to be length agnostic
                            fluxind = fluxind_hold[len(fluxind_hold)-1]

                            # since we're using a depositional velocity, LBOUND will be '0' and VDEP0 and FIXEDMR will be '0.' with fixed amount of spaces
                            nsp_new.write('0     ')

                            # Now SGFLUX is set to be the new flux value 
                            newvdval = "{:.3E}".format(fluxes[fluxind])
                            nsp_new.write(newvdval+' ')
                            
                            # Now DISTH is set to '0.' with fixed spaces
                            nsp_new.write('0.      0.        0.      ')

                        else: # Need to preserve the lines LBOUND through DISTH (line indicies 8-12)
                            # First is the LBOUND with 5 spaces 
                            nsp_new.write(currline[8]+'     ')

                            # Then VDEP0 with 8 characters total
                            nsp_new.write(currline[9])
                            add_spaces = 8 - len(currline[9])
                            for space in range(add_spaces):
                                nsp_new.write(' ')

                            # Then FIXEDMR with 8 characters total
                            nsp_new.write(currline[10])
                            add_spaces = 8 - len(currline[10])
                            for space in range(add_spaces):
                                nsp_new.write(' ')

                            # Then SGFLUX with 10 characters total
                            nsp_new.write(currline[11])
                            add_spaces = 10 - len(currline[11])
                            for space in range(add_spaces):
                                nsp_new.write(' ')
                            
                            # Then DITSH with 8 characters total
                            nsp_new.write(currline[12])
                            add_spaces = 8 - len(currline[12])
                            for space in range(add_spaces):
                                nsp_new.write(' ')

                    # if the Escape flux is changing, we can now explicitly set that at MBOUND
                    if currgas in self.escape_species_gridsweep:

                        speciesind = self.escape_species_gridsweep.index(currgas)
                        loss_type = self.escape_species_losstype[speciesind]

                        if loss_type == 'TOA' or loss_type == 'toa':

                            # Get the index of the new value in fluxes array
                            fluxind_hold = np.where(np.array(all_affected_species) == currgas)[0] # escape will always come last so need to be length agnostic
                            fluxind = fluxind_hold[len(fluxind_hold)-1]

                            # MBOUND explicitly set to 2
                            nsp_new.write('2      ')

                            # Set the new flux as SMFLUX, always put one space after but NOTE if the exponent is >=100 there could be a character error in fortran file reading
                            newsmval = "{:.4E}".format(fluxes[fluxind])
                            nsp_new.write(newsmval+' ')

                            # VEFF0 set to 0. with fixed spaces, and then you're done
                            nsp_new.write('0.   ')
                            nsp_new.write('\n')

                        elif loss_type == 'Veff' or loss_type == 'veff':

                            # Get the index of the new value in fluxes array
                            fluxind_hold = np.where(np.array(all_affected_species) == currgas)[0] # escape will always come last so need to be length agnostic
                            fluxind = fluxind_hold[len(fluxind_hold)-1]

                            # MBOUND explicitly set to 0 and the SMFLUX to 0 
                            nsp_new.write('0      0.      ')

                            # Set the effusion velocity
                            newveffval = "{:.3E}".format(fluxes[fluxind])
                            nsp_new.write(newveffval+' ')
                            nsp_new.write('\n')

                        else: # Need to preserve the lines MBOUND through VEFF0 (line indicies 13-15)

                            # First MBOUND with 6 spaces
                            nsp_new.write(currline[13]+'      ')

                            # Then SMFLUX with 8 characters
                            nsp_new.write(currline[14])
                            add_spaces = 8 - len(currline[14])
                            if add_spaces < 0:
                                add_spaces == 0
                            for space in range(add_spaces):
                                nsp_new.write(' ')

                            # Finally VEFF0 with 5 characters
                            nsp_new.write(currline[15])
                            add_spaces = 8 - len(currline[15])
                            if add_spaces < 0:
                                add_spaces == 0
                            for space in range(add_spaces):
                                nsp_new.write(' ')
                            nsp_new.write('\n')

                    else: # Need to preserve the lines MBOUND through VEFF0 (line indicies 13-15)

                        # First MBOUND with 6 spaces
                        nsp_new.write(currline[13]+'      ')

                        # Then SMFLUX with 8 characters
                        nsp_new.write(currline[14])
                        add_spaces = 8 - len(currline[14])
                        if add_spaces < 0:
                            add_spaces == 0
                        for space in range(add_spaces):
                            nsp_new.write(' ')

                        # Finally VEFF0 with 5 characters
                        nsp_new.write(currline[15])
                        add_spaces = 8 - len(currline[15])
                        if add_spaces < 0:
                            add_spaces == 0
                        for space in range(add_spaces):
                            nsp_new.write(' ')
                        nsp_new.write('\n')

                else: # currgas not an affected species
                    nsp_new.write(l)    
            else: # Out of the LL species block
                nsp_new.write(l)

        nsp_new.close()
        nsp.close()

        # Delete old species and rename fixed version to be species.dat
        subprocess.run('rm '+pipelineobj.photochem_InputsDir+'species.dat', shell=True)
        subprocess.run('mv '+pipelineobj.photochem_InputsDir+'species_new.dat '+pipelineobj.photochem_InputsDir+'species.dat', shell=True)

    # To be used after 'replace_fluxes', basically if SO2 or CO2 is being added in, certain daughter products containing S and/or C ...
    # ... may have leftover deposition velocities set that need to be changed to 0 to allow daughter products to build up appropriately
    def change_deposition_daughter_species(self, pipelineobj):

        nsp = open(pipelineobj.photochem_InputsDir+'species.dat', 'r')
        nsp_new = open(pipelineobj.photochem_InputsDir+'species_new.dat', 'w')

        lines = nsp.readlines()
        for l in lines:
            hold = l.split()
            if len(hold) > 0:
                if 'S' in hold[0] or 'C' in hold[0]: 
                    # You have a S bearing species that isn't SO2 or a C species that isn't CO2 or C
                    if hold[0] not in ['SO2', 'CO2', 'CO']:
                        nsp_new.write(hold[0])
                        add_spaces = 11-len(hold[0])
                        for space in range(add_spaces):
                            nsp_new.write(' ')
                        nsp_new.write(hold[1]+'  ') # will be 'LL' and then 2 spaces
                        # Now writing the 'O H C S N CL' block, each has a space after with 4 spaces after CL to get to LBOUND
                        nsp_new.write(hold[2]+' '+hold[3]+' '+hold[4]+' '+hold[5]+' '+hold[6]+' '+hold[7]+'    ')
                        nsp_new.write('0     0.      0.      0.        0.      0      0.      0.     \n') # zeros across

                    else:
                        nsp_new.write(l)

                else: # Anything else stays the same
                    nsp_new.write(l)

            else:
                nsp_new.write(l)

        nsp_new.close()
        nsp.close()

        # Delete old species and rename fixed version to be species.dat
        subprocess.run('rm '+pipelineobj.photochem_InputsDir+'species.dat', shell=True)
        subprocess.run('mv '+pipelineobj.photochem_InputsDir+'species_new.dat '+pipelineobj.photochem_InputsDir+'species.dat', shell=True)

    # Euclidean distance metric to find out which previous run is closest to the current one
    # Changed to a percent change 
    def euclidean_distance(self, a, b):
        return np.sqrt(np.sum((np.log1p(a) - np.log1p(b)) ** 2))
    
    # Find the closest model from a previous sweep to use as a starting point
    # Compare to the newfluxes for the current model that overlap with the previous sweep
    def find_closest_prev_model(self, input_options, newfluxes):

        # Get the outgassed and escaped species from the previous run
        prev_gases = list(input_options.columns)
        prev_outgassed = []
        prev_escaped = []
        for prev in prev_gases:
            if len(prev.split('_Outgass')) == 2:
                prev_outgassed.append(prev.split('_Outgass')[0])
            elif len(prev.split('_Escape')) == 2:
                prev_escaped.append(prev.split('_Escape')[0])

        # Find which gasses will overlap between the previous run and the current run
        overlap_outgass = list(np.intersect1d(np.array(self.outgass_species_gridsweep), np.array(prev_outgassed)))
        overlap_escape = list(np.intersect1d(np.array(self.escape_species_gridsweep), np.array(prev_escaped)))

        # Extract the fluxes need to be matched for this current model, in the correct order
        compare_to_curr_run = np.zeros(len(overlap_outgass)+len(overlap_escape))
        fluxind = 0
        for currgas in self.outgass_species_gridsweep:
            if currgas in overlap_outgass:
                compareind = overlap_outgass.index(currgas)
                compare_to_curr_run[compareind] = newfluxes[fluxind]
            fluxind += 1
        
        for currgas in self.escape_species_gridsweep:
            if currgas in overlap_escape:
                compareind = overlap_escape.index(currgas) + len(overlap_outgass)
                compare_to_curr_run[compareind] = newfluxes[fluxind]
            fluxind += 1

        # Find closest model from previous sweep
        closestmodel = None
        closestdist = np.inf

        for i in input_options.index:
            # Get overlapping fluxes from previous sweep in correct order
            compare_prev_run = []
            for currgas in overlap_outgass:
                compare_prev_run.append(input_options.loc[i][currgas+'_OutgassRate'])
            for currgas in overlap_escape:
                compare_prev_run.append(input_options.loc[i][currgas+'_EscapeRate'])

            # Find distance of model
            if self.Starting_Point == 'Euclidean':
                curr_dist = self.euclidean_distance(compare_to_curr_run, np.array(compare_prev_run))

            # Determine if it's the closest model
            if curr_dist < closestdist:
                closestdist = curr_dist
                closestmodel = i 

        #print('Closest Model: '+closestmodel)

        return closestmodel
    
    # If you are running multinest and the climate copycase is different than the photochemical one
    # will need to change the T and EDD columns of in.dist
    def update_indist_T_EDD_Profiles(self, pipelineobj, climatecopycase):

        # WARNING: this may be hard coded if folks don't write their parameters.inc file the same way as templates etc
        # I don't think it will be a problem but just noting in case
        # This retrieves the NQ and NZ from parameters.inc
        #
        # This is for the current path to the photochemical copy case (will have MORE NQ than the climate one)
        NQFile = open(pipelineobj.photochem_InputsDir+'parameters.inc', 'r')
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

        ## Now do the same thing with the climate copy case to retrieve the T/EDD columns
        climateprofile_path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+climatecopycase+'/PhotochemInputs/'
        NQFile_clim = open(climateprofile_path+'parameters.inc', 'r')
        NQFileLines_clim = NQFile_clim.readlines()
        for l in NQFileLines_clim:
            if len(l.split('NQ=')) > 1:
                NQ_clim = l.split('NQ=')[1]
                NQ_clim = int(NQ_clim.split(',')[0])
                NZ_clim = l.split('NZ=')[1]
                NZ_clim = int(NZ_clim.split(',')[0])
                break

        # Find number of blocks of mixing ratios until T/EDD columns
        NQblocks_clim = np.ceil(NQ_clim/10)

        ###################
        # Now find the T and EDD profile we want to use from the climatecopycase
        # Keep track of the lines that are starting and ending the current block in the in.dist file
        blockstart_clim = 0
        blockend_clim = NZ_clim

        # this should handle all the mixing ratio blocks
        for i in range(int(NQblocks_clim)):
            blockstart_clim = blockend_clim
            blockend_clim = blockend_clim+NZ_clim

        # now the T/EDD/DEN/O3/CO2 block
        T_edd_block_clim = ascii.read(climateprofile_path+'in.dist', data_start=blockstart_clim, data_end=blockend_clim)
        blockstart_clim = blockend_clim
        new_T = T_edd_block_clim['col1']
        new_edd = T_edd_block_clim['col2']

        #########################
        ### Now to create new in.dist

        fnew = open(pipelineobj.photochem_InputsDir+'NEWin.dist', 'w')

        # Keep track of the lines that are starting and ending the current block in the in.dist file
        blockstart = 0
        blockend = NZ

        # this should handle all the mixing ratio blocks
        for i in range(int(NQblocks)):
            curr_nq_block = ascii.read(pipelineobj.photochem_InputsDir+'in.dist', data_start=blockstart, data_end=blockend)
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
        T_edd_block = ascii.read(pipelineobj.photochem_InputsDir+'in.dist', data_start=blockstart, data_end=blockend)
        blockstart = blockend
        T_edd_block['col1'] = new_T
        T_edd_block['col2'] = new_edd # Need to convert from m**2/ to cm**2/s; climate model will be between 1/2 and 1000, photochem will be 1e4-1e6 ish

        for line in range(len(T_edd_block)):
            fnew.write('   ')
            for col in T_edd_block.columns:
                fnew.write("{:.8E}".format(T_edd_block[col][line])+'   ')
            fnew.write('\n')

        # now the rest of in.dist
        olddist_txt = open(pipelineobj.photochem_InputsDir+'in.dist', 'r')
        alllines = olddist_txt.readlines()
        olddist_txt.close()
        for line in range(len(alllines)):
            if line >= blockstart:
                fnew.write(alllines[line])

        fnew.close()

        subprocess.run('rm '+pipelineobj.photochem_InputsDir+'in.dist', shell=True)
        subprocess.run('mv '+pipelineobj.photochem_InputsDir+'NEWin.dist '+pipelineobj.photochem_InputsDir+'in.dist', shell=True)


    # expand the in dist for more species
    def expand_indist(self, pipelineobj, climatecopycase, NZ=200):

        # Read in old species.dat file
        climateprofile_path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+climatecopycase+'/PhotochemInputs/'
        nsp = open(climateprofile_path+'species.dat', 'r')

        done = False
        oldgases = []
        for l in nsp.readlines():
            if l.split()[0][0] == '*':
                if done == True:
                    break
            else:
                oldgases.append(l.split()[0])
                done = True

        # Read in new species.dat file
        nsp = open(pipelineobj.photochem_InputsDir+'species.dat', 'r')

        done = False
        newgases = []
        for l in nsp.readlines():
            if l.split()[0][0] == '*':
                if done == True:
                    break
            else:
                newgases.append(l.split()[0])
                done = True

        NQold = len(oldgases)
        NQnew = len(newgases)

        # Make sure this is right
        assert NQnew > NQold

        NQblocksold = np.ceil(NQold/10) # Number of blocks in old indist
        NQblocksnew = np.ceil(NQnew/10) # Number of blocks to have in new indist 

        fnew = open(pipelineobj.photochem_InputsDir+'NEWin.dist', 'w')

        oldgasmixings = {}
        # Keep track of the lines that are starting and ending the current block in the in.dist file
        blockstart = 0
        blockend = NZ
        spec_counter = 0
        for i in range(int(NQblocksold)):
            curr_nq_block = ascii.read(pipelineobj.photochem_InputsDir+'in.dist', data_start=blockstart, data_end=blockend)
            blockstart = blockend
            blockend = blockend+NZ

            for c in curr_nq_block.colnames:
                oldgasmixings[oldgases[spec_counter]] = curr_nq_block[c]
                spec_counter += 1

        zeroscol = oldgasmixings['HCO'] # need a column with 1e-99s 

        for i in range(int(NQblocksnew)):
            gases_to_add = newgases[i*10:i*10+10]
            write_table = Table()
            for gind, gas in enumerate(gases_to_add):
                if gas in oldgases:
                    write_table.add_column(oldgasmixings[gas], name='col'+str(gind))
                else:
                    write_table.add_column(zeroscol, name='col'+str(gind))

            for line in range(len(write_table)):
                    fnew.write('   ')
                    for col in write_table.columns:
                        if write_table[col][line] < 9e-99:
                            val = 9e-99
                        else:
                            val = write_table[col][line]
                        fnew.write("{:.8E}".format(val)+'   ')
                    fnew.write('\n')
            
        olddist_txt = open(pipelineobj.photochem_InputsDir+'in.dist', 'r')
        alllines = olddist_txt.readlines()
        olddist_txt.close()
        for line in range(len(alllines)):
            if line >= blockstart:
                fnew.write(alllines[line])

        fnew.close()
        subprocess.run('rm '+pipelineobj.photochem_InputsDir+'in.dist', shell=True)
        subprocess.run('mv '+pipelineobj.photochem_InputsDir+'NEWin.dist '+pipelineobj.photochem_InputsDir+'in.dist', shell=True)


    # For a given suite of inputs, run the photochem/climate pipeline 
    # InputFlux - list of input fluxes for both outgassing and escaping species built up in a run function 
    ##    samples for all outgassed species in their order followed by escaping species, plus a number to indicate the unique casename 
    def run_one_model(self, InputFlux, verbose=True):

        # Model ID number for file naming will be the last value of the input string
        hold = len(InputFlux)-1
        modelID = InputFlux[hold]
        fluxes = InputFlux[:hold]

        # Initialize model pipeline object
        
        # If you want to use starting points from a previous try, set self.Restart_Run to be the previous sweep name
        if type(self.Restart_Run) == str:
            if self.Starting_Point == 'Exact':
                currmodel = VPLModelingPipeline('RunNumber'+str(modelID),  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+self.Restart_Run+'/RunNumber'+str(modelID)+'/PhotochemInputs/', 
                                                verbose, find_molecules_of_interest=False, hitran_year=self.hitran_year, planet=self.planet)
            else:

                # Read in the csv file that compiled the rates used in the previous run
                input_options = pd.read_csv('/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+self.Restart_Run+'/RatesInSweep_ForFutureInputOptions.dat', delimiter=' ', index_col=['ModelNumber'])
                
                # Find the closest model
                use_starting_point = self.find_closest_prev_model(input_options, fluxes)

                currmodel = VPLModelingPipeline('RunNumber'+str(modelID),  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+use_starting_point+'/PhotochemInputs/', 
                                                verbose, find_molecules_of_interest=False, hitran_year=self.hitran_year, planet=self.planet)
                
                if self.multinest_fit_data == True:
                    input_options_climate = pd.read_csv('/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+self.multinest_input_options_climate+'/RatesInSweep_ForFutureInputOptions.dat', delimiter=' ', index_col=['ModelNumber'])
                    climate_starting_point = self.find_closest_prev_model(input_options_climate, fluxes)

                    # If the climate profile is DIFFERENT than the closest matching photochemical files, need to change the T column of in.dist AFTER pipeline vars are set
                    currmodel.multinest_climate_copycase = climate_starting_point

        else:
            currmodel = VPLModelingPipeline('RunNumber'+str(modelID), self.photochemInitial, verbose, find_molecules_of_interest=False, hitran_year=self.hitran_year, planet=self.planet)

        # Set relevant values of object 
        self.set_pipeline_vars('RunNumber'+str(modelID), currmodel)

        # Save fluxes for later
        currmodel.fluxes_used_in_sweep = fluxes

        # Need to replace flux values in species.dat for this run 
        # A pipeline object will copy new files if currmodel.photochem_InputsDir is not equal to currmodel.photochemInitialInput
        # set them equal after setting up the current models files before running automatic pipeline (this is handled in replace_fluxes function)
        self.replace_fluxes(currmodel, fluxes)

        if self.multinest_fit_data == True and climate_starting_point != use_starting_point:
            print('Model '+str(modelID)+' used different climate starting point ('+climate_starting_point+'), photochem from '+use_starting_point)
            
            # Need to change the in.dist and PLANET.dat files COMPLETELY
            subprocess.run('rm '+currmodel.photochem_InputsDir+'in.dist', shell=True)
            subprocess.run('rm '+currmodel.photochem_InputsDir+'PLANET.dat', shell=True)
            shutil.copyfile('/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+climate_starting_point+'/PhotochemInputs/in.dist', currmodel.photochem_InputsDir+'in.dist')
            shutil.copyfile('/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+climate_starting_point+'/PhotochemInputs/PLANET.dat', currmodel.photochem_InputsDir+'PLANET.dat')

            # Expand indist
            self.expand_indist(currmodel, climate_starting_point)
            #self.update_indist_T_EDD_Profiles(currmodel, climate_starting_point)
        else:
            print('Model '+str(modelID)+' has the same climate and photochem starting points ('+climate_starting_point+')')

        # If this is multinest, probably need to make sure S and C daughter species dont have deposition velocities default set 
        if self.multinest_fit_data == True:
            self.change_deposition_daughter_species(currmodel)
        
        # Run the Photochem-Climate-SMART pipeline
        converged = currmodel.run_automatic()
        #print(converged)

        # Clean abs files out as they take up the most space
        subprocess.run('rm -rf '+currmodel.LBLABC_AbsFilesDir+'*.abs', shell=True)

        # Delete copied atmos directory
        subprocess.run('rm -rf '+currmodel.atmosDir, shell=True)

        return currmodel


    # Run a Grid sweep
    # Use every combination of samples for evey outgassing and escape species of interest 
    def run_grid_sweep(self):

        ### First generate samples for all species of interest
        ### For outgassing species:
        for curr_samp in range(len(self.outgass_sample_type_gridsweep)):
            if self.outgass_sample_type_gridsweep[curr_samp] == 'Linear':
                curr_species = self.outgass_species_gridsweep[curr_samp]

                # Ensure Min / Max flux units are correct
                self.outgass_species_MinMax_gridsweep[curr_species] = self.fix_flux_units(self.outgass_species_MinMax_gridsweep[curr_species]*self.outgass_species_units_gridsweep, 
                                                                                          curr_species, 'outgass', source_type=self.outgass_species_sourcetype[curr_samp]).value

                # Linear sampling at user requested resolution
                self.outgass_samples_gridsweep[curr_species] = np.linspace(self.outgass_species_MinMax_gridsweep[curr_species][0], 
                                                                            self.outgass_species_MinMax_gridsweep[curr_species][1], 
                                                                            self.outgass_sample_resolution_gridsweep[curr_samp])
                
            elif self.outgass_sample_type_gridsweep[curr_samp] == 'Log':
                curr_species = self.outgass_species_gridsweep[curr_samp]

                # Ensure Min / Max flux units are correct
                self.outgass_species_MinMax_gridsweep[curr_species] = self.fix_flux_units(self.outgass_species_MinMax_gridsweep[curr_species]*self.outgass_species_units_gridsweep, 
                                                                                          curr_species, 'outgass', source_type=self.outgass_species_sourcetype[curr_samp]).value

                # Log sampling at user requested resolution
                self.outgass_samples_gridsweep[curr_species] = np.logspace(np.log10(self.outgass_species_MinMax_gridsweep[curr_species][0]), 
                                                                            np.log10(self.outgass_species_MinMax_gridsweep[curr_species][1]), 
                                                                            self.outgass_sample_resolution_gridsweep[curr_samp])

            elif self.outgass_sample_type_gridsweep[curr_samp] == 'UserDef':
                curr_species = self.outgass_species_gridsweep[curr_samp]

                # Fix flux units
                self.outgass_samples_gridsweep[curr_species] = self.fix_flux_units(self.outgass_samples_gridsweep[curr_species]*self.outgass_species_units_gridsweep, 
                                                                                          curr_species, 'outgass', source_type=self.outgass_species_sourcetype[curr_samp]).value
                
        ### For escaping species:
        for curr_samp in range(len(self.escape_sample_type_gridsweep)):
            if self.escape_sample_type_gridsweep[curr_samp] == 'Linear':
                curr_species = self.escape_species_gridsweep[curr_samp]

                # Ensure Min / Max flux units are correct
                self.escape_species_MinMax_gridsweep[curr_species] = self.fix_flux_units(self.escape_species_MinMax_gridsweep[curr_species]*self.escape_species_units_gridsweep[curr_samp], 
                                                                                          curr_species, 'escape', loss_type=self.escape_species_losstype[curr_samp]).value

                # Linear sampling at user requested resolution
                self.escape_samples_gridsweep[curr_species] = np.linspace(self.escape_species_MinMax_gridsweep[curr_species][0], 
                                                                            self.escape_species_MinMax_gridsweep[curr_species][1], 
                                                                            self.escape_sample_resolution_gridsweep[curr_samp])
                
            elif self.escape_sample_type_gridsweep[curr_samp] == 'Log':
                curr_species = self.escape_species_gridsweep[curr_samp]

                # Ensure Min / Max flux units are correct
                self.escape_species_MinMax_gridsweep[curr_species] = self.fix_flux_units(self.escape_species_MinMax_gridsweep[curr_species]*self.escape_species_units_gridsweep[curr_samp], 
                                                                                          curr_species, 'escape', loss_type=self.escape_species_losstype[curr_samp]).value

                # Log sampling at user requested resolution
                self.escape_samples_gridsweep[curr_species] = np.logspace(np.log10(self.escape_species_MinMax_gridsweep[curr_species][0]), 
                                                                            np.log10(self.escape_species_MinMax_gridsweep[curr_species][1]), 
                                                                            self.escape_sample_resolution_gridsweep[curr_samp])

            elif self.escape_sample_type_gridsweep[curr_samp] == 'UserDef':
                curr_species = self.escape_species_gridsweep[curr_samp]

                # Fix flux units
                self.escape_samples_gridsweep[curr_species] = self.fix_flux_units(self.escape_samples_gridsweep[curr_species]*self.escape_species_units_gridsweep[curr_samp], 
                                                                                          curr_species, 'escape', loss_type=self.escape_species_losstype[curr_samp]).value
                
        ### Define all input combinations
        ### Syntax will be samples for all outgassed species in their order followed by escaping species 
        outgass_samps = [self.outgass_samples_gridsweep[species] for species in self.outgass_samples_gridsweep.keys()]
        esc_samps = [self.escape_samples_gridsweep[species] for species in self.escape_samples_gridsweep.keys()]
        all_samps = outgass_samps+esc_samps
        inputs = [[]]
        for i in range(len(all_samps)):
            newset = []
            ns_ind = 0
            for j in all_samps[i]:
                for f in range(len(inputs)):
                    newset.append(copy.deepcopy(inputs[f]))
                    newset[ns_ind].append(j)
                    ns_ind += 1
            inputs = newset

        # add a unique ID number for each run
        for i in range(len(inputs)):
            inputs[i].append(i)

        self.gridsweep_inputstrings = inputs

        with Pool() as p:
            models = p.map(self.run_one_model, self.gridsweep_inputstrings)

        final_pressures = [m.updated_atm_pressure for m in models]
        convergence = [m.global_convergence for m in models]
        run_names = [m.casename for m in models]
        output_col_names = ['ModelNumber', 'FinalState', 'LastPsurf']
        data_for_table = []
        data_for_table.append(run_names)
        data_for_table.append(convergence)
        data_for_table.append(final_pressures)

        # Add rates used for each run
        for species in self.outgass_samples_gridsweep.keys():
            output_col_names.append(species+'_OutgassRate')
        for species in self.escape_samples_gridsweep.keys():
            output_col_names.append(species+'_EscapeRate')

        for i in range(len(models[0].fluxes_used_in_sweep)):
            hold = [m.fluxes_used_in_sweep[i] for m in models]
            data_for_table.append(hold)

        tab = Table(data_for_table, names=output_col_names)

        # Save output info
        ascii.write(tab, self.master_out+'ParameterSweep_RunStats.dat', delimiter=' ', format='fixed_width')

        '''
        print('final pressures are ')
        for press in final_pressures:
            print(press)

        print('\n')
        print('Convergence of each model:')
        for conv in convergence:
            print(conv)
        '''


    # Compile the data from a gridsweep into a python dictionary
    # Need to make sure restart_run is false when initializing a class object from scratch
    def compile_run_output(self, photochem=True):

        # put run statistics into dictionary from output file of run
        stats = ascii.read(self.master_out+'ParameterSweep_RunStats.dat')
        d = {}
        for run in range(len(stats['ModelNumber'])):
            d[stats['ModelNumber'][run]] = {}
            d[stats['ModelNumber'][run]]['FinalState'] = stats['FinalState'][run]
            d[stats['ModelNumber'][run]]['SurfacePress[bar]'] = stats['LastPsurf'][run]
            for rate in stats.colnames:
                if rate not in ['ModelNumber', 'FinalState', 'LastPsurf']:
                    d[stats['ModelNumber'][run]][rate] = stats[rate][run]

        # if photochem is True, compile the data from final PTZ mixingratios photochem output
        if photochem == True:
            for runlabel in stats['ModelNumber']:
                if d[runlabel]['FinalState'] == True:
                    ptz = ascii.read(self.master_out+runlabel+'/FINAL_PTZ_mixingratios_out.dist')
                    d[runlabel]['PTZMixingRatiosOut'] = {}
                    for col in ptz.colnames:
                        d[runlabel]['PTZMixingRatiosOut'][col] = list(ptz[col])

        # save as json
        f = open(self.sweepname+'_OutputDict.json', 'w')
        dh = json.dumps(d)
        json.dump(dh, f)
        f.close()


    # Compile data (e.g., why things failed, etc) from a failed run
    # Will also compile which runs suceeded, which timed out, etc. 
    ## Input:
    # Num_of_Models - the number of models in the master out, this just makes code easier
    def compile_info_failed_run(self, Num_of_Models=80, dict_output=True):

        # To get atm profiles and convergence info in a dict (for later plotting purposes)
        if dict_output == True:
            d = {}

        model_ID = []
        final_state = [] # 'Converged', 'Failed', or 'Timeout'
        final_pressure = [] # Last surface pressure of model

        fail_reason = [] # Reason for failure, or NaN / None
        climate_ran = [] # If climate has starting running
        clim2col_cnvtype = []

        # Set up data calls for outgassing & escape rates
        rate_cols = []
        outgass_rates = []
        escape_rates = []
        for species in self.outgass_species_gridsweep:
            rate_cols.append(species+'_OutgassRate')
            outgass_rates.append([])
        for species in self.escape_species_gridsweep:
            rate_cols.append(species+'_EscapeRate')
            escape_rates.append([])

        for i in range(1,Num_of_Models+1):

            # Get the model ID ('RunNumber#')
            model_ID_hold = 'RunNumber'+str(i)
            model_ID.append(model_ID_hold)
            path_hold = self.master_out+model_ID_hold+'/'

            if dict_output == True:
                d[model_ID_hold] = {}

            # If the run failed, try to find out why
            if os.path.exists(path_hold+'FINAL_out_FAILED.out'):
                final_state.append('Failed')
                clim2col_cnvtype.append('NA')

                if dict_output == True:
                    d[model_ID_hold]['FinalState'] = 'Failed'

                f = open(path_hold+model_ID_hold+'_SavingInfoOut.txt', 'r')
                lines = f.readlines()
                f.close()
                hold = lines[len(lines)-2]
                hold = hold.split()
                fail_reason_hold = 'Unclear'
                for lihold in lines:
                    hold4 = lihold.split()
                    if 'SGBSL' in hold4:
                        fail_reason_hold = 'SGBSL Error'
                        break
                if fail_reason_hold == 'Unclear':
                    if 'Max' in hold:
                        # If this was printed out, either failed running photochem or trying to converge on a surface pressure 
                        if 'iterations' in hold or 'Iterations' in hold: 
                            if 'inner' in hold and 'convergence' in hold:
                                fail_reason_hold = 'Failed trying to find photochem convergence with a new pressure'
                            else:
                                hold2 = lines[len(lines)-4]
                                hold2 = hold2.split()
                                if 'Photochem' in hold2:
                                    fail_reason_hold = 'Failed running photochem (without new pressure)'
                                elif 'Surf' in hold2 and 'Pressure' in hold2:
                                    hold3 = lines[len(lines)-6]
                                    hold3 = hold3.split()
                                    hold3 = "{:.5f}".format(float(hold3[len(hold3)-1]))
                                    fail_reason_hold = 'Failed trying to find new pressure, max change: '+hold3
                                else:
                                    fail_reason_hold = 'Unclear'
                                    

                    elif 'Climate' in hold and 'convergence' in hold: # failed in climate
                        fail_reason_hold = 'Failed trying to find climate convergence'

                    else:
                        fail_reason_hold = 'Unclear'

                fail_reason.append(fail_reason_hold)

            # The run was successful
            elif os.path.exists(path_hold+'FINAL_out.out'):
                if os.path.exists(path_hold+'RunSMART_'+model_ID_hold+'.run'):
                    final_state.append('Converged')
                    fail_reason.append('NA')

                    if os.path.exists(path_hold+'RunVPLClimate_2column_'+model_ID_hold+'.script'):
                        f = open(path_hold+model_ID_hold+'_SavingInfoOut.txt', 'r')
                        lines = f.readlines()
                        f.close()

                        cnvtypehold = 'Tier1'
                        for l in reversed(lines):
                            hold = l.split('2 col cnv type')
                            if len(hold) > 1:
                                cnvtypehold = l.split()[len(l.split())-1]
                                break
                        
                        clim2col_cnvtype.append(cnvtypehold)
                    
                    else:
                        clim2col_cnvtype.append('NA')
                        
                else:
                    if os.path.exists(path_hold+'RunVPLClimate_2column_'+model_ID_hold+'.script'):
                        final_state.append('Unconv2col')
                    else:
                        final_state.append('Converged')
                    fail_reason.append('NA')
                    clim2col_cnvtype.append('NA')

                if dict_output == True:
                    if os.path.exists(path_hold+'RunSMART_'+model_ID_hold+'.trnst'):
                        d[model_ID_hold]['FinalState'] = 'Converged'
                    else:
                        d[model_ID_hold]['FinalState'] = 'Unconv2col'

                    ptz = ascii.read(path_hold+'FINAL_PTZ_mixingratios_out.dist')
                    d[model_ID_hold]['PTZFile'] = {}

                    for ptzcol in ptz.colnames:
                        d[model_ID_hold]['PTZFile'][ptzcol] = list(ptz[ptzcol])

            # the run timed out
            else:
                final_state.append('Timeout')
                fail_reason.append('NA')
                clim2col_cnvtype.append('NA')

                if dict_output == True:
                    d[model_ID_hold]['FinalState'] = 'Timeout'

                    if os.path.exists(path_hold+'atmos/PHOTOCHEM/OUTPUT/PTZ_mixingratios_out.dist'):
                        try:
                            ptz = ascii.read(path_hold+'atmos/PHOTOCHEM/OUTPUT/PTZ_mixingratios_out.dist')
                            d[model_ID_hold]['PTZFile'] = {}

                            for ptzcol in ptz.colnames:
                                d[model_ID_hold]['PTZFile'][ptzcol] = list(ptz[ptzcol])
                        except:
                            pass

            # Now find the final pressure
            f = open(path_hold+'PhotochemInputs/PLANET.dat', 'r')
            lines = f.readlines()
            f.close()
            for i in range(len(lines)):
                hold = lines[i].split()
                if 'surface' in hold and 'pressure' in hold:
                    psurf_hold = float(hold[0])
                    break

            final_pressure.append(psurf_hold)

            if dict_output == True:
                d[model_ID_hold]['Psurf'] = psurf_hold

            # Now retrieve the outgassing and escape rates
            f = open(path_hold+'PhotochemInputs/species.dat', 'r')
            lines = f.readlines()
            f.close()

            for species in range(len(self.outgass_species_gridsweep)):
                gas_hold = self.outgass_species_gridsweep[species]
                for l in lines:
                    if l.split()[0][0] != '*':
                        if l.split()[0] == gas_hold:
                            if gas_hold == 'SO2':
                                outgass_rates[species].append(float(l.split()[10]))
                            else:
                                outgass_rates[species].append(float(l.split()[11]))

                            if dict_output == True:
                                d[model_ID_hold][gas_hold+'_OutgassRate'] = float(l.split()[11])

                            break

            for species in range(len(self.escape_species_gridsweep)):
                gas_hold = self.escape_species_gridsweep[species]
                esctypehold = self.escape_species_losstype[species]
                for l in lines:
                    if l.split()[0][0] != '*':
                        if l.split()[0] == gas_hold:
                            if esctypehold == 'TOA' or esctypehold == 'toa':
                                escape_rates[species].append(float(l.split()[14]))
                            elif esctypehold == 'Vdep' or esctypehold == 'vdep':
                                escape_rates[species].append(float(l.split()[9]))
                            elif esctypehold == 'Veff' or esctypehold == 'veff':
                                escape_rates[species].append(float(l.split()[15]))

                            if dict_output == True:
                                if esctypehold == 'TOA' or esctypehold == 'toa':
                                    d[model_ID_hold][gas_hold+'_EscapeRate'] = float(l.split()[14])
                                elif esctypehold == 'Vdep' or esctypehold == 'vdep':
                                    d[model_ID_hold][gas_hold+'_EscapeRate'] = float(l.split()[9])
                                elif esctypehold == 'Veff' or esctypehold == 'veff':
                                    d[model_ID_hold][gas_hold+'_EscapeRate'] = float(l.split()[15])

                            break

            # Find out if climate has run yet:
            for dirs, subdirs, files in os.walk(path_hold):
                files = files
                break

            for fih in files:
                if len(fih.split('climate')) > 1:
                    climate_ran_hold = True
                    break
                else:
                    climate_ran_hold = False
            
            climate_ran.append(climate_ran_hold)


        # Compile the information
        dat = [model_ID, final_state, clim2col_cnvtype, final_pressure, climate_ran, fail_reason]
        col_names = ['ModelNumber', 'FinalState', '2colConv', 'LastPsurf', 'ClimateRan', 'FailReason']
        for col in rate_cols:
            col_names.append(col)
        for col in outgass_rates:
            dat.append(col)
        for col in escape_rates:
            dat.append(col)
        tab = Table(dat, names=col_names)
        ascii.write(tab, self.master_out+'ParameterSweep_RunStats_failedrun.dat', delimiter=' ', format='fixed_width')

        if dict_output == True:
            f = open(self.master_out+'DataCompilation_wAtmProfiles.json', 'w')
            dh = json.dumps(d)
            json.dump(dh, f)
            f.close()

    # Compile the SMART spectra from converged runs for plotting
    def compile_smart_spectra(self, Num_of_Models=150):

        # Start dict to save output
        d = {}

        # Loop through models
        for i in range(Num_of_Models):

            # Get the model ID ('RunNumber#')
            model_ID_hold = 'RunNumber'+str(i)
            path_hold = self.master_out+model_ID_hold+'/'

            # smart spectra only created if the run converges
            if os.path.exists(path_hold+'FINAL_out.out'):
                d[model_ID_hold] = {}

                # Get the final surface pressure used
                f = open(path_hold+'PhotochemInputs/PLANET.dat', 'r')
                lines = f.readlines()
                f.close()
                for i in range(len(lines)):
                    hold = lines[i].split()
                    if 'surface' in hold and 'pressure' in hold:
                        psurf_hold = float(hold[0])
                        break

                d[model_ID_hold]['Psurf'] = psurf_hold

                # Read in the trnst spectrum
                trnst = ascii.read(path_hold+model_ID_hold+'_SMART.trnst')

                # Save trnst spectrum
                d[model_ID_hold]['Trnst_Wavlength_um'] = list(trnst['col1'])
                d[model_ID_hold]['Trnst_Depth'] = list(trnst['col4'])

                # Read in emission spectrum
                emiss = ascii.read(path_hold+model_ID_hold+'_SMART_toa.rad')

                # Save emission spectrum
                d[model_ID_hold]['TOA_Wavelength_um'] = list(emiss['col1'])
                d[model_ID_hold]['TOA_StarFlux_Wm-2um-1'] = list(emiss['col3'])
                d[model_ID_hold]['TOA_PlanetFlux_Wm-2um-1'] = list(emiss['col4'])

        f = open(self.master_out+'SpectraCompilation.json', 'w')
        dh = json.dumps(d)
        json.dump(dh, f)
        f.close()

    # Compile data outgassing/escape rates from previous runs to determine which will be closest starting points for future runs
    # Specifically used if you want to start a new run using the final state of a previous sweep but not the exact same sweep 
    ## Input:
    # Num_of_Models - the number of models in the master out, this just makes code easier
    def compile_restart_input_options(self, add_to_file=False, include_2col=True, atm_type=None):


        model_ID = []
        
        # Set up data calls for outgassing & escape rates
        rate_cols = []
        outgass_rates = []
        escape_rates = []
        for species in self.outgass_species_gridsweep:
            rate_cols.append(species+'_OutgassRate')
            outgass_rates.append([])
        for species in self.escape_species_gridsweep:
            rate_cols.append(species+'_EscapeRate')
            escape_rates.append([])

        for dirs, sdirs, fs in os.walk(self.master_out):
            break

        for model_ID_hold in sdirs:

            # Get the model ID ('RunNumber#')
            #model_ID_hold = 'Run'+str(i)
            path_hold = self.master_out+model_ID_hold+'/'

            # Check for convergence
            for dirs, subdirs, fis in os.walk(path_hold):
                break

            conv2col = True
            if include_2col == True:
                conv2col = False
                fio = open(path_hold+'/'+model_ID_hold+'_SavingInfoOut.txt', 'r')
                lines = fio.readlines()
                fio.close()
                for l in reversed(lines):
                    hold = l.split('2 column Climate convergence found')
                    if len(hold) > 1:
                        conv2col = True

            if 'FINAL_out.dist' in fis and conv2col == True:

                ptz = ascii.read(path_hold+'FINAL_PTZ_mixingratios_out.dist')

                add_to_db = True
                if atm_type == 'H2O-O2':
                    if (ptz['O'][0] + ptz['O2'][0] + ptz['H2O'][0] + ptz['O3'][0]) < 0.9:
                        add_to_db = False
                elif atm_type == 'CO2':
                    if ptz['C'][0] > 0.05:
                        add_to_db = False
                elif atm_type == 'SO2-H2O':
                    if ptz['CH4'][0] > 0.01:
                        add_to_db = False

                if add_to_db == True:
                    model_ID.append(self.sweepname+'/'+model_ID_hold)

                    # Now retrieve the outgassing and escape rates
                    f = open(path_hold+'PhotochemInputs/species.dat', 'r')
                    lines = f.readlines()
                    f.close()

                    for species in range(len(self.outgass_species_gridsweep)):
                        gas_hold = self.outgass_species_gridsweep[species]
                        sourcetypehold = self.outgass_species_sourcetype[species]
                        for l in lines:
                            if l.split()[0][0] != '*':
                                if l.split()[0] == gas_hold:
                                    if sourcetypehold == 'Flx':
                                        outgass_rates[species].append(float(l.split()[11]))
                                    elif sourcetypehold == 'FixMR':
                                        outgass_rates[species].append(float(l.split()[10]))
                                    break

                    for species in range(len(self.escape_species_gridsweep)):
                        gas_hold = self.escape_species_gridsweep[species]
                        esctypehold = self.escape_species_losstype[species]
                        for l in lines:
                            if l.split()[0][0] != '*':
                                if l.split()[0] == gas_hold:
                                    if esctypehold == 'TOA' or esctypehold == 'toa':
                                        escape_rates[species].append(float(l.split()[14]))
                                    elif esctypehold == 'Vdep' or esctypehold == 'vdep':
                                        escape_rates[species].append(float(l.split()[9]))
                                    elif esctypehold == 'Veff' or esctypehold == 'veff':
                                        escape_rates[species].append(float(l.split()[15]))

                                    break

        # Compile the information
        dat = [model_ID]
        col_names = ['ModelNumber']
        for col in rate_cols:
            col_names.append(col)
        for col in outgass_rates:
            dat.append(col)
        for col in escape_rates:
            dat.append(col)
        #tab = Table(dat, names=col_names)

        if add_to_file == False:
            tab = Table(dat, names=col_names)
            ascii.write(tab, self.master_out+'RatesInSweep_ForFutureInputOptions.dat', delimiter=' ')
        
        else:
            prevfile = ascii.read(add_to_file, delimiter=' ')

            for d in range(len(dat[0])):
                newrow = []
                for name in dat:
                    newrow.append(name[d])
                prevfile.add_row(newrow)
            
            ascii.write(prevfile, add_to_file, delimiter=' ', overwrite=True)


    # Same as the above but compiles bulk compositions 
    def compile_BulkComp_T_restart_input_options(self, add_to_file=False, include_2col=True, atm_type='H2O-O2'):


        model_ID = []
        surfpressure = []
        
        # Set up data calls for outgassing & escape rates
        species_cols = []
        for species in self.outgass_species_gridsweep:
            species_cols.append(species)
        for species in self.escape_species_gridsweep:
            species_cols.append(species)
        
        species_cols = list(set(species_cols))
        species_vmrs = [[] for s in range(len(species_cols))]

        for dirs, sdirs, fs in os.walk(self.master_out):
            break

        for model_ID_hold in sdirs:

            # Get the model ID ('RunNumber#')
            #model_ID_hold = 'Run'+str(i)
            path_hold = self.master_out+model_ID_hold+'/'

            # Check for convergence
            for dirs, subdirs, fis in os.walk(path_hold):
                break

            conv2col = True
            if include_2col == True:
                conv2col = False
                fio = open(path_hold+'/'+model_ID_hold+'_SavingInfoOut.txt', 'r')
                lines = fio.readlines()
                fio.close()
                for l in reversed(lines):
                    hold = l.split('2 column Climate convergence found')
                    if len(hold) > 1:
                        conv2col = True

            if 'FINAL_out.dist' in fis and conv2col == True:

                # Now retrieve the VMRS
                ptz = ascii.read(path_hold+'FINAL_PTZ_mixingratios_out.dist')

                add_to_db = True
                if atm_type == 'H2O-O2':
                    if (ptz['O'][0] + ptz['O2'][0] + ptz['H2O'][0] + ptz['O3'][0]) < 0.9:
                        add_to_db = False
                elif atm_type == 'CO2':
                    if ptz['C'][0] > 0.05:
                        add_to_db = False
                elif atm_type == 'SO2-H2O':
                    if ptz['CH4'][0] > 0.01:
                        add_to_db = False

                if add_to_db == True:
                    model_ID.append(self.sweepname+'/'+model_ID_hold)

                    for s in range(len(species_cols)):
                        species_vmrs[s].append(ptz[species_cols[s]][0])

                    # Now get surface pressure
                    pdat = open(path_hold+'PhotochemInputs/PLANET.dat', 'r')
                    lines = pdat.readlines()
                    for l in lines:
                        if len(l.split('surface pressure')) > 1:
                            p = float(l.split()[0])
                            break
                    surfpressure.append(p)

        # Compile the information
        dat = [model_ID, surfpressure]
        col_names = ['ModelNumber', 'SurfPress']
        for s in range(len(species_cols)):
            col_names.append(species_cols[s])
            dat.append(species_vmrs[s])

        if add_to_file == False:
            tab = Table(dat, names=col_names)
            ascii.write(tab, self.master_out+'VMRSSurfP_RatesInSweep_ForFutureInputOptions.dat', delimiter=' ')
        
        else:
            prevfile = ascii.read(add_to_file, delimiter=' ')

            for d in range(len(dat[0])):
                newrow = []
                for name in dat:
                    newrow.append(name[d])
                prevfile.add_row(newrow)
            
            ascii.write(prevfile, add_to_file, delimiter=' ', overwrite=True)



    '''
    # Find the last created climate run file 
    def find_latest_climate_filename(self, directory):
        # Regular expression to match filenames like "filename_TryXX.run"
        pattern = r"filename_Try(\d+)\.run"

        max_try_number = -1
        latest_filename = None

        # Loop through all files in the directory
        for file in os.listdir(directory):
            match = re.match(pattern, file)
            if match:
                try_number = int(match.group(1))  # Extract the Try number
                if try_number > max_try_number:
                    max_try_number = try_number
                    latest_filename = file
        
        return latest_filename


    def compile_last_available_atm_profiles(self, Num_of_Models=80):
        # Example usage
        directory = "/path/to/your/files"  # Replace with your directory path
        latest_file = self.find_latest_climate_filename(directory)

        if latest_file:
            print(f"The latest file is: {latest_file}")
        else:
            print("No matching files found.")

    '''

    # log likelihood for mcmc run
    def mcmc_lnlike(self, x):
        
        # Run the pipeline to get pressure, convergence, etc
        rng = np.random.default_rng()  # Automatically uses entropy from OS
        modelID = rng.integers(1e5)
        while os.path.exists(self.master_out+'RunNumber'+str(modelID)):
            rng = np.random.default_rng()  # Automatically uses entropy from OS
            modelID = rng.integers(1e5)

        inputfluxes = []
        for flx in x:
            inputfluxes.append(flx)
        inputfluxes.append(modelID)

        model = self.run_one_model(inputfluxes, verbose=False)

        #updated_atm_pressure, global_convergence, casename

        # Attempt to save some info to a text file for quick assessment
        outstr = str(modelID)
        outstr = outstr+' '+str(model.global_convergence)
        outstr = outstr+' '+str(model.updated_atm_pressure)
        outstr = outstr+' '+"{:.4E}".format(x[0])
        outstr = outstr+' '+"{:.4E}".format(x[1])
        outstr = outstr+' '+"{:.4E}".format(x[2])
        outstr = outstr+' '+"{:.3E}".format(x[3])
        outstr = outstr+' '+"{:.3E}".format(x[4])
        outstr = outstr+'\n'

        fsimoutputs = open(self.master_out+'EmceeSimulationOutputs.txt', 'a')
        fsimoutputs.write(outstr)
        fsimoutputs.close()

        # If the model blew up, etc, this run is discarded, likelihood set to neg infinity
        if model.global_convergence == False:
            L = -np.inf
        
        else:
            # Fit chi sq to the likelihood of 0.1 bar atmosphere +/- sigma bars
            sigma = 0.09 
            s = ((model.updated_atm_pressure-0.1)**2)/(sigma**2)
            L = -0.5*s
        
        return L

    # Priors for MCMC
    def mcmc_priors(self, x):

        # starting point
        prior = 0

        # water prior 
        if x[0] < 44552887.2545331 or x[0] > 9.47899801e11:
            prior = -np.inf
        
        # O prior
        if x[1] < 0 or x[1] > 10:
            prior = -np.inf
        
        # O2 prior
        if x[2] < 0 or x[2] > 10:
            prior = -np.inf

        # O3 prior
        if x[3] > 0.41 or x[3] < 0.003:
            prior = -np.inf
        
        # H2O2 prior
        if x[4] > 1 or x[4] < 0.00001:
            prior = -np.inf
        
        return prior
        
    # Prob function for MCMC
    def mcmc_lnprob(self, x):

        if self.mcmc_priors(x) == -np.inf:
            return -np.inf

        else:
            return self.mcmc_lnlike(x)


    # Run MCMC to find a matching pressure 
    def match_surf_pressure_MCMC(self):

        self.mcmc_pressure_only = True

        # Set up output file
        fsimoutputs = open(self.master_out+'EmceeSimulationOutputs.txt', 'w')
        fsimoutputs.write('ID Converged Psurf H2O O O2 O3 H2O2\n')
        fsimoutputs.close()

        # Set MCMC relevant parameters
        self.mcmc_ndim = 5 # H2O outgassing, O TOA loss, O2 TOA loss, H2O2 vdep, O3 vdep
        self.mcmc_nwalkers = 20
        self.mcmc_nsteps = 1000
        self.mcmc_burnin = 100 # Start with no burnin

        # Initial guesses based off of stable parameter sweep run of pressure 0.006 bar
        h2o_outgass_initguess = 228290000000.0 # molec/cm2s
        o_toaloss_initguess = 0.01 # cm/s
        o2_toaloss_initguess = 0.01 # cm/s
        o3_vdep_initguess = 0.01 # cm/s
        h2o2_vdep_initguess = 0.005 # cm/s

        # Pick initial positions for every walker minimally perturbed by 1e-4 * the mag of the initial guess
        x0 = [h2o_outgass_initguess, o_toaloss_initguess, o2_toaloss_initguess, o3_vdep_initguess, h2o2_vdep_initguess]
        pos = []
        for i in range(self.mcmc_nwalkers):
            hold = []
            hold.append(self.fix_flux_units((x0[0]+1e7*np.random.randn())*(1 / (u.cm**2 * u.s)), 'H2O', 'outgass').value)
            hold.append(self.fix_flux_units((x0[1]+1e-5*np.random.randn())*(u.cm / (u.s)), 'O', 'escape', loss_type='Veff').value)
            hold.append(self.fix_flux_units((x0[2]+1e-5*np.random.randn())*(u.cm / (u.s)), 'O2', 'escape', loss_type='Veff').value)
            hold.append(self.fix_flux_units((x0[3]+1e-5*np.random.randn())*(u.cm / (u.s)), 'O3', 'escape', loss_type='Vdep').value)
            hold.append(self.fix_flux_units((x0[4]+1e-5*np.random.randn())*(u.cm / (u.s)), 'H2O2', 'escape', loss_type='Vdep').value)

            pos.append(np.array(hold))
        
        print('Positions: '+str(len(pos)))

        lnProb = partial(self.mcmc_lnprob)
        backendfile = self.sweepname+'.h5'
        backend = emcee.backends.HDFBackend(self.master_out+backendfile)
        backend.reset(self.mcmc_nwalkers, self.mcmc_ndim)

        with Pool() as pool:
            sampler = emcee.EnsembleSampler(self.mcmc_nwalkers, self.mcmc_ndim, lnProb, backend=backend, pool=pool)
            sampler.run_mcmc(pos, self.mcmc_nsteps, progress=False)
#        samples = sampler.get_chain(discard=self.mcmc_burnin, flat=True)


    # Priors to use for PyMultiNest
    # PyMultiNest always sets a prior value to 0 - 1 and this will modify to be within flat prior
    def multinest_prior(self, cube, ndim, nparams):

        # H2O outgassing rate prior
        wat_lowlim = 44552887.2545331
        wat_hilim = 9.47899801e11
        cube[0] = (cube[0]*(wat_hilim - wat_lowlim)) + wat_lowlim

        # CO2 Flux prior
        co2_lowlim = 0
        co2_hilim = 5.49528534e10 # Upper 97.5 percentile value from Trent
        cube[1] = cube[1]*co2_hilim

        # SO2 Fixed MR prior
        cube[2] = cube[2]*0.01 # SO2 cannot exceed 1% fixed MR at the bottom layer
        
        # O effusion velocity prior
        o_lowlim = 0.001
        o_hilim = 1
        cube[3] = (cube[3]*(o_hilim - o_lowlim)) + o_lowlim
        
        # O2 effusion velocity prior
        o2_lowlim = 0.001
        o2_hilim = 0.1
        cube[4] = cube[4]*0.04

        # CO2 effusion velocity prior
        cube[5] = cube[5]*0.1

        # CO deposition velocty prior
        cube[6] = cube[6]*0.03  

        return cube
    
    def get_JWST_measurement(self, toaspec, bandpass=15):
        rp = (1.096*u.Rearth)
        rs = (0.1192*u.Rsun)
        a = 0.0158*u.AU
        fac = (rp/a)**2
        fac = fac.decompose()

        MIRI_filters, wl_filter, filt_names = open_MIRI_filters()

        if bandpass == 15:
            krn = np.interp(toaspec['col1'], wl_filter, MIRI_filters[:,5], right=0, left=0)
            krn = -krn/np.trapz(krn, toaspec['col1'])
            t1model = -np.trapz(toaspec['col3']*krn, toaspec['col1'])
        elif bandpass == 12.8:
            krn = np.interp(toaspec['col1'], wl_filter, MIRI_filters[:,4], right=0, left=0)
            krn = -krn/np.trapz(krn, toaspec['col1'])
            t1model = -np.trapz(toaspec['col3']*krn, toaspec['col1'])

        #flx = (t1spect['col2']*u.mJy).to(u.W/u.m**2/u.Hz)
        #flx = flx*(const.c/((t1spect['col1']*u.um)**2))
        #flx = flx.to(u.W/u.m**2/u.um)
        #flx = flx.value
        #flx = flx*((((12.4669*u.pc)**2)/((a)**2)).decompose().value) # Rescale to distance assumed by SMART
        #flx = np.interp(toaspec['col1'], t1spect['col1'], flx)

        t115umT = 1867
        flx = planck(t115umT, toaspec['col1']*1e-6)*np.pi/(((a/rs)**2).decompose())

        t1um = -np.trapz(flx*krn, toaspec['col1'])
        starfac = 1/(t1um/t1model)
        #print(starfac)

        hc = (const.h*const.c).to(u.J*u.um).value

        if bandpass == 15:
            krn = np.interp(toaspec['col1'][toaspec['col1'] > 1], wl_filter, MIRI_filters[:,5], right=0, left=0)
        elif bandpass == 12.8:
            krn = np.interp(toaspec['col1'][toaspec['col1'] > 1], wl_filter, MIRI_filters[:,4], right=0, left=0)
        krn = -krn/np.trapz(krn, toaspec['col1'][toaspec['col1'] > 1])

        planet = np.trapz(toaspec['col4'][toaspec['col1'] > 1]*hc/toaspec['col1'][toaspec['col1'] > 1]*krn, toaspec['col1'][toaspec['col1'] > 1])
        star = np.trapz(toaspec['col3'][toaspec['col1'] > 1]*hc/toaspec['col1'][toaspec['col1'] > 1]*krn, toaspec['col1'][toaspec['col1'] > 1])

        measurement = planet/star*fac*starfac*1e6

        return measurement.value, starfac, fac
    
    def emission_likeli(self, day, night, modid, include_12um=False):

        measd = self.get_JWST_measurement(day)[0]
        measn = self.get_JWST_measurement(night)[0]

        if include_12um == True:
            measd_12um = self.get_JWST_measurement(day, bandpass=12.8)

        # Measurements Zieba et al 2023, MG #1, ED #1, TJB #3, ZH Sin
        Fc_D = np.array([421, 392])#, 405, 410, 333.9])
        #Fc_D_err = np.array([[94, 63, 71, 110, 78.3],[94, 75, 71, 110, 78.8]]) 
        Fc_D_err = np.array([[94, 63],[94, 75]]) 

        if include_12um == True:
            Fc_D_12um = np.array([553])
            Fc_D_err_12um = np.array([[62],[72]])

        # Night measurements # MG #1, ED #1, TJB #3, ZH Sin
        Fc_N = np.array([62])#, 125, 89, 170])
        #Fc_N_err = np.array([[43, 90, 91, 92],[60, 90, 91, 31]])
        Fc_N_err = np.array([[43],[60]])

        '''
        if model.updated_atm_pressure > 0.1:
            sigma = 0.7 # okay up to 0.8 bar
        else:
            sigma = 0.06 # okay down to 0.04 bar
        '''
        L = 0
        for dm in range(len(Fc_D)):

            if measd < Fc_D[dm]:
                sigma = Fc_D_err[0][dm]
            else:
                sigma = Fc_D_err[1][dm]

            try:
                li = ((measd - Fc_D[dm])**2)/(sigma**2)
            
            except UnboundLocalError:
                print('UNBOUND ERROR: ID number '+str(modid)+', Dayside: '+str(measd))
                li = 1e10

            L = L + li
        
        for nm in range(len(Fc_N)):

            if measn < Fc_N[nm]:
                sigma = Fc_N_err[0][nm]
            else:
                sigma = Fc_N_err[1][nm]

            try:
                li = ((measn - Fc_N[nm])**2)/(sigma**2)
            except UnboundLocalError:
                print('UNBOUND ERROR: ID number '+str(modid)+', Nightside: '+str(measn))
                li = 1e10

            L = L + li

        if include_12um == True:
            for dm in range(len(Fc_D_12um)):

                if measd_12um < Fc_D_12um[dm]:
                    sigma = Fc_D_err_12um[0][dm]
                else:
                    sigma = Fc_D_err_12um[1][dm]

                try:
                    li = ((measd_12um - Fc_D_12um[dm])**2)/(sigma**2)
                
                except UnboundLocalError:
                    print('UNBOUND ERROR: ID number '+str(modid)+', 12.8 um Dayside: '+str(measd_12um))
                    li = 1e10

                L = L + li

        L = -0.5*L
        return L, measd, measn
    
    # Delete files that aren't needed:
    def clean_up_single_run_dir(self, rundir, runnum):

        for ds, sds, fis in os.walk(rundir):
            break

        fis_to_save = ['FINAL_out.dist', 'FINAL_out.out', 'FINAL_PTZ_mixingratios_out.dist', 'FINAL_photochem_run_output.run',
                       'FINAL_out_FAILED.dist', 'FINAL_out_FAILED.out', 'FINAL_PTZ_mixingratios_out_FAILED.dist', 'FINAL_photochem_run_output_FAILED.run',
                       'MixingRs_'+runnum+'.dat', 'PT_profile_dayside_'+runnum+'.pt', 'PT_profile_nightside_'+runnum+'.pt', 'PT_profile_'+runnum+'.pt', 
                       runnum+'_SavingInfoOut.txt', 'RunLBLABC_d_H2O_'+runnum+'.script']
        
        for f in fis:
            if f not in fis_to_save:
                if len(f.split('SMART')) > 1:
                    pass
                else:
                    os.remove(rundir+'/'+f)

    
    # log likelihood for PyMultiNest
    def multinest_loglike(self, cube, ndim, nparams):

        # Run the pipeline to get pressure, convergence, etc
        rng = np.random.default_rng()  # Automatically uses entropy from OS
        modelID = rng.integers(1e5)
        modelID = modelID+self.rank
        try:
            os.mkdir(self.master_out+'RunNumber'+str(modelID)+'/')
            made = True
        except FileExistsError:
            made = False
        while made == False:
            rng = np.random.default_rng()  # Automatically uses entropy from OS
            modelID = rng.integers(1e5)
            modelID = modelID + self.rank
            try:
                os.mkdir(self.master_out+'RunNumber'+str(modelID)+'/')
                made = True
            except FileExistsError:
                made = False

        # watflx = cube[0]
        # oflx = cube[1]
        # o2flx = cube[2]
        # o3flx = cube[3]
        # h2o2flx = cube[4]

        watflx = cube[0]
        co2flx = cube[1]
        so2fixmr = cube[2]
        oveff = cube[3]
        o2veff = cube[4]
        co2veff = cube[5]
        covdep = cube[6]

        #inputfluxes = [watflx, oflx, o2flx, o3flx, h2o2flx, modelID]
        
        # Updated version which fixes h2o2 and O3 vdep at 0.02, and allows co2 flux/escape and co vdep to be fitted
        inputfluxes = [watflx, co2flx, so2fixmr, oveff, o2veff, 0.02, 0.02, co2veff, covdep, modelID]

        # Run the model
        model = self.run_one_model(inputfluxes, verbose=True)

        ''' This function will match to a specific surface pressue
        # If the model blew up, etc, this run is discarded, likelihood set to neg infinity
        if model.global_convergence == False:
            L = -np.inf
        
        else:
            # Fit chi sq to the likelihood of 0.1 bar atmosphere +/- sigma bars
            if model.updated_atm_pressure > 0.1:
                sigma = 0.7 # okay up to 0.8 bar
            else:
                sigma = 0.06 # okay down to 0.04 bar
            s = ((model.updated_atm_pressure-0.1)**2)/(sigma**2)
            L = -0.5*s
        '''

        ## This function will match to the emission (day/night) and transmission data 
        if model.global_convergence == False:
            L = -1e90
            measd = 0
            measn = 0
        
        else:
            try:
                dayside_emiss = ascii.read(self.master_out+'RunNumber'+str(modelID)+'/RunNumber'+str(modelID)+'_dayside_SMART_toa.rad')
                nightside_emiss = ascii.read(self.master_out+'RunNumber'+str(modelID)+'/RunNumber'+str(modelID)+'_nightside_SMART_toa.rad')
                #transmiss = ascii.read(self.master_out+'RunNumber'+str(modelID)+'/RunNumber'+str(modelID)+'_SMART.trnst')

                L, measd, measn = self.emission_likeli(dayside_emiss, nightside_emiss, modelID, include_12um=self.multinest_use_12um)
            
            except:
                print('ASCII COULDNT READ IN SPECTRA ERROR, '+str(modelID))
                L = -1e90
                measd = 0
                measn = 0
        
        # Attempt to save some info to a text file for quick assessment
        try:
            outstr = str(modelID)
            outstr = outstr+' '+str(model.global_convergence)
            outstr = outstr+' '+str(model.multinest_climate_copycase)
            outstr = outstr+' '+str(model.updated_atm_pressure)
            outstr = outstr+' '+"{:.4E}".format(model.last_gross_err)
            outstr = outstr+' '+"{:.4E}".format(L)
            outstr = outstr+' '+"{:.4E}".format(measd)
            outstr = outstr+' '+"{:.4E}".format(measn)
            outstr = outstr+' '+"{:.4E}".format(watflx)
            outstr = outstr+' '+"{:.4E}".format(co2flx)
            outstr = outstr+' '+"{:.4E}".format(so2fixmr)
            outstr = outstr+' '+"{:.4E}".format(oveff)
            outstr = outstr+' '+"{:.4E}".format(o2veff)
            outstr = outstr+' '+"{:.4E}".format(co2veff)
            outstr = outstr+' '+"{:.3E}".format(covdep)
            outstr = outstr+'\n'

            fsimoutputs = open(self.master_out+'EmceeSimulationOutputs.txt', 'a')
            fsimoutputs.write(outstr)
            fsimoutputs.close()
        except:
            print('FILE WRITING ERROR: DID NOT WORK '+str(modelID))

        if L > 0 and L < 1e-200:
            L = 0

        ### File clean up
        # self.master_out+'RunNumber'+str(modelID)+'/'
        self.clean_up_single_run_dir(self.master_out+'RunNumber'+str(modelID), 'RunNumber'+str(modelID))

        return L


    # Run a multinest fit
    def match_data_multinest(self):

        # To skip climate and spectra
        self.mcmc_pressure_only = True
        self.multinest_fit_data = True

        # Set up output file
        self.rank = MPI.COMM_WORLD.Get_rank()
        if self.rank == 0:
            if not os.path.exists(self.master_out+'EmceeSimulationOutputs.txt'):
                fsimoutputs = open(self.master_out+'EmceeSimulationOutputs.txt', 'w')
                fsimoutputs.write('ID Converged Copyfrom Psurf GrossErr Likeli DayEm NightEm H2Oflx CO2flx SO2fixMR Oveff O2veff CO2veff COvdep\n')
                fsimoutputs.close()

        # Need to take the closest matching climate and 2 col climate profiles
        self.Starting_Point = 'Euclidean'

        parameters = ['H2OFlx', 'CO2Flx', 'SO2FixMR', 'OVeff', 'O2Veff', 'CO2Veff', 'COVdep'] # CHANGED to not select O3 or H2O2 vdep, instead selecting CO2 / CO params 
        nparams = len(parameters)

        #lnlike = partial(self.multinest_loglike)
        lnlike = lambda cube, ndim, nparams: self.multinest_loglike(cube, ndim, nparams)
        #prior = partial(self.multinest_prior)
        prior = lambda cube, ndim, nparams: self.multinest_prior(cube, ndim, nparams)

        pymultinest.run(lnlike, prior, nparams, outputfiles_basename='chain/T1cMN_', resume=True, verbose=True, evidence_tolerance=10, n_live_points=800)




