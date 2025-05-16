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

    def __init__(self, sweepname, photochemInitial, restart_run=False, starting_point='Exact', hitran_year='2020'):

        self.sweepname = sweepname # Naming convention for directory structure
        self.photochemInitial = photochemInitial # Input files for photochem to copy and change 
        self.hitran_year = hitran_year # hitran year, 2016 or 2020 (default should be latter)
        self.supernode = True # If you are using the supernode on hyak, needs to be True
        if hitran_year == '2020':
            self.lblabc_qtxt_dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/hitranQtips2020/' # For the hitran distribution you want
        elif hitran_year == '2016':
            self.lblabc_qtxt_dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/hitranQtips/'

        self.R_p = 1.097*u.Rearth # Currently radius of Trappist-1c 

        self.atmos_Dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/' # Path to dir containing atmos, will be copied for runs 

        self.mcmc_pressure_only = False

        #########  Parameters to set if you want to do a grid sweep  #########

        self.outgass_species_gridsweep = ['H2O'] # Species to vary outgassing rates of
        self.outgass_species_molarmass = {} # Molar masses in g/mol
        self.outgass_species_molarmass['H2O'] = [18.015]*(u.g/u.mol)# Molar masses in g/mol
        self.escape_species_gridsweep = ['O', 'O2', 'O3', 'H2O2'] # Species to vary escape rates of
        self.escape_species_losstype = ['Veff', 'Veff', 'Vdep', 'Vdep'] # Vdep (depositional velocity at surface) or TOA (flux at top of atmosphere)
        self.escape_species_molarmass = {}
        self.escape_species_molarmass['O'] = [15.999]*(u.g/u.mol) 
        self.escape_species_molarmass['O2'] = [31.998]*(u.g/u.mol) 
        self.escape_species_molarmass['O3'] = [47.997]*(u.g/u.mol) 
        self.escape_species_molarmass['H2O2'] = [34.014]*(u.g/u.mol) 

        self.outgass_sample_type_gridsweep = ['Log'] # How to sample outgassed molecules: 'Linear', 'Log', or 'UserDef' 
        self.escape_sample_type_gridsweep = ['Linear', 'Linear', 'Linear', 'UserDef'] #['UserDef', 'UserDef', 'UserDef', 'UserDef']
        # Linear - sample every flux on a linear grid (np.linspace) with some defined resolution
        # Log - sample every flux on a log grid (np.logspace) with some defined resolution
        # UserDef - User defined arrays of samples for every flux to vary 

        # Need to set Min / Max ranges for each molecule to vary in the form of a dictionary if using Linear or Log sampling
        self.outgass_species_MinMax_gridsweep = {}
        self.outgass_species_MinMax_gridsweep['H2O'] = [44552887.2545331, 9.47899801e11]
        #[1.00028455e+11, 1.00089313e+11]
        #[9.97550516e+10, 9.99980399e+10]
        #[9.96337789e+10, 1.00608103e+11]
        #[8.48538802e+10, 9.91501611e+10] 
        #[1.07177663e+11, 1.99798009e+11]#[2.25884849e+10, 2.72793861e+11]# full bound: [44552887.2545331, 9.47899801e11] 
        #[90000000000.0, 100000000000.0] #[1.65329797e8, 3.00359578e12] # min, max

        self.escape_species_MinMax_gridsweep = {}
        self.escape_species_MinMax_gridsweep['O'] = [0,10]
        self.escape_species_MinMax_gridsweep['O2'] = [0,10]
        self.escape_species_MinMax_gridsweep['O3'] = [0,0.5]
        self.escape_species_MinMax_gridsweep['H2O2'] = [0,1]

        # Sample resolution if using Linear / Log sampling 
        self.outgass_sample_resolution_gridsweep = [7] # number of samples for each outgassed species
        self.escape_sample_resolution_gridsweep = [3,3,3,0]

        # Need to pass samples for user defined option
        self.outgass_samples_gridsweep = {}
        #self.outgass_samples_gridsweep['H2O'] = [78000000000.0]

        self.escape_samples_gridsweep = {}
        #self.escape_samples_gridsweep['O'] = [0.01, 0.1]#[0, 1e27, 1e29]#[0, 1e27, 1e28, 1e29] #[1e28, 1e29, 1e30] #[0, 1e26, 1e27] #[0, 1e23, 5e23, 1e24, 5e24, 1e25, 5e25, 1e26]
        #self.escape_samples_gridsweep['O2'] = [0.01, 0.1]#[1e26, 5e26, 1e27]
        #self.escape_samples_gridsweep['O3'] = [0.01, 0.2, 0.4] 
        self.escape_samples_gridsweep['H2O2'] = [0.02]
        


        # Units for either Min/Max values, or the user defined samples 
        self.outgass_species_units_gridsweep = 1 / (u.cm**2 * u.s) # molecules / cm2*s (can convert from mass/time with molar mass or mol/time)
        #self.escape_species_units_gridsweep = 1 / u.s # Molecules per second
        self.escape_species_units_gridsweep = [u.cm/u.s, u.cm/u.s, u.cm / u.s, u.cm/u.s] # Molecules per second

        #######################################################################

        #########  Setting output paths  #########

        self.Restart_Run = restart_run
        self.Starting_Point = starting_point # If restart_run is a string, it gives the initial files to use separately ...
        # ... if this is 'Exact', that means the runs inputs will be the exact same ...
        # ... otherwise this will point to a run statistics file to use to determine the closest available input files
        self.master_out = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+self.sweepname+'/'

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

        if self.mcmc_pressure_only == True:
            pipelineobj.include_2column_climate = False
            pipelineobj.run_spectra = False
        else:
            pipelineobj.include_2column_climate = True
            pipelineobj.run_spectra = True

        # Testing if climate executable needs to be copied
        if self.supernode == True:
            pipelineobj.vplclimate_executable = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimateModel/vpl_climate_supernode'

        # Molecules for the type of atmosphere we're interested in 

        pipelineobj.molecule_dict = {} # key-value pairs of molecules of interest (keys, str) and their hitran codes (value, int)
        gas_names = ['O2', 'H2O', 'O3']
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
    def fix_flux_units(self, var, species, fluxtype, loss_type='TOA'):

        if (fluxtype == 'escape' and loss_type == 'TOA') or fluxtype == 'outgass':
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

    # Euclidean distance metric to find out which previous run is closest to the current one
    def euclidean_distance(self, a, b):
        return np.sqrt(np.sum((a - b) ** 2))
    
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

        print('Closest Model: '+closestmodel)

        return closestmodel



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
                                                verbose, find_molecules_of_interest=False, hitran_year=self.hitran_year)
            else:

                # Read in the csv file that compiled the rates used in the previous run
                input_options = pd.read_csv('/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+self.Restart_Run+'/RatesInSweep_ForFutureInputOptions.dat', delimiter=' ', index_col=['ModelNumber'])
                
                # Find the closest model
                use_starting_point = self.find_closest_prev_model(input_options, fluxes)

                currmodel = VPLModelingPipeline('RunNumber'+str(modelID),  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+use_starting_point+'/PhotochemInputs/', 
                                                verbose, find_molecules_of_interest=False, hitran_year=self.hitran_year)

        else:
            currmodel = VPLModelingPipeline('RunNumber'+str(modelID), self.photochemInitial, verbose, find_molecules_of_interest=False, hitran_year=self.hitran_year)

        # Set relevant values of object 
        self.set_pipeline_vars('RunNumber'+str(modelID), currmodel)

        # Save fluxes for later
        currmodel.fluxes_used_in_sweep = fluxes

        # Need to replace flux values in species.dat for this run 
        # A pipeline object will copy new files if currmodel.photochem_InputsDir is not equal to currmodel.photochemInitialInput
        # set them equal after setting up the current models files before running automatic pipeline (this is handled in replace_fluxes function)
        self.replace_fluxes(currmodel, fluxes)
        
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
                                                                                          curr_species, 'outgass').value

                # Linear sampling at user requested resolution
                self.outgass_samples_gridsweep[curr_species] = np.linspace(self.outgass_species_MinMax_gridsweep[curr_species][0], 
                                                                            self.outgass_species_MinMax_gridsweep[curr_species][1], 
                                                                            self.outgass_sample_resolution_gridsweep[curr_samp])
                
            elif self.outgass_sample_type_gridsweep[curr_samp] == 'Log':
                curr_species = self.outgass_species_gridsweep[curr_samp]

                # Ensure Min / Max flux units are correct
                self.outgass_species_MinMax_gridsweep[curr_species] = self.fix_flux_units(self.outgass_species_MinMax_gridsweep[curr_species]*self.outgass_species_units_gridsweep, 
                                                                                          curr_species, 'outgass').value

                # Log sampling at user requested resolution
                self.outgass_samples_gridsweep[curr_species] = np.logspace(np.log10(self.outgass_species_MinMax_gridsweep[curr_species][0]), 
                                                                            np.log10(self.outgass_species_MinMax_gridsweep[curr_species][1]), 
                                                                            self.outgass_sample_resolution_gridsweep[curr_samp])

            elif self.outgass_sample_type_gridsweep[curr_samp] == 'UserDef':
                curr_species = self.outgass_species_gridsweep[curr_samp]

                # Fix flux units
                self.outgass_samples_gridsweep[curr_species] = self.fix_flux_units(self.outgass_samples_gridsweep[curr_species]*self.outgass_species_units_gridsweep, 
                                                                                          curr_species, 'outgass').value
                
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

        for i in range(Num_of_Models):

            # Get the model ID ('RunNumber#')
            model_ID_hold = 'RunNumber'+str(i)
            model_ID.append(model_ID_hold)
            path_hold = self.master_out+model_ID_hold+'/'

            if dict_output == True:
                d[model_ID_hold] = {}

            # If the run failed, try to find out why
            if os.path.exists(path_hold+'FINAL_out_FAILED.out'):
                final_state.append('Failed')

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
                final_state.append('Converged')
                fail_reason.append('NA')

                if dict_output == True:
                    d[model_ID_hold]['FinalState'] = 'Converged'

                    ptz = ascii.read(path_hold+'FINAL_PTZ_mixingratios_out.dist')
                    d[model_ID_hold]['PTZFile'] = {}

                    for ptzcol in ptz.colnames:
                        d[model_ID_hold]['PTZFile'][ptzcol] = list(ptz[ptzcol])

            # the run timed out
            else:
                final_state.append('Timeout')
                fail_reason.append('NA')

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
        dat = [model_ID, final_state, final_pressure, climate_ran, fail_reason]
        col_names = ['ModelNumber', 'FinalState', 'LastPsurf', 'ClimateRan', 'FailReason']
        for col in rate_cols:
            col_names.append(col)
        for col in outgass_rates:
            dat.append(col)
        for col in escape_rates:
            dat.append(col)
        tab = Table(dat, names=col_names)
        ascii.write(tab, self.master_out+'ParameterSweep_RunStats_failedrun.dat', delimiter=' ', format='fixed_width')

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
    def compile_restart_input_options(self, Num_of_Models=80, add_to_file=False):

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

        for i in range(Num_of_Models):

            # Get the model ID ('RunNumber#')
            model_ID_hold = 'RunNumber'+str(i)
            path_hold = self.master_out+model_ID_hold+'/'

            # Check for convergence
            for dirs, subdirs, fis in os.walk(path_hold):
                break

            if 'FINAL_out.dist' in fis:

                model_ID.append(self.sweepname+'/'+model_ID_hold)

                # Now retrieve the outgassing and escape rates
                f = open(path_hold+'PhotochemInputs/species.dat', 'r')
                lines = f.readlines()
                f.close()

                for species in range(len(self.outgass_species_gridsweep)):
                    gas_hold = self.outgass_species_gridsweep[species]
                    for l in lines:
                        if l.split()[0][0] != '*':
                            if l.split()[0] == gas_hold:
                                outgass_rates[species].append(float(l.split()[11]))

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
    def multinest_prior(self, cube, ndim, nparams):

        # H2O outgassing rate prior
        wat_lowlim = 44552887.2545331
        wat_hilim = 9.47899801e11
        cube[0] = (cube[0]*(wat_hilim - wat_lowlim)) + wat_lowlim

        # O effusion velocity prior
        cube[1] = cube[1]*10
        
        # O2 effusion velocity prior
        cube[2] = cube[2]*10

        # O3 deposition velocity prior
        cube[3] = cube[3]*0.5

        # H2O2 deposition velocty prior
        #cube[4] = cube[4] # Go from 0 to 1 

        # Could add CO2 here?

        return cube
    
    # log likelihood for PyMultiNest
    def multinest_loglike(self, cube, ndim, nparams):

        # Run the pipeline to get pressure, convergence, etc
        rng = np.random.default_rng()  # Automatically uses entropy from OS
        modelID = rng.integers(1e5)
        modelID = modelID+self.rank
        while os.path.exists(self.master_out+'RunNumber'+str(modelID)):
            rng = np.random.default_rng()  # Automatically uses entropy from OS
            modelID = rng.integers(1e5)
            modelID = modelID + self.rank

        watflx = cube[0]
        oflx = cube[1]
        o2flx = cube[2]
        o3flx = cube[3]
        h2o2flx = cube[4]

        inputfluxes = [watflx, oflx, o2flx, o3flx, h2o2flx, modelID]

        # Run the model
        model = self.run_one_model(inputfluxes, verbose=True)

        # Attempt to save some info to a text file for quick assessment
        outstr = str(modelID)
        outstr = outstr+' '+str(model.global_convergence)
        outstr = outstr+' '+str(model.updated_atm_pressure)
        outstr = outstr+' '+"{:.4E}".format(watflx)
        outstr = outstr+' '+"{:.4E}".format(oflx)
        outstr = outstr+' '+"{:.4E}".format(o2flx)
        outstr = outstr+' '+"{:.3E}".format(o3flx)
        outstr = outstr+' '+"{:.3E}".format(h2o2flx)
        outstr = outstr+'\n'

        fsimoutputs = open(self.master_out+'EmceeSimulationOutputs.txt', 'a')
        fsimoutputs.write(outstr)
        fsimoutputs.close()

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
        
        return L


    # Run a multinest fit
    def match_data_multinest(self):

        # To skip climate and spectra
        self.mcmc_pressure_only = True

        # Set up output file
        self.rank = MPI.COMM_WORLD.Get_rank()
        if self.rank == 0:
            fsimoutputs = open(self.master_out+'EmceeSimulationOutputs.txt', 'w')
            fsimoutputs.write('ID Converged Psurf H2O O O2 O3 H2O2\n')
            fsimoutputs.close()

        # Need to take the closest matching climate and 2 col climate profiles
        #self.Starting_Point = 'Euclidean'

        parameters = ['H2OFlx', 'OVeff', 'O2Veff', 'O3Vdep', 'H2O2Vdep']
        nparams = len(parameters)

        #lnlike = partial(self.multinest_loglike)
        lnlike = lambda cube, ndim, nparams: self.multinest_loglike(cube, ndim, nparams)
        #prior = partial(self.multinest_prior)
        prior = lambda cube, ndim, nparams: self.multinest_prior(cube, ndim, nparams)

        pymultinest.run(lnlike, prior, nparams, outputfiles_basename='chain/Test_Run_Multinest_', resume=False, verbose=True, evidence_tolerance=1, n_live_points=800)




