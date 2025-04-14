from Pipeline import *
import pandas as pd
from multiprocessing import Pool
import numpy as np

def run_one_sim(casename):

    inputsdir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/MCMCTest/'+casename+'/PhotochemInputs/'
    masterout = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/MCMCClimateCompare/'

    pipelineobj = VPLModelingPipeline(casename, 
                                  inputsdir, 
                                  True, find_molecules_of_interest=False, hitran_year='2020')

    pipelineobj.vplclimate_executable = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimateModel/vpl_climate_supernode'
    os.mkdir(masterout+casename+'/')
    shutil.copytree('/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/', masterout + casename + '/atmos/')

    pipelineobj.photochemDir = masterout+casename+'/atmos/PHOTOCHEM/' # path to PHOTOCHEM/ dir
    pipelineobj.atmosDir = masterout+casename+'/atmos/' # path to atmos/ dir
    pipelineobj.lblabcDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/' # path to lblabc/ dir (such that lblabcDir/lblabc is the executable to call)
    pipelineobj.OutPath = masterout+casename+'/' # path for the raw model run outputs (NOT for created data products like dictionaries)
    pipelineobj.DataOutPath = masterout+casename+'/' # path for created data products like dictionaries
    pipelineobj.AtmProfPath = masterout+casename+'/' # path to put atmospheric profile files (.pt files really)
    pipelineobj.BackupPhotochemRuns = False # Make backups of individual photochem runs
    pipelineobj.photochemBackupDir = masterout+casename+'/PhotochemBackup/' # path to save output from each photochem run
    pipelineobj.LBLABC_AbsFilesDir = masterout+casename+'/ABSFiles/' # path to put the created lbl .abs files in 
    pipelineobj.lblabc_RunScriptDir = masterout+casename+'/' # path to put lbl runscripts in
    pipelineobj.vplclimate_RunScriptDir = masterout+casename+'/' # path to put vpl climate runscripts in
    pipelineobj.photochem_InputsDir = masterout+casename+'/PhotochemInputs/' # The path to create new photochem inputs in
    pipelineobj.xsec_Path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/xsec/' # The path where cross section files can be found
    pipelineobj.SMART_RunScriptDir = masterout+casename+'/'

    # Adjust the atmospheric pressure
    pipelineobj.adjust_atmospheric_pressure = True
    pipelineobj.suppress_IOerrors = True
    pipelineobj.run_spectra = True
    pipelineobj.include_2column_climate = False

    pipelineobj.run_automatic()


simlist = pd.read_csv('./EmceeSimulationOutputs.txt', delimiter=' ')
simnumbers = list(simlist.index)
simstorun = np.random.choice(simlist, 190, replace=False)

cases = ['RunNumber'+str(s) for s in simstorun]

with Pool() as p:
    models = p.map(run_one_sim, cases)