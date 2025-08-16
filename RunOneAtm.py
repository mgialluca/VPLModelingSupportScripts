from Pipeline import *
import copy
import subprocess
import os
from multiprocessing import Pool

master = '/gscratch/vsm/gialluca/T1dAtms/'

def set_pipeline_vars(casename, pipelineobj, gas_names, master_out=master):

    # Paths are the main thing to set, because they will be massive amounts of running/files, want to keep each sweep colocated in one master dir
    atmos_Dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/megan_atmos/atmos/'
    # Create casename dir
    if not os.path.exists(master_out+casename+'/'):
        os.mkdir(master_out+casename+'/')
    
    if not os.path.exists(master_out+casename+'/atmos/'):
        shutil.copytree(atmos_Dir,  master_out + casename + '/atmos/')

    pipelineobj.photochemDir = master_out+casename+'/atmos/PHOTOCHEM/' # path to PHOTOCHEM/ dir
    pipelineobj.atmosDir = master_out+casename+'/atmos/' # path to atmos/ dir
    pipelineobj.lblabcDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/lblabc/' # path to lblabc/ dir (such that lblabcDir/lblabc is the executable to call)
    pipelineobj.OutPath = master_out+casename+'/' # path for the raw model run outputs (NOT for created data products like dictionaries)
    pipelineobj.DataOutPath = master_out+casename+'/' # path for created data products like dictionaries
    pipelineobj.AtmProfPath = master_out+casename+'/' # path to put atmospheric profile files (.pt files really)
    pipelineobj.BackupPhotochemRuns = False # Make backups of individual photochem runs
    pipelineobj.photochemBackupDir = master_out+casename+'/PhotochemBackup/' # path to save output from each photochem run
    pipelineobj.LBLABC_AbsFilesDir = master_out+casename+'/ABSFiles/' # path to put the created lbl .abs files in 
    pipelineobj.lblabc_RunScriptDir = master_out+casename+'/' # path to put lbl runscripts in
    pipelineobj.vplclimate_RunScriptDir = master_out+casename+'/' # path to put vpl climate runscripts in
    pipelineobj.photochem_InputsDir = master_out+casename+'/PhotochemInputs/' # The path to create new photochem inputs in
    pipelineobj.xsec_Path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/xsec/' # The path where cross section files can be found
    pipelineobj.SMART_RunScriptDir = master_out+casename+'/'

    # Adjust the atmospheric pressure
    pipelineobj.adjust_atmospheric_pressure = False
    pipelineobj.include_2column_climate = True
    pipelineobj.suppress_IOerrors = True
    pipelineobj.run_spectra = True
    pipelineobj.dayside_starting_PT = None
    pipelineobj.nightside_starting_PT = None
    pipelineobj.NewPressure_Psurf_tolerance = 0.035

    pipelineobj.c_NumberSolarZeniths = 4

    pipelineobj.MCMC_pressure_only = False
    pipelineobj.MultiNest_DataFit = False
    pipelineobj.rerun_smart_for_2col = True
    
    if pipelineobj.MCMC_pressure_only == True:
        pipelineobj.include_2column_climate = False
        pipelineobj.run_spectra = False
    
    if pipelineobj.MultiNest_DataFit == True:
        pipelineobj.multinest_climate_copycase = 'ClimTestMulti/Run99'
        copycase = pipelineobj.multinest_climate_copycase.split('/')[1]
        pipelineobj.dayside_starting_PT = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+pipelineobj.multinest_climate_copycase+'/PT_profile_dayside_'+copycase+'.pt'
        pipelineobj.nightside_starting_PT = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+pipelineobj.multinest_climate_copycase+'/PT_profile_nightside_'+copycase+'.pt'
        pipelineobj.run_spectra = True

        climoutcopy = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+pipelineobj.multinest_climate_copycase+'/vpl_2col_climate_output_'+copycase+'.run'
        fi = open(climoutcopy, 'r')
        lines = fi.readlines()
        fi.close()

        # Want to get the last output trop heating rate and avg flux, should be last two lines
        # so loop in reversed order, break loop after to conserve efficiency

        nightside_found = False
        for i in reversed(range(len(lines))):
            hold = lines[i].split()
            if len(hold) > 2:
                if hold[0] == 'surface:':
                    if nightside_found == False:
                        pipelineobj.surface_temp_nightside = float(hold[8])
                        nightside_found = True
                    else:
                        pipelineobj.surface_temp_dayside = float(hold[8])
                        # After retrieving surface temp for nightside, will have all values, break
                        break

    pipelineobj.vplclimate_executable = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimateModel/vpl_climate_supernode'

    # Molecules for the type of atmosphere we're interested in 

    pipelineobj.molecule_dict = {} # key-value pairs of molecules of interest (keys, str) and their hitran codes (value, int)
    #gas_names = ['O2', 'H2O', 'O3', 'CO2', 'CO', 'CH4', 'N2O']
    pipelineobj.molecule_dict['Gas_names'] = gas_names
    for m in range(len(gas_names)):
        pipelineobj.molecule_dict[gas_names[m]] = pipelineobj.hitran_lookup.loc[gas_names[m]]['HitranNumber']
        pipelineobj.molecule_dict[gas_names[m]+'_RmixCol'] = m+2


'''
# To run one atmosphere:
case = 'Baseline'

pipelineobj = VPLModelingPipeline(case, 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/atmos/PHOTOCHEM/INPUTFILES/TEMPLATES/ModernEarth/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020', planet='Earth')
    
set_pipeline_vars(case, pipelineobj)

converged = pipelineobj.run_automatic()
'''


def run_starting_points(inputs):

    master = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/T1dAtms/'
    case, gasnames = inputs
    initin = master+case+'/PhotochemInputs/'
    planet = 'T1d'

    '''
    if case == 'T1bSt':
        initin = master+'b/'
        planet = 'T1b'
    elif case == 'T1cSt':
        initin = master+'c/'
        planet = 'T1c'
    elif case == 'T1dSt':
        initin = master+'d/'
        planet = 'T1d'
    elif case == 'T1eSt':
        initin = master+'e/'
        planet = 'T1e'
    elif case == 'T1fSt':
        initin = master+'f/'
        planet = 'T1f'
    elif case == 'T1gSt':
        #initin = master+'g/'
        initin = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/UpdatedStarts/T1gSt/PhotochemInputs/'
        planet = 'T1g'
    elif case == 'T1hSt':
        initin = master+'h/'
        planet = 'T1h'
    '''

    pipelineobj = VPLModelingPipeline(case, 
                                    initin, 
                                    True, find_molecules_of_interest=False, hitran_year='2020', planet=planet)
    
    set_pipeline_vars(case, pipelineobj, gasnames)

    '''
    for sdshol, dshol, fishol in os.walk(master+case+'/'):
        fishol = fishol
        break

    for fhold in fishol:
        subprocess.run('rm '+master+case+'/'+fhold, shell=True)

    subprocess.run('rm -rf '+pipelineobj.LBLABC_AbsFilesDir+'*.abs', shell=True)
    '''

    converged = pipelineobj.run_automatic()

    return pipelineobj



inputs = [['O2CO201', ['O2', 'H2O', 'O3', 'CO2', 'CO']], 
          ['O2CO2', ['O2', 'H2O', 'O3', 'CO2', 'CO']], 
          ['O2H2O01', ['O2', 'H2O', 'O3']], 
          ['O2H2O', ['O2', 'H2O', 'O3']], 
          ['O2SO2', ['O2', 'H2O', 'O3', 'SO2', 'SO3']], 
          ['O2SO201', ['O2', 'H2O', 'O3', 'SO2', 'SO3']],
          ['PureO2', ['O2', 'H2O', 'O3']],
          ['PureO2p1', ['O2', 'H2O', 'O3']]]

with Pool() as p:
    models = p.map(run_starting_points, inputs)

