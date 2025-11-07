from Pipeline import *
import copy
import subprocess
import os
from multiprocessing import Pool
import os

# Currently set up for: TRAPPIST-1G


def set_pipeline_vars(casename, pipelineobj, master_out, gas_names = ['O2', 'H2O', 'O3']):

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
    pipelineobj.adjust_atmospheric_pressure = True
    pipelineobj.include_2column_climate = False # Might need to change
    pipelineobj.suppress_IOerrors = True
    pipelineobj.run_spectra = True
    pipelineobj.dayside_starting_PT = None
    pipelineobj.nightside_starting_PT = None
    pipelineobj.NewPressure_Psurf_tolerance = 0.035

    pipelineobj.c_NumberSolarZeniths = 4

    pipelineobj.MCMC_pressure_only = False
    pipelineobj.MultiNest_DataFit = False
    pipelineobj.rerun_smart_for_2col = False # Might need to change
    
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
    #, 'CO2', 'CO', 'CH4', 'N2O']
    pipelineobj.molecule_dict['Gas_names'] = gas_names
    for m in range(len(gas_names)):
        pipelineobj.molecule_dict[gas_names[m]] = pipelineobj.hitran_lookup.loc[gas_names[m]]['HitranNumber']
        pipelineobj.molecule_dict[gas_names[m]+'_RmixCol'] = m+2


# so2amount - amount of SO2 to add as a fixed MR, should be 1e-4 (100 ppm), 1e-3 (1000 ppm, or 0.1%), or 1e-2 (1%)
# atmtype - original atmtype, H2O-O2 or CO2
def edit_species(pipelineobj, so2amount, atmtype):

     # First need to create the inputs directory and copy the initial master files to change
    # Pipeline already has a function that does this:
    pipelineobj.setup_intial_photochem_dir()

    # Set pipeline objects "photochemInitialInput" dir to the inputs dir to avoid this happening again in the auto run
    pipelineobj.photochemInitialInput = pipelineobj.photochem_InputsDir

    nsp = open(pipelineobj.photochem_InputsDir+'species.dat', 'r')
    nsp_new = open(pipelineobj.photochem_InputsDir+'species_new.dat', 'w')

    lines = nsp.readlines()
    for l in lines:
        hold = l.split()
        if len(hold) > 0:
            if 'S' in hold[0] and hold[0] != 'SO2': # You have a S bearing species that isn't SO2 (that will be done separately)
                if ('C' in hold[0] and atmtype == 'CO2') or 'C' not in hold[0]: # if you have a C/S species the atm type allows CO2 so you need to set this deposition to 0, otherwise it will remain 1
                    nsp_new.write(hold[0])
                    add_spaces = 11-len(hold[0])
                    for space in range(add_spaces):
                        nsp_new.write(' ')
                    nsp_new.write(hold[1]+'  ') # will be 'LL' and then 2 spaces
                    # Now writing the 'O H C S N CL' block, each has a space after with 4 spaces after CL to get to LBOUND
                    nsp_new.write(hold[2]+' '+hold[3]+' '+hold[4]+' '+hold[5]+' '+hold[6]+' '+hold[7]+'    ')
                    nsp_new.write('0     0.      0.      0.        0.      0      0.      0.     \n') # zeros across
                
                else: # You need this to remain excluded if it is a C bearing molecule in a H2O-O2 atmosphere
                    nsp_new.write(l)

            elif hold[0] == 'SO2': # If this is SO2, add the appropriate fixed MR amount
                nsp_new.write('SO2        LL  2 0 0 1 0 0    1     0.      '+"{:.1E}".format(so2amount)+' 0.        0.      0      0.      0.    \n')

            else: # Anything else stays the same
                nsp_new.write(l)

        else:
            nsp_new.write(l)

    nsp_new.close()
    nsp.close()

    # Delete old species and rename fixed version to be species.dat
    subprocess.run('rm '+pipelineobj.photochem_InputsDir+'species.dat', shell=True)
    subprocess.run('mv '+pipelineobj.photochem_InputsDir+'species_new.dat '+pipelineobj.photochem_InputsDir+'species.dat', shell=True)


def run_one_model(inputs):

    originalpath, atmtype, so2amount, case = inputs

    splpath = originalpath.split('/')
    ogmaster = splpath[len(splpath)-3]
    ogcase = splpath[len(splpath)-2]
    initin = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+ogmaster+'/'+ogcase+'/PhotochemInputs/'
    newmaster = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/AdjSO2/'+ogmaster+'/'

    pipelineobj = VPLModelingPipeline(case, 
                                    initin, 
                                    True, find_molecules_of_interest=False, hitran_year='2020', planet='T1g')
    
    if atmtype == 'H2O-O2':
        gases = ['O2', 'H2O', 'O3', 'H2', 'SO2', 'SO'] # DOUBLE CHECK 
    elif atmtype == 'CO2':
        gases = ['O2', 'H2O', 'O3', 'CO2', 'CO', 'SO2', 'SO'] # DOUBLE CHECK 


    set_pipeline_vars(case, pipelineobj, newmaster, gas_names=gases)
    edit_species(pipelineobj, so2amount, atmtype)

    '''
    for sdshol, dshol, fishol in os.walk(master+case+'/'):
        fishol = fishol
        break

    for fhold in fishol:
        subprocess.run('rm '+master+case+'/'+fhold, shell=True)

    subprocess.run('rm -rf '+pipelineobj.LBLABC_AbsFilesDir+'*.abs', shell=True)
    '''

    converged = pipelineobj.run_automatic()
    pipelineobj.so2input = so2amount
    pipelineobj.original_atmtype = atmtype

    return pipelineobj


def populate_tracking_json(planet): # Want a function to run once that checks the final database and tracks what needs to be tested

    f = open('/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/'+planet+'_FinalSpectra_Database.json', 'r')
    fj = json.load(f)
    fj = json.loads(fj)
    f.close()

    tracking = {}

    for i in fj['H2O-O2']['AtmIDs']:

        originalpath = fj['H2O-O2']['Atm'+str(i)]['OriginalPath']
        splpath = originalpath.split('/')
        ogmaster = splpath[len(splpath)-3]
        ogcase = splpath[len(splpath)-2]

        if ogmaster not in tracking.keys():
            tracking[ogmaster] = {}
            tracking[ogmaster]['AtmIDs'] = []
            curr_id = 1
        
        else:
            curr_id = tracking[ogmaster]['AtmIDs'][len(tracking[ogmaster]['AtmIDs'])-1]+1

        if not os.path.exists('/gscratch/vsm/gialluca/VPLModelingTools_Dev/AdjSO2/'+ogmaster+'/'):
            os.mkdir('/gscratch/vsm/gialluca/VPLModelingTools_Dev/AdjSO2/'+ogmaster+'/')

        tracking[ogmaster]['Atm'+str(curr_id)] = {}
        tracking[ogmaster]['Atm'+str(curr_id)]['SO2Amount'] = 1e-4 
        tracking[ogmaster]['Atm'+str(curr_id)]['RunCompleted'] = False
        tracking[ogmaster]['Atm'+str(curr_id)]['Inputs'] = [originalpath, 'H2O-O2', 1e-4, 'RunNumber'+str(curr_id)]
        tracking[ogmaster]['AtmIDs'].append(curr_id)

        curr_id += 1
        tracking[ogmaster]['Atm'+str(curr_id)] = {}
        tracking[ogmaster]['Atm'+str(curr_id)]['SO2Amount'] = 1e-3 
        tracking[ogmaster]['Atm'+str(curr_id)]['RunCompleted'] = False
        tracking[ogmaster]['Atm'+str(curr_id)]['Inputs'] = [originalpath, 'H2O-O2', 1e-3, 'RunNumber'+str(curr_id)]
        tracking[ogmaster]['AtmIDs'].append(curr_id)

        curr_id += 1
        tracking[ogmaster]['Atm'+str(curr_id)] = {}
        tracking[ogmaster]['Atm'+str(curr_id)]['SO2Amount'] = 1e-2 
        tracking[ogmaster]['Atm'+str(curr_id)]['RunCompleted'] = False
        tracking[ogmaster]['Atm'+str(curr_id)]['Inputs'] = [originalpath, 'H2O-O2', 1e-2, 'RunNumber'+str(curr_id)]
        tracking[ogmaster]['AtmIDs'].append(curr_id)


    for i in fj['CO2']['AtmIDs']:

        originalpath = fj['CO2']['Atm'+str(i)]['OriginalPath']
        splpath = originalpath.split('/')
        ogmaster = splpath[len(splpath)-3]
        ogcase = splpath[len(splpath)-2]

        if ogmaster not in tracking.keys():
            tracking[ogmaster] = {}
            tracking[ogmaster]['AtmIDs'] = []
            curr_id = 1
        
        else:
            curr_id = tracking[ogmaster]['AtmIDs'][len(tracking[ogmaster]['AtmIDs'])-1]+1

        if not os.path.exists('/gscratch/vsm/gialluca/VPLModelingTools_Dev/AdjSO2/'+ogmaster+'/'):
            os.mkdir('/gscratch/vsm/gialluca/VPLModelingTools_Dev/AdjSO2/'+ogmaster+'/')

        tracking[ogmaster]['Atm'+str(curr_id)] = {}
        tracking[ogmaster]['Atm'+str(curr_id)]['SO2Amount'] = 1e-4 
        tracking[ogmaster]['Atm'+str(curr_id)]['RunCompleted'] = False
        tracking[ogmaster]['Atm'+str(curr_id)]['Inputs'] = [originalpath, 'CO2', 1e-4, 'RunNumber'+str(curr_id)]
        tracking[ogmaster]['AtmIDs'].append(curr_id)

        curr_id += 1
        tracking[ogmaster]['Atm'+str(curr_id)] = {}
        tracking[ogmaster]['Atm'+str(curr_id)]['SO2Amount'] = 1e-3 
        tracking[ogmaster]['Atm'+str(curr_id)]['RunCompleted'] = False
        tracking[ogmaster]['Atm'+str(curr_id)]['Inputs'] = [originalpath, 'CO2', 1e-3, 'RunNumber'+str(curr_id)]
        tracking[ogmaster]['AtmIDs'].append(curr_id)

        curr_id += 1
        tracking[ogmaster]['Atm'+str(curr_id)] = {}
        tracking[ogmaster]['Atm'+str(curr_id)]['SO2Amount'] = 1e-2 
        tracking[ogmaster]['Atm'+str(curr_id)]['RunCompleted'] = False
        tracking[ogmaster]['Atm'+str(curr_id)]['Inputs'] = [originalpath, 'CO2', 1e-2, 'RunNumber'+str(curr_id)]
        tracking[ogmaster]['AtmIDs'].append(curr_id)


    f = open('/gscratch/vsm/gialluca/VPLModelingTools_Dev/AdjSO2/'+planet+'_Tracking.json', 'w')
    dh = json.dumps(tracking)
    json.dump(dh, f)
    f.close()



#### Rerunning with SO2 fixed MRs 

# First need the tracking document

planet = 'T1g' # CHANGES PLANET TO PLANET 
avail_cores = 40 # in case we use a 40 core node

if not os.path.exists('/gscratch/vsm/gialluca/VPLModelingTools_Dev/AdjSO2/'+planet+'_Tracking.json'):
    populate_tracking_json(planet)

f = open('/gscratch/vsm/gialluca/VPLModelingTools_Dev/AdjSO2/'+planet+'_Tracking.json', 'r')
fj = json.load(f)
track = json.loads(fj)
f.close()

# Build up the inputs to run 
all_inputs_curr_iteration = []

for i in range(avail_cores):

    found_atm_to_run = False

    for mast in track.keys(): # Possible original master dirs with atmospheres left to test
        if found_atm_to_run == True:
            break

        for a in track[mast]['AtmIDs']:

            if track[mast]['Atm'+str(a)]['RunCompleted'] == False:
                found_atm_to_run = True
                all_inputs_curr_iteration.append(track[mast]['Atm'+str(a)]['Inputs'])
                track[mast]['Atm'+str(a)]['RunCompleted'] = True
                break 

    # If no atm has been found you are towards the end 
    if found_atm_to_run == False:
        print('At the end of atmospheres for this planet, running '+str(len(all_inputs_curr_iteration))+' cases')
        break

assert len(all_inputs_curr_iteration) <= avail_cores # Just to be safe 

f = open('/gscratch/vsm/gialluca/VPLModelingTools_Dev/AdjSO2/'+planet+'_Tracking.json', 'w')
dh = json.dumps(track)
json.dump(dh, f)
f.close()


# Run the models
with Pool() as p:
    models = p.map(run_one_model, all_inputs_curr_iteration)