from Pipeline import *
import copy
import subprocess
from multiprocessing import Pool

master = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/T1cComparison/'

def set_pipeline_vars(casename, pipelineobj, master_out=master):

    # Paths are the main thing to set, because they will be massive amounts of running/files, want to keep each sweep colocated in one master dir
    atmos_Dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/megan_atmos/atmos/'
    # Create casename dir
    if not os.path.exists(master_out+casename+'/'):
        os.mkdir(master_out+casename+'/')
    
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
    pipelineobj.suppress_IOerrors = True
    pipelineobj.run_spectra = True
    pipelineobj.dayside_starting_PT = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/AtmProfiles/PT_profile_dayside_T100mbar.pt'
    pipelineobj.nightside_starting_PT = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/AtmProfiles/PT_profile_nightside_T100mbar.pt'

    pipelineobj.vplclimate_executable = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimateModel/vpl_climate/vpl_climate_supernode'

    # Molecules for the type of atmosphere we're interested in 

    pipelineobj.molecule_dict = {} # key-value pairs of molecules of interest (keys, str) and their hitran codes (value, int)
    gas_names = ['O2', 'H2O', 'O3', 'CO2', 'CO']
    pipelineobj.molecule_dict['Gas_names'] = gas_names
    for m in range(len(gas_names)):
        pipelineobj.molecule_dict[gas_names[m]] = pipelineobj.hitran_lookup.loc[gas_names[m]]['HitranNumber']
        pipelineobj.molecule_dict[gas_names[m]+'_RmixCol'] = m+2

def edit_speciesdat(pipelineobj, co2in, h2oin):

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
            if hold[0] == 'H2O':
                nsp_new.write('H2O        LL  1 2 0 0 0 0    1     0.      '+str(h2oin)+'  0.        0.      0      0.      0.      \n')
            elif hold[0] == 'CO2':
                nsp_new.write('CO2        LL  2 0 1 0 0 0    1     0.      '+"{:.1e}".format(co2in)+' 0.        0.      0      0.      0.\n')
            else:
                nsp_new.write(l)
        else:
            nsp_new.write(l)
    
    nsp_new.close()
    nsp.close()

    # Delete old species and rename fixed version to be species.dat
    subprocess.run('rm '+pipelineobj.photochem_InputsDir+'species.dat', shell=True)
    subprocess.run('mv '+pipelineobj.photochem_InputsDir+'species_new.dat '+pipelineobj.photochem_InputsDir+'species.dat', shell=True)


def run_one_model(inputstring):

    co2input, h2oinput, modelid = inputstring

    '''
    co2_label = str(int(np.ceil(co2input*1e6)))+'ppm'
    if h2oinput == 0.001:
        h2o_label = '1e-1percent'
    else:
        h2o_label = str(int(h2oinput*1e2))+'percent'
    '''
    case = 'Run'+str(modelid)

    pipelineobj = VPLModelingPipeline(case, 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/T100mbar/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')
    
    set_pipeline_vars(case, pipelineobj)
    edit_speciesdat(pipelineobj, co2input, h2oinput)

    pipelineobj.co2input = co2input
    pipelineobj.h2oinput = h2oinput

    # Run the Photochem-Climate-SMART pipeline
    converged = pipelineobj.run_automatic()
    #print(converged)

    # Clean abs files out as they take up the most space
    #subprocess.run('rm -rf '+currmodel.LBLABC_AbsFilesDir+'*.abs', shell=True)

    # Delete copied atmos directory
    subprocess.run('rm -rf '+pipelineobj.atmosDir, shell=True)

    return pipelineobj
    

H2O_mixing = [0.001, 0.01, 0.05, 0.1]
ppms = np.linspace(5, 100, 20)
CO2_ppm = [p*1e-6 for p in ppms]
all_samps = [CO2_ppm, H2O_mixing]
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

for i in range(len(inputs)):
    inputs[i].append(i)


with Pool() as p:
    models = p.map(run_one_model, inputs)

convergence = [m.global_convergence for m in models]
run_names = [m.casename for m in models]
co2s = [m.co2input for m in models]
h2os = [m.h2oinput for m in models]

data_for_table = [run_names, convergence, co2s, h2os]
output_col_names = ['RunName', 'Converged', 'CO2Input', 'H2OInput']

tab = Table(data_for_table, names=output_col_names)

# Save output info
ascii.write(tab, master+'ParameterSweep_RunStats.dat', delimiter=' ', format='fixed_width')
