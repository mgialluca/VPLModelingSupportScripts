from Pipeline import *
import copy
import subprocess
import os
from multiprocessing import Pool

master = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimMN/'
#master = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimTestMulti/'

def set_pipeline_vars(casename, pipelineobj, master_out=master):

    # Paths are the main thing to set, because they will be massive amounts of running/files, want to keep each sweep colocated in one master dir
    atmos_Dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/megan_atmos/atmos/'
    # Create casename dir
    #if not os.path.exists(master_out+casename+'/'):
    #    os.mkdir(master_out+casename+'/')
    
    #shutil.copytree(atmos_Dir,  master_out + casename + '/atmos/')

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
    pipelineobj.include_2column_climate = True
    pipelineobj.suppress_IOerrors = True
    pipelineobj.run_spectra = False
    pipelineobj.dayside_starting_PT = None
    pipelineobj.nightside_starting_PT = None
    pipelineobj.NewPressure_Psurf_tolerance = 0.035

    pipelineobj.adjust_atmospheric_pressure = True
    pipelineobj.suppress_IOerrors = True
    pipelineobj.MCMC_pressure_only = False
    pipelineobj.MultiNest_DataFit = False
    
    if pipelineobj.MCMC_pressure_only == True:
        pipelineobj.include_2column_climate = False
        pipelineobj.run_spectra = False
    else:
        pipelineobj.include_2column_climate = True
        pipelineobj.run_spectra = True

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
    gas_names = ['O2', 'H2O', 'O3']
    pipelineobj.molecule_dict['Gas_names'] = gas_names
    for m in range(len(gas_names)):
        pipelineobj.molecule_dict[gas_names[m]] = pipelineobj.hitran_lookup.loc[gas_names[m]]['HitranNumber']
        pipelineobj.molecule_dict[gas_names[m]+'_RmixCol'] = m+2

def edit_speciesdat(pipelineobj, h2oin, oin, o2in, o3in, h2o2in):

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
                nsp_new.write('H2O        LL  1 2 0 0 0 0    2     0.      0.      '+"{:.4E}".format(h2oin)+'  0.      0      0.      0. \n')
            elif hold[0] == 'O':
                nsp_new.write('O          LL  1 0 0 0 0 0    0     0.      0.      0.        0.      0      0.      '+"{:.3E}".format(oin)+'\n')
            elif hold[0] == 'O2':
                nsp_new.write('O2         LL  2 0 0 0 0 0    0     0.      0.      0.        0.      0      0.      '+"{:.3E}".format(o2in)+'\n')
            elif hold[0] == 'O3':
                nsp_new.write('O3         LL  3 0 0 0 0 0    0     '+"{:.3E}".format(o3in)+'  0.      0.        0.      0      0.      0. \n')
            elif hold[0] == 'H2O2':
                nsp_new.write('H2O2       LL  2 2 0 0 0 0    0     '+"{:.3E}".format(h2o2in)+' 0.      0.        0.      0      0.      0. \n')
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

    h2oinput, oin, o2in, o3in, h2o2in, modelid = inputstring

    '''
    co2_label = str(int(np.ceil(co2input*1e6)))+'ppm'
    if h2oinput == 0.001:
        h2o_label = '1e-1percent'
    else:
        h2o_label = str(int(h2oinput*1e2))+'percent'
    '''
    case = 'RunNumber'+str(modelid)
    #case = inputstring#+'T2'
    master = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/ClimMN/'

    pipelineobj = VPLModelingPipeline(case, 
                                  master+case+'/PhotochemInputs/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')
    
    set_pipeline_vars(case, pipelineobj)
    #edit_speciesdat(pipelineobj, h2oinput, oin, o2in, o3in, h2o2in)
    print(master+case)
    for sdshol, dshol, fishol in os.walk(master+case+'/'):
        print(sdshol)
        print(dshol)
        print(fishol)
        break

    if 'FINAL_PTZ_mixingratios_out.dist' in fishol and 'vpl_2col_climate_output_'+case+'.run' in fishol:
        pipelineobj.global_convergence = True
        pipelineobj.clim2col_restarting = True

    
    else:
    #elif 'FINAL_PTZ_mixingratios_out_FAILED.dist' in fis:
        for fhold in fishol:
            subprocess.run('rm -rf '+master+case+'/'+fhold)
        
        atmos_Dir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/megan_atmos/atmos/'
        shutil.copytree(atmos_Dir,  master+case+'/atmos/')
        subprocess.run('rm -rf '+pipelineobj.LBLABC_AbsFilesDir+'*.abs', shell=True)


    
    pipelineobj.h2oinput = h2oinput
    pipelineobj.oinput = oin
    pipelineobj.o2input = o2in
    pipelineobj.o3input = o3in
    pipelineobj.h2o2input = h2o2in
    
    
    # Run the Photochem-Climate-SMART pipeline
    converged = pipelineobj.run_automatic()
    #print(converged)

    # Clean abs files out as they take up the most space
    subprocess.run('rm -rf '+pipelineobj.LBLABC_AbsFilesDir+'*.abs', shell=True)

    # Delete copied atmos directory
    subprocess.run('rm -rf '+pipelineobj.atmosDir, shell=True)

    return pipelineobj
    
'''
#H2O_mixing = [0.001, 0.01, 0.05, 0.1]
H2O_mixing = [0.12, 0.15, 0.17, 0.2, 0.3]
#ppms = np.linspace(5, 100, 20)
#CO2_ppm = [p*1e-6 for p in ppms]
CO2_ppm = [100e-6]
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


H2O_mixing = [0.001, 0.01]
ppms = np.linspace(110, 500, 40)
CO2_ppm = [p*1e-6 for p in ppms]
all_samps = [CO2_ppm, H2O_mixing]
inputs2= [[]]
for i in range(len(all_samps)):
    newset = []
    ns_ind = 0
    for j in all_samps[i]:
        for f in range(len(inputs2)):
            newset.append(copy.deepcopy(inputs2[f]))
            newset[ns_ind].append(j)
            ns_ind += 1
    inputs2 = newset

modelinputs = []
for i in range(len(inputs)):
    modelinputs.append(inputs[i])

for i in range(len(inputs2)):
    modelinputs.append(inputs2[i])

for i in range(len(modelinputs)):
    modelinputs[i].append(i)
'''

'''
modestorun = np.load('modes_to_try.npy')

modelinputs = []
for i in modestorun:
    hold = []
    hold.append(i[0])
    hold.append(i[1])
    hold.append(i[2])
    hold.append(i[3])
    hold.append(i[4])
    modelinputs.append(hold)

for i in range(len(modelinputs)):
    modelinputs[i].append(i)
'''
'''
need2conv = ['Run0',
 'Run47',
 'Run81',
 'Run86',
 'Run75',
 'Run84',
 'Run89',
 'Run90',
 'Run61',
 'Run63',
 'Run17',
 'Run51',
 'Run87',
 'Run50',
 'Run25',
 'Run41',
 'Run94',
 'Run18',
 'Run10',
 'Run31',
 'Run48',
 'Run21',
 'Run68',
 'Run57',
 'Run78',
 'Run56',
 'Run88',
 'Run55',
 'Run19',
 'Run93',
 'Run11',
 'Run43',
 'Run79',
 'Run24',
 'Run60',
 'Run27',
 'Run62',
 'Run97',
 'Run92',
 'Run20',
 'Run34',
 'Run49',
 'Run73',
 'Run23',
 'Run96',
 'Run67',
 'Run44',
 'Run26',
 'Run8',
 'Run9',
 'Run95',
 'Run30',
 'Run22',
 'Run52',
 'Run82',
 'Run64',
 'Run80',
 'Run77',
 'Run59',
 'Run83',
 'Run40',
 'Run91',
 'Run74']

with Pool() as p:
    models = p.map(run_one_model, need2conv)

convergence = [m.global_convergence for m in models]
run_names = [m.casename for m in models]
h2os = [m.h2oinput for m in models]
oins = [m.oinput for m in models]
o2s = [m.o2input for m in models]
o3s = [m.o3input for m in models]
h2o2s = [m.h2o2input for m in models]
psurf = [m.updated_atm_pressure for m in models]

data_for_table = [run_names, convergence, psurf, h2os, oins, o2s, o3s, h2o2s]
output_col_names = ['RunName', 'Converged', 'Psurf', 'H2OInput', 'OInput', 'O2Input', 'O3Input', 'H2O2Input']

tab = Table(data_for_table, names=output_col_names)

# Save output info
ascii.write(tab, master+'ParameterSweep_RunStats.dat', delimiter=' ', format='fixed_width')
'''

#testmodel = run_one_model([3.4339E+11, 3.146E-02, 2.359E-02, 1.854E-01, 4.376E-01, 99])

sims = ascii.read('/gscratch/vsm/gialluca/VPLModelingTools_Dev/MNEmiss2/EmceeSimulationOutputs.txt')
likeli = sims['Likeli']
likelinonan = likeli[np.where(np.isnan(likeli) != True)]

inds = []
for i in np.sort(likelinonan)[-40:]:
    inds.append(list(likeli).index(i))

inds = set(inds)
inputs = []
for i in inds:
    if sims['ID'][i] != 17825:
        string = []
        string.append(sims['H2O'][i])
        string.append(sims['O'][i])
        string.append(sims['O2'][i])
        string.append(sims['O3'][i])
        string.append(sims['H2O2'][i])
        string.append(sims['ID'][i])

        inputs.append(string)

with Pool() as p:
    models = p.map(run_one_model, inputs)

convergence = [m.global_convergence for m in models]
run_names = [m.casename for m in models]
h2os = [m.h2oinput for m in models]
oins = [m.oinput for m in models]
o2s = [m.o2input for m in models]
o3s = [m.o3input for m in models]
h2o2s = [m.h2o2input for m in models]
psurf = [m.updated_atm_pressure for m in models]

data_for_table = [run_names, convergence, psurf, h2os, oins, o2s, o3s, h2o2s]
output_col_names = ['RunName', 'Converged', 'Psurf', 'H2OInput', 'OInput', 'O2Input', 'O3Input', 'H2O2Input']

tab = Table(data_for_table, names=output_col_names)

# Save output info
ascii.write(tab, master+'ParameterSweep_RunStats.dat', delimiter=' ', format='fixed_width')
