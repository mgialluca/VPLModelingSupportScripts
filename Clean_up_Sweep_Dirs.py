
import os
import subprocess
import shutil
from multiprocessing import Pool

# sweepdir goes on /gscratch/vsm/gialluca/VPLModelingTools_Dev/
def clean_one_sweep_dir(sweepdir):

    path = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/'+sweepdir

    for ds, sds, fis in os.walk(path):
        break 

    # sds are runnumber specific directories 
    for sd in sds:
        pathhold = path+'/'+sd
        for dshold, sdshold, fishold in os.walk(pathhold):
            break

        fis_to_save = ['vpl_climate_output_'+sd+'.run', sd+'_SavingInfoOut.txt', 'PT_profile_'+sd+'.pt', 'MixingRs_'+sd+'.dat', 'photochem_run_output_'+sd+'.run',
                   'FINAL_out.dist', 'FINAL_out.out', 'FINAL_PTZ_mixingratios_out.dist', 'PT_profile_nightside_'+sd+'.pt', 'PT_profile_nightside_'+sd+'.pt',
                   'smart_run_output_'+sd+'.run', 'smart_run_output_nightside_'+sd+'.run', 'smart_run_output_dayside_'+sd+'.run', 'FINAL_out_FAILED.dist', 'FINAL_out_FAILED.out', 
                   'FINAL_PTZ_mixingratios_out_FAILED.dist', 'vpl_2col_climate_output_'+sd+'.run']

        for f in fishold:
            if f not in fis_to_save:
                if len(f.split('vpl_2column_')) > 1:
                    pass
                elif len(f.split('SMART')) > 1:
                    pass
                else:
                    os.remove(pathhold+'/'+f)

        if 'atmos' in sdshold:
            shutil.rmtree(pathhold+'/atmos')
        
        if 'ABSFiles' in sdshold:
            shutil.rmtree(pathhold+'/ABSFiles')



to_clean = ['AdjSO2/Bz1Re', 'AdjSO2/Bzoom1', 'AdjSO2/Bzoom2', 'AdjSO2/Cco2', 'AdjSO2/Cco2wL', 'AdjSO2/Cco2wL2', 'AdjSO2/Cinit', 'AdjSO2/Czoom1', 'AdjSO2/Czoom2', 'AdjSO2/CzoomRT', 
            'AdjSO2/Dco2', 'AdjSO2/Dco2rego', 'AdjSO2/Dco2wL', 'AdjSO2/Dco2wL2', 'AdjSO2/Dinit', 'AdjSO2/Dzoom1', 'AdjSO2/Dzoom2', 'AdjSO2/Eco2', 'AdjSO2/Eco2wL', 'AdjSO2/Eco2wL2', 'AdjSO2/Einit', 'AdjSO2/Ezoom1', 
            'AdjSO2/Ezoom2', 'AdjSO2/Ezoom3', 'AdjSO2/Fco2', 'AdjSO2/Fco2wL', 'AdjSO2/Fco2wL2', 'AdjSO2/Gco2lo', 'AdjSO2/Gco2wL', 'AdjSO2/Gco2wL2', 'AdjSO2/Hco2', 'AdjSO2/Hco2lo', 'AdjSO2/Hco2wL', 'AdjSO2/Hinit']

with Pool() as p:
    models = p.map(clean_one_sweep_dir, to_clean)