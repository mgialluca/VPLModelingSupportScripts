
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



to_clean = ['Bco2', 'Bco2wL', 'Bco2wL2', 'Bco2wL3', 'Bz1Re', 'Bzoom2', 'Cco2', 'Cco2wL', 'Cco2wL2', 'Cinit', 'Czoom1', 'Czoom2', 'CzoomRT', 
            'Dco2', 'Dco2rego', 'Dco2wL', 'Dco2wL2', 'Dinit', 'Dzoom1', 'Dzoom2', 'Eco2', 'Eco2lo', 'Eco2wL', 'Eco2wL2', 'Einit', 'Ezoom1', 
            'Ezoom2', 'Ezoom3', 'Fco2', 'Fco2lo', 'Fco2wL', 'Fco2wL2', 'Gco2', 'Gco2lo', 'Gco2wL', 'Gco2wL2', 'Hco2', 'Hco2lo', 'Hco2wL']

with Pool() as p:
    models = p.map(clean_one_sweep_dir, to_clean)