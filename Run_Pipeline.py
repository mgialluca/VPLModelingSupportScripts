
from Pipeline import *

'''
pipelineobj = VPLModelingPipeline('RunNumber189', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/O2Sens4TestH2O/RunNumber189/PhotochemInputs/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')

pipelineobj.OutPath = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/O2Sens4TestH2O/RunNumber189/'
trn = '14_Subtry3'
pipelineobj.photochem_InputsDir = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/O2Sens4TestH2O/RunNumber189/PhotochemInputs/'
oldptz = '/gscratch/vsm/gialluca/VPLModelingTools_Dev/O2Sens4TestH2O/RunNumber189/atmos/PHOTOCHEM/OUTPUT/PTZ_mixingratios_out.dist'
subprocess.run('cp '+pipelineobj.photochem_InputsDir+'in.dist '+pipelineobj.photochem_InputsDir+'save_from_before_thickatmclimate_in.dist', shell=True)
dat = pipelineobj.get_final_climate_output_temp_profile(trynum=trn)
pipelineobj.update_indist_T_EDD(oldptz, dat)
'''

pipelineobj = VPLModelingPipeline('ConvergeHigherO', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/Input_ConstMix/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')




pipelineobj.run_automatic()
