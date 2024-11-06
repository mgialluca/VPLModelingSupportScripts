
from Pipeline import *

test_object = VPLModelingPipeline('TestPipelineAutoRun_EarthlikeH2O_T1c', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/Bodies/T1c/TestPipeline_T1cH2OOutgas_Initial', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')

test_object.run_automatic()