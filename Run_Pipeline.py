
from Pipeline import *

test_object = VPLModelingPipeline('TestSO2withClim', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/Bodies/Earth_Test_SO2/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')

test_object.run_automatic()