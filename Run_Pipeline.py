
from Pipeline import *

test_object = VPLModelingPipeline('PsurfTest', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/InputDir_ForTesting/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')

test_object.adjust_atmospheric_pressure = True
test_object.max_iterations_master = 6

test_object.run_automatic()