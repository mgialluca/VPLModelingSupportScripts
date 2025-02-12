
from Pipeline import *

test_object = VPLModelingPipeline('TestSMART', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/Bodies/T1c/T1cOutgas_Testing/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')

test_object.adjust_atmospheric_pressure = False
test_object.max_iterations_master = 10

test_object.run_automatic()