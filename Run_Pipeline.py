
from Pipeline import *

test_object = VPLModelingPipeline('NewPresTest', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/InputDir_ForTesting/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')

test_object.adjust_atmospheric_pressure = True

test_object.run_automatic()