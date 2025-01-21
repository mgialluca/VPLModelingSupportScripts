
from Pipeline import *

test_object = VPLModelingPipeline('bar01PresTest', 
                                  '/gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/Bodies/T1c/T1cOutgas_Testing/', 
                                  True, find_molecules_of_interest=False, hitran_year='2020')

test_object.adjust_atmospheric_pressure = True

test_object.run_automatic()