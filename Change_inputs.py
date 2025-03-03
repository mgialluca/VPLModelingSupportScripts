import subprocess
import os

for i in range(32, 50):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/InputDir_1bar/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    
for i in range(82, 100):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/InputDir_1bar/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    
for i in range(132, 150):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/InputDir_1bar/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    

for i in range(100, 105):

    subprocess.run('rm -rf /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber105/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
   
for i in range(107, 110):

    subprocess.run('rm -rf /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber106/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
   
for i in range(111, 113):

    subprocess.run('rm -rf /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber110/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
   
for i in [50, 52, 54, 55, 56, 57]:

    subprocess.run('rm -rf /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber51/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
   
for i in range(59, 63):

    subprocess.run('rm -rf /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber58/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
   
for i in [0, 1, 2, 3, 5, 6, 7]:

    subprocess.run('rm -rf /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber4/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
   
for i in [9]:

    subprocess.run('rm -rf /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber8/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
   
for i in [12, 14, 16]:

    subprocess.run('rm -rf /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber13/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T3/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
   