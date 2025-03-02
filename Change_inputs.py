import subprocess
import os

for i in range(1, 9):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber0/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    

for i in [11, 12]:

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber10/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    

for i in [13, 16]:

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber14/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    

for i in range(32, 50):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber31/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
   
for i in range(50, 59):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber0/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    

for i in range(59, 63):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber63/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    

for i in range(82, 100):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber81/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
   
for i in range(100, 106):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber106/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    

for i in range(107, 109):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber106/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    

for i in range(111, 114):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber110/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    
for i in range(132, 150):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/VPLModelingSupportScripts/InputDir_ForTesting/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1T2/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
   