import subprocess
import os

for i in range(132, 150):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber131/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    

for i in [127, 128, 77, 78, 27, 28]:

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/SmallerRange/RunNumber40/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    

for i in range(100, 110):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber110/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    
for i in range(82, 100):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber81/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    
for i in range(50, 63):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber63/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    

for i in range(32, 50):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber31/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    

for i in range(0, 13):

    subprocess.run('rm /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber'+str(i)+'/PhotochemInputs/*', shell=True)
    subprocess.run('cp -r /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber13/PhotochemInputs/* /gscratch/vsm/gialluca/VPLModelingTools_Dev/Vdep4e-1/RunNumber'+str(i)+'/PhotochemInputs/', shell=True)
    