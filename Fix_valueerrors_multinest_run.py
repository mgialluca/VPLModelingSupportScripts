import re
from pymultinest.analyse import Analyzer
import os
import subprocess

def fix_files(basefilename):

    indir = basefilename.split('/')[0]
    for ds, sds, fis in os.walk(indir):
        break

    # Get all the files to fix
    to_fix = []
    for f in fis:
        if len(f.split(basefilename.split('/')[1])) > 1:
            to_fix.append(f)

    pattern = re.compile(r'(-?\d+\.\d+)(-)(\d+)')

    for f in to_fix:
        fin = open(indir+f, 'r')
        fout = open(indir+f+'_fixed', 'w')
        lines = fin.readlines()
        fin.close()

        for line in lines:
            fixed = pattern.sub(r'\1e\2\3', line)
            fout.write(fixed)
        
        fout.close()

        subprocess.run('rm '+indir+f, shell=True)
        subprocess.run('mv '+indir+f+'_fixed '+indir+f, shell=True)

