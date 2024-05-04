#!/usr/bin/python3

from SimplexUtil import ReadGauss, GetZmat, WriteQChem, WriteLog
from SimplexUtil import units, infile, initjob

from copy import deepcopy

import os
import sys


comfile = sys.argv[1]
outfile = ''
if len(sys.argv) > 2:
    outfile = sys.argv[2]


# Get Z-matrix and DOFs #

scale = 1.0 # 5.0
zmat, dofs = ReadGauss(comfile)
if len(sys.argv) > 2:
    zmat = GetZmat(outfile)


# Generate Q-Chem Input #

WriteQChem(zmat, infile.format(0))

nruns = 1
for dof in dofs:
    i = dof[0]
    j = dof[1]
    geom = deepcopy(zmat)
    geom[i][2*j+2] += scale*units[j]
    WriteQChem(geom, infile.format(nruns))
    nruns += 1

for i in range(nruns):
    os.system(initjob.format(i))


# Log #

WriteLog(it=0, verts=range(nruns))

