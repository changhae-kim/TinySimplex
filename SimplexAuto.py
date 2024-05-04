#!/usr/bin/python3

from SimplexUtil import ReadLog, GetRate, WriteLog, GetZmat, AverageZmat, DevZmat, AddScaledZmat, WriteQChem
from SimplexUtil import ndims, units, nconcur, infile, outfile, logfile, autojob

from copy import deepcopy

import os
import threading


# Constants #

maxcycle = 100
alphas = [-0.5, 0.5, 1.0, 2.0]
tol = 0.5
minconv = 3


# Begin Cycles #

for i in range(maxcycle):


    # Get Vertices and Sort by Values #

    it, verts = ReadLog()

    vals = []
    for vert in verts:
        val = GetRate(outfile.format(vert))
        vals.append(val)

    vals, order = (list(x) for x in zip(*reversed(sorted(zip(vals, verts)))))
    WriteLog(order=order, vals=vals)


    # Check Stopping Criteria and Get Candidate Vertices #

    zmats = []
    for vert in order:
        zmat = GetZmat(infile.format(vert))
        zmats.append(zmat)

    conv = 0
    devs = DevZmat(AverageZmat(*zmats), *zmats)
    for i in range(ndims):
        if max(devs[i]) < tol*units[i] and -min(devs[i]) < tol*units[i]:
            conv += 1

    if conv >= minconv:
        f = open(logfile, 'at')
        f.write('Converged!\n')
        f.close()
        exit()

    cen = AverageZmat(*zmats[:-1])
    vec = AddScaledZmat(cen, zmats[-1], -1.0)


    # Reflect, Expand, and Contract #

    rec_vals = []

    nruns = max(verts)
    while os.path.exists(infile.format(nruns)):
        nruns += 1
    nrecs = 0

    threads = []
    for alpha in alphas:
        geom = AddScaledZmat(cen, vec, alpha)
        WriteQChem(geom, infile.format(nruns+nrecs))
        threads.append( threading.Thread(target=os.system, args=(autojob.format(nruns+nrecs),)) )
        nrecs += 1

    rec_verts = range(nruns, nruns+nrecs)
    WriteLog(rec_verts=rec_verts)

    for thread in threads:
        thread.start()

    for thread in threads:
        thread.join()

    for vert in range(nruns, nruns+nrecs):
        val = GetRate(outfile.format(vert))
        rec_vals.append(val)


    # Sort RECs by Values #

    rec_vals, rec_order = (list(x) for x in zip(*reversed(sorted(zip(rec_vals, rec_verts)))))
    WriteLog(rec_order=rec_order, rec_vals=rec_vals)


    # Determine Shrink #

    if rec_vals[0] > vals[-2]:

        verts = deepcopy(order)
        verts[-1] = rec_order[0]
        WriteLog(it=it+1, verts=verts)
        continue

    else:

        nruns += nrecs

        verts = [vert for vert in range(nruns, nruns+len(order))]
        verts[-1] = order[0]
        WriteLog(it=it+1, verts=verts)

        nshrnk = 0
        threads = []
        for zmat in zmats[1:]:
            geom = AverageZmat(zmat, zmats[0])
            WriteQChem(geom, infile.format(nruns+nshrnk))
            threads.append(threading.Thread( target=os.system, args=(autojob.format(nruns+nshrnk),) ))
            nshrnk += 1

        nbatches = nshrnk//nconcur
        if nshrnk%nconcur > 0:
            nbatches += 1

        for i in range(nbatches):

            tlist = threads[ nconcur*i : min(nconcur*(i+1), nshrnk) ]

            for thread in tlist:
                thread.start()

            for thread in tlist:
                thread.join()

        continue

