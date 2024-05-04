from copy import deepcopy
from math import exp


# Constants #

ndims = 3
elems = ['C', 'H', 'O', 'N']
units = [0.01, 1.0, 1.0]

C = 0.02304572860 # us eV^2
kT = 0.02569257028 # eV ~ 298.15 K

nconcur = 4
fname = 'simplex{0:04d}'
infile = '{:s}.in'.format(fname)
outfile = '{:s}.out'.format(fname)
initjob = 'sqthis -J {0:s} -t 10-00 -c 1 --mem-per-cpu 4000 qchem.20180724 -nt 1 {1:s} {2:s}'.format(fname, infile, outfile)
orecjob = 'sqthis -J {0:s} -t 10-00 -c 4 --mem-per-cpu 4000 qchem.20180724 -nt 4 {1:s} {2:s}'.format(fname, infile, outfile)
autojob = 'qchem.20180724 -nt 4 {1:s} {2:s}'.format(fname, infile, outfile)
logfile = 'simplex.log'

keys = ['it', 'verts', 'order', 'vals', 'rec_verts', 'rec_order', 'rec_vals']

sprem = '''
$rem
mem_total               4000
mem_static              500
method                  B3LYP
basis                   6-31+G*
unrestricted            false
SCF_algorithm           RCA_DIIS
SCF_convergence         8
thresh_RCA_switch       5
thresh                  14
incfock                 0
incdft                  0
CIS_n_roots             5
CIS_singlets            true
CIS_triplets            true
CIS_convergence         6
RPA                     1
$end
'''

socrem = '''
$rem
mem_total               4000
mem_static              500
method                  B3LYP
basis                   6-31+G*
purecart                2222
unrestricted            false
SCF_algorithm           RCA_DIIS
SCF_convergence         8
thresh_RCA_switch       5
thresh                  14
incFock                 0
incDFT                  0
CIS_n_roots             5
CIS_singlets            true
CIS_triplets            true
CIS_convergence         6
RPA                     1
calc_SOC                true
print_orbitals          true
molden_format           true
$end
'''

optrem = '''
$rem
jobtype                 opt
mem_total               4000
mem_static              500
method                  B3LYP
basis                   6-31+G*
unrestricted            false
SCF_algorithm           RCA_DIIS
SCF_convergence         8
thresh_RCA_switch       5
thresh                  14
incfock                 0
incdft                  0
$end
'''

# Gaussian #

def ReadGauss(comfile):

    f = open(comfile, 'rt')

    zmat = []
    dofs = []
    for line in f:
        data = line.split()
        if (len(data) > 0) and (data[0] in elems):
            row = [data[0]]
            for i in range(min(len(zmat), ndims)):
                row.append(int(data[2*i+1]))
                if data[2*i+2].endswith('!'):
                    dofs.append([len(zmat), i])
                    row.append(float(data[2*i+2][:-1]))
                else:
                    row.append(float(data[2*i+2]))
            zmat.append(row)

    f.close()

    return zmat, dofs


# Q-Chem #

def WriteQChem(zmat, infile, rem=sprem):

    f = open(infile, 'wt')

    f.write('$molecule\n')
    f.write('0 1\n')
    for row in zmat:
        f.write(row[0])
        for x in row[1:]:
            f.write(' {:.10g}'.format(x))
        f.write('\n')
    f.write('$end\n')
    f.write(rem)

    f.close()


def GetZmat(outfile):

    f = open(outfile, 'rt')

    rw = 0
    zmat = []
    while len(zmat) == 0:
        for line in f:
            if rw == 0:
                if line.lstrip().startswith('Z-mat'):
                    rw = 1
            elif rw == 1:
                if line.lstrip().startswith('$mole'):
                    rw = 2
            elif rw == 2:
                if line.lstrip().startswith('0 1'):
                    rw = 3
            elif rw == 3:
                if line.lstrip().startswith('$end'):
                    break
                else:
                    data = line.split()
                    row = [data[0]]
                    for i in range(min(len(zmat), ndims)):
                        row.append(int(data[2*i+1]))
                        row.append(float(data[2*i+2]))
                    zmat.append(row)
        if rw == 0:
            rw = 1
            f.seek(0)
        elif rw == 2:
            break

    f.close()

    return zmat


def GetSCF(outfile):

    f = open(outfile, 'rt')

    scf = 0.0
    for line in f:
        if line.lstrip().startswith('Total energy in'):
            scf = float(line.split('=')[-1])
            break
    if scf == 0.0:
        print('SCF Failure: {:s}'.format(outfile))

    f.close()

    return scf


def GetRate(outfile):

    f = open(outfile, 'rt')
    lines = f.readlines()
    f.close()

    pos = -1
    for i in range(pos, len(lines)):
        if lines[i].lstrip().startswith('TDDFT Excitation Energies'):
            pos = i
            break
    if pos == -1:
        print('TDDFT Failure: {:s}'.format(outfile))
        return 0.

    nTmax = 2
    nSmax = 1
    T = []
    S = []
    osc = []
    for i in range(pos, len(lines)):
        if (len(T) < nTmax) and lines[i].rstrip().endswith('Triplet'):
            T.append(float(lines[i-2].split('=')[-1]))
        elif (len(S) < nSmax) and lines[i].rstrip().endswith('Singlet'):
            S.append(float(lines[i-2].split('=')[-1]))
            osc.append(float(lines[i+2].split(':')[-1]))
        elif (len(S) >= nSmax) and (len(T) >= nTmax):
            break

    rate = 0.
    for n in range(len(S)):
        numer = S[n]*S[n]*osc[n]/C
        denom = 1.
        for Tn in T:
                denom += 3.*exp((S[n]-Tn)/kT)
        rate += numer/denom

    return rate


# Z-matrix #

def AverageZmat(*zmats):

    mean = deepcopy(zmats[0])
    for i in range(len(mean)):
        for j in range(min(i, ndims)):
            summ = 0.0
            for zmat in zmats:
                summ += zmat[i][2*j+2]
            mean[i][2*j+2] = summ/len(zmats)

    return mean


def AddScaledZmat(zmat1, zmat2, alpha):

    affine = deepcopy(zmat1)
    for i in range(len(affine)):
        for j in range(min(i, ndims)):
            affine[i][2*j+2] += alpha*zmat2[i][2*j+2]

    return affine


def DevZmat(ref, *zmats):

    devs = [[],[],[]]

    for zmat in zmats:

        dmax = [0.0, 0.0, 0.0]
        diff = AddScaledZmat(zmat, ref, -1.0)

        for i in range(len(diff)):
            for j in range(min(i, ndims)):
                if abs(diff[i][2*j+2]) > abs(dmax[j]):
                    dmax[j] = diff[i][2*j+2]

        for i in range(ndims):
            devs[i].append(dmax[i])

    return devs


# Log #

def WriteLog(**kwargs):

    status = 0
    for key in keys:

        if key in kwargs:
            if status == 0:
                status = 1
            elif status == 2:
                print('incorrect order of log entries')

            args = kwargs.get(key)

            f = open(logfile, 'at')
            f.write('{:s}:'.format(key))
            if key == 'it':
                f.write(' {:d}'.format(args))
            elif key in ['verts', 'order', 'rec_verts', 'rec_order']:
                for x in args:
                    f.write(' {:d}'.format(x))
            elif key in ['vals', 'rec_vals']:
                for x in args:
                    f.write(' {:g}'.format(x))
            f.write('\n')
            f.close()

        elif status == 1:
            status = 2
            

def ReadLog():

    f = open(logfile, 'rt')
    lines = f.readlines()
    f.close()

    pos = 0
    count = 0
    for line in lines[pos:]:
        if line.lstrip().startswith('it'):
            pos = count
        count += 1

    logs = []
    for line in lines[pos:]:
        key = line.lstrip().split(':')[0]
        if key in keys:
            data = []
            raw = line.split()[1:]

            if key in ['it', 'verts', 'order', 'rec_verts', 'rec_order']:
                for x in raw:
                    data.append(int(x))
            elif key in ['vals', 'rec_vals']:
                for x in raw:
                    data.append(float(x))

            if key == 'it':
                logs.append(data[0])
            else:
                logs.append(data)

    return logs

