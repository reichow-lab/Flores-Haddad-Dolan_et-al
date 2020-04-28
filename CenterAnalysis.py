#!/home/bassamh/miniconda2/envs/PY3/bin/python

from sys import argv
from math import sqrt
import statistics as stat

script, CentersFile, LipNum = argv

# Generate Dictionary containing all of the x,y data for each lipid, this dictionary will be used to calculate all 66 pairwise distances

Chains  =   ['A','B','C','D','E','F','G','H','I','J','K','L']

DenCen  =   {}

LipNum  =   int(LipNum)

with open(CentersFile) as centers:

    for i in range(0,len(Chains),1):

        DenCen[Chains[i]]   =   {}

        for j in range(1,(LipNum + 1),1):

            DenCen[Chains[i]][j]    =   {}

    i,j =   0,1

    for line in centers:

        if j == (LipNum + 1):

            i   +=  1
            j   =   1

        val =   line.split()

        DenCen[Chains[i]][j]['x']   =   float(val[0])
        DenCen[Chains[i]][j]['y']   =   float(val[1])

        j   +=  1

# Generate dictionary of all the pairwise distances

PWDist  =   {}

for i in range(0,len(Chains),1):

    PWDist[Chains[i]]   =   {}

# Calculate pair-wise distances, and populate the dictionary

for k in range(0,len(Chains),1):

    for i in range(1,(LipNum + 1),1):

        for j in range(i+1,(LipNum + 1),1):

            dist    =    sqrt(((DenCen[Chains[k]][i]['x']-DenCen[Chains[k]][j]['x'])**2)+((DenCen[Chains[k]][i]['y']-DenCen[Chains[k]][j]['y'])**2))

            key     =   "{}-{}".format(i,j)

            PWDist[Chains[k]][key]  =   dist

# Generate dictionary of the variance in the pairwise distances

PWSDev      =   {}

PWSDev['AVG'] =   {}
PWSDev['STD'] =   {}

for i in range(1,(LipNum + 1),1):

    for j in range(i+1,(LipNum + 1),1):

        key =   "{}-{}".format(i,j)

        tmp =   []

        for k in range(0,len(Chains),1):

            tmp.append(PWDist[Chains[k]][key])

        PWSDev['STD'][key]   =    sqrt(stat.variance(tmp))

        PWSDev['AVG'][key]   =    stat.mean(tmp)

# Write out the variances to a text file

with open((CentersFile + '_PWSDeviation'), 'w') as ofile:

    ofile.write('Pair' + '\t' + 'Avg (\u212B)' + '\t\t\t' + 'SDev (\u212B)' + '\n')

    tmp =   []

    for key in PWSDev['STD']:

        ofile.write(str(key) + '\t' + str(PWSDev['AVG'][key]) + '\t' + str(PWSDev['STD'][key]) + '\n')

        tmp.append(PWSDev['STD'][key])

    avg =   stat.mean(tmp)

    ofile.write('The average standard deviation in Lipid Center assignment is ' + str(avg) + '\u212B')
