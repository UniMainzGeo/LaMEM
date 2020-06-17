# reads scaling law files and compares the order of parameters

############ INPUT ###########
# outFiles: .dat files to be read and plottet
# phases: names of phases to plot
# params: names of parameters to plot
# (if a phase does not have a parameter it will plot as NaN)
outFiles = ['ScalLaw_1Obs.dat','ScalLaw_4Obs.dat']
phases   = ['Air','Sediments','UpperCrust','MiddleCrust','LowerCrust','Mantle','Dacite','Andesite']
params   = ['delta(rho)','G','nu','Bn','n','En','Vn','eta','alpha','Cp','k','T']

####### End of Input #########

import matplotlib.pyplot as plt
import numpy as np

numPh    = len(phases)
numPa    = len(params)

# loop through runs
for run in outFiles:
    # create list  of dictionaries for phases
    Phases = []
    # create example phase
    example = {}
    for iParam in params:
        example[iParam] = [float('NaN'), float('NaN')]
    for phase in phases:
        ph         = example.copy()
        ph['Name'] = phase
        Phases.append(ph.copy())
    # read file
    file     = open(run, 'r')
    lines    = file.readlines()
    file.close
    # loop through the file and find the first line to not start with a '#'
    for iLine,line in enumerate(lines):
        if line[0] != '#':
            start = iLine
            break
    # remove header from lines
    lines    = lines[start:len(lines)]
    # loop through lines
    for iLine, line in enumerate(lines):
        # split line into columns
        line = line.split()
        # which phase does this parameter belong to
        try:
            ind = phases.index(line[6])
        except:
            print line[6]+' was not in the list of phases \n'
            continue
        # write information to correct phase
        Phases[ind][line[0]] = [iLine,float(line[2])]

    # print phase per phase
    fig = plt.figure()
    # loop over phases
    for iPh, ph in enumerate(phases):
        # get all the numbers
        nums = np.zeros((numPa))
        for iPa,pa in enumerate(params):
            nums[iPa] = abs(Phases[iPh][pa][1])
        plt.plot(range(numPa),nums,'o-',label=ph)

    plt.xticks(range(numPa),params)
    plt.xlabel('Parameter', fontsize=16)
    plt.ylabel('abs(exponent)', fontsize=16)
    plt.legend()
    fig.suptitle(run, fontsize=20)
    fig.savefig(run+'.png')
    plt.close()