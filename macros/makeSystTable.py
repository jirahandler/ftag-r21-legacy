#!/usr/bin/env python
import re
import numpy as np
import sys

f = open(sys.argv[1],'r')

lines = f.readlines()

f.close()

sysList = [
        'generator_sys',
        'd0_smearing',
        'z0_smearing',
        'd0_smearing_Zmm',
        'z0_smearing_Zmm',
        'fake',
        'sHad',
        'matInt',
        'tailRateSys',
        'corlSmr'
        ]

sysListNames = [
        'generator sys',
        'd0 smearing',
        'z0 smearing',
        'd0 smearing Zmm',
        'z0 smearing Zmm',
        'fake',
        'sHad',
        'matInt',
        'tailRateSys',
        'correl. smearing'
        ]

currentbin_pt = ''
currentbin_eta = ''

bin_index = 0

binDict = {}

for line in lines:
	if 'bin(' in line:

		binedges = re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", line)
		currentbin_pt = [binedges[0],binedges[1]]
		currentbin_eta = [binedges[2],binedges[3]]
		bin_index+=1
		binDict[bin_index] = [currentbin_pt,currentbin_eta]

	if('central_value(' in line):

		centralvalues = re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", line)
		binDict[bin_index].append(centralvalues[0])
		binDict[bin_index].append(centralvalues[1])

	for systname in sysList:
		if(systname in line):
			if( systname == 'd0_smearing' and ('d0_smearing_Zmm' in line) ):
				continue
			if( systname == 'z0_smearing' and ('z0_smearing_Zmm' in line) ):
				continue
			binDict[bin_index].append([systname,re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", line)])

topline = ''
sfLine = ''
statErrLine = ''
systErrLine = ''
TotErrLine = ''

systLines = []

for syst in sysList:
	systLines.append('')

for binind in binDict:

	whichEtaBin = int(sys.argv[2]) # 1 for eta between 0 and 1.2, and 2 for eta between 1.2 and 2.5

	if whichEtaBin == 1:
		etaBinCondition = float(binDict[binind][1][0]) < 1.0
	if whichEtaBin == 2:
		etaBinCondition = float(binDict[binind][1][0]) > 1.0

	if etaBinCondition:

		topline+='& $['+str(int(float(binDict[binind][0][0])))+';'+str(int(float(binDict[binind][0][1])))+']$'

		nominalval = float(binDict[binind][2])
		sfLine+='& '+str( binDict[binind][2]  )
		statErrLine+='& '+str( binDict[binind][3]  )

		statErr = float(binDict[binind][3])
		sysErrSum = 0


		for i, syst in enumerate(sysList):
			if('d0' in syst or 'z0' in syst):
				systLines[i]+='& '+str(binDict[binind][4+i][1][1])+'\% '
				sysErrSum+= ( nominalval* (1.0/100.0)* float( binDict[binind][4+i][1][1] ) )**2
			else:
				systLines[i]+='& '+str(binDict[binind][4+i][1][0])+'\% '
				sysErrSum+= ( nominalval* (1.0/100.0)* float( binDict[binind][4+i][1][0] ) )**2

		sysErr = np.round( np.sqrt(sysErrSum), 2)

		totalErr = np.round( np.sqrt(sysErrSum+statErr**2), 2)

		systErrLine+='& '+str(sysErr)
		TotErrLine+='& '+str(totalErr)

print '\\begin{tabular}{|l|c|c|c|c|c|c|c|c|}'
print ' \hline '

print '$p_{T}$ [GeV] ', topline, '\\\\ \\hline'
print 'SF  			', sfLine, '\\\\ \\hline'
print 'Total error	', TotErrLine, '\\\\ \\hline'
print 'Stat error 	', statErrLine, '\\\\ '
print 'Syst error 	', systErrLine, '\\\\ \\hline'

for i, syst in enumerate(sysList):
	print sysListNames[i],' ',systLines[i], '\\\\ '

print ' \\hline \n \\end{tabular} '
