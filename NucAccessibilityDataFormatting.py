#NucAccessibilityFormatter V1
#Bastian Stark
#July 31, 2024
#Short Python script meant to take DR_SASA .atmasa output files as input and apply offsets to align nucleosome positions and average/list the atom accessibilities. Which nucleotides and which atoms are used for analysis can be modified in the "whitelist" lists at beginning of script.

import os

#modify these whitlists to change which nucleotides and which atoms to include
nucleotide_whitelist = ['DC']
atom_whitelist = ['N4']

#opens files and extracts relevant data from .atmasa files for further processing
atomslist = []
directory = os.getcwd()
for file in os.listdir(directory):
    fullfilename = file.split(".")
    filename = fullfilename[0]
    if file.endswith(".atmasa"):
        atoms = []
        with open(file, 'r') as PDBfile:
            for line in PDBfile:
                atom = line.split()
                if atom[0] in atom_whitelist and atom[1] in nucleotide_whitelist:
                    info = [atom[0], atom[1], atom[2], atom[3], atom[5]]
                    #opens offset file for alignment. All nucleosome structures will be aligned to position 0 at the dyad
                    offsetfile = open("offsets.txt", 'r')
                    for line2 in offsetfile:
                        offset = line2.split('\t')
                        structure, chain, positionOffset = offset[0], offset[1], offset[2]
                        if structure == filename and chain == info[2]:
                            info[3] = str(int(info[3]) - int(positionOffset))
                            atoms.append(info)
        for item in atoms:
            atomslist.append(item)

#generates positions for output data structure
atomslist.sort()
atomPostionDict = {}
n = -74
m = 74
while n <= m:
    atomPostionDict[n] = []
    n += 1
for key in atomPostionDict:
    for atom in atomslist:
        if atom[3] == str(key):
            atomPostionDict[key].append(float(atom[4]))

#writes output file for list of accessibilities
outfile1 = open(f"AtomAccessibilitiesList_{atom_whitelist}.txt", "w")
for key in atomPostionDict:
    for item in atomPostionDict[key]:
        outfile1.write(str(key) + '\t' + str(item) + '\n')

#generates positions for output data structure
accessibilityDict = {}
i = -74
j = 74
while i <= j:
    accessibilityDict[i] = 0
    i += 1
for keyA in accessibilityDict:
    for keyB in atomPostionDict:
        if keyA == keyB:
            if len(atomPostionDict[keyB]) > 0:
                accessibilityDict[keyA] = sum(atomPostionDict[keyB])/len(atomPostionDict[keyB])
            else:
                accessibilityDict[keyA] = 0

#writes output file for averaged accessibilities
outfile2 = open(f"AtomAccessibilitiesAveraged_{atom_whitelist}.txt", "w")
for key in accessibilityDict:
    outfile2.write(str(key) + '\t' + str(accessibilityDict[key]) + '\n')
