#!/usr/local/bin/python3

import os
from sys import argv 
import numpy as np


def gor_prediction(helix, strand, coil, profile, size, structures):
    SS_SEQ = ''
    PROBS = [0.0 for ss in structures]

    i, j, k = 0, 8, 9
    while k < size:
        #print(j,17, '\t', i,k)
        PROBS[0] = np.sum(helix[j:] * profile[i:k])
        PROBS[1] = np.sum(strand[j:] * profile[i:k])
        PROBS[2] = np.sum(coil[j:] * profile[i:k])
        SSres = max(PROBS)
        if SSres == PROBS[0]: SS_SEQ += 'H'
        else:
            if  SSres == PROBS[1]:  SS_SEQ += 'E'
            else:
                SS_SEQ += '-'

        j, k = j-1, k+1
        

    i, j, k = 0, 8, 17
    while k < len(profile):
        PROBS[0] = np.sum(helix * profile[i:k])
        PROBS[1] = np.sum(strand * profile[i:k])
        PROBS[2] = np.sum(coil * profile[i:k])
        
        SSres = max(PROBS)
        if SSres == PROBS[0]: SS_SEQ += 'H'
        else:
            if  SSres == PROBS[1]:  SS_SEQ += 'E'
            else:
                SS_SEQ += '-'

        i, j, k = i+1, j+1, k+1

    s, k = 0, len(profile)
    while j < k:
        #print(s,17, '\t', i, len(fasta)) 
        PROBS[0] = np.sum(helix[s:] * profile[i:k])
        PROBS[1] = np.sum(strand[s:] * profile[i:k])
        PROBS[2] = np.sum(coil[s:] * profile[i:k])
        
        SSres = max(PROBS)
        if SSres == PROBS[0]: SS_SEQ += 'H'
        else:
            if  SSres == PROBS[1]:  SS_SEQ += 'E'
            else:
                SS_SEQ += '-'

        i, j, s = i+1, j+1, s+1

    return(SS_SEQ)


if __name__ == '__main__':
    try:
        PWD = argv[1] # /Users/davidebuzzao/Bioinformatica_Bologna/Second_Year/LB2/ModuleB/Project/Data/TrainingSet/
        ID_FILE = argv[2] # /Users/davidebuzzao/Bioinformatica_Bologna/Second_Year/LB2/ModuleB/Project/Data/TrainingSet/Blasted.txt (1250 IDs)
        MATRICES_DIR = argv[3] #/Users/davidebuzzao/Bioinformatica_Bologna/Second_Year/LB2/ModuleB/Project/Data/TrainingSet/Matrix
        OUTDIR = argv[4] # /Users/davidebuzzao/Bioinformatica_Bologna/Second_Year/LB2/ModuleB/Project/predictions/
    except:
        print('Program Usage: Script.py <PWD> <ID FILE> <OUTPUT DIRECTORY>')
        raise SystemExit
    else:
        RESIDUES = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        SS = ['H', 'E', '-'] 
        SIZE = 17
        
        HELIX_INFO = np.load(MATRICES_DIR + 'H.info.npy')
        STRAND_INFO = np.load(MATRICES_DIR + 'E.info.npy')
        COIL_INFO = np.load(MATRICES_DIR + 'C.info.npy')
        
        with open(ID_FILE) as ID:
            ID_LIST = ID.read().splitlines()

        for domain in ID_LIST:
            PROFILE = np.load(PWD + 'PSIBLAST/BINprofiles/' + domain + '.npy')
            PREDICTION_FILE = OUTDIR + domain + '.gor'
            SS_SEQ = gor_prediction(HELIX_INFO, STRAND_INFO, COIL_INFO, PROFILE, SIZE, SS)
            DSSP = open(PWD + 'DSSP/ss/' + domain + '.ss').read().splitlines()
            print('>' + domain + '.predict' + '\n' + SS_SEQ)
            print(DSSP[0] + '.dssp' + '\n' + DSSP[1])
            with open(OUTDIR + domain + '.gor', 'w') as FILEOUT:
                FILEOUT.write('>' + domain + '\n' + SS_SEQ)
