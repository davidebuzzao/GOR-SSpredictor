#!/usr/local/bin/python3

import os
from sys import argv 
import numpy as np

def gor_training(residues, structures, fasta, dssp, profile, size, helix, strand, coil, tot, ss_count):

    i, j, k = 0, 8, 9
    while k < size:
        #print(j,17, '\t', i,k)
        if dssp[j] not in structures: 
            j, k = j-1, k+1
        else:
            if dssp[j] == 'H':
                helix[j:] += profile[i:k]
                ss_count[0] += 1
            elif dssp[j] == 'E':
                strand[j:] += profile[i:k]
                ss_count[1] += 1
            else:
                coil[j:] += profile[i:k]
                ss_count[2] += 1
            
            tot[j:] += profile[i:k]
            j, k = j-1, k+1
        

    i, j, k = 0, 8, 17
    while k < len(fasta):
        if dssp[j] not in structures: 
            i, j, k = i+1, j+1, k+1
        else:
            if dssp[j] == 'H':
                helix += profile[i:k]
                ss_count[0] += 1
            elif dssp[j] == 'E':
                strand += profile[i:k]
                ss_count[1] += 1
            else:
                coil += profile[i:k]
                ss_count[2] += 1
            
            tot += profile[i:k]
            i, j, k = i+1, j+1, k+1


    s, k = 0, len(fasta)
    while j < k:
        #print(s,17, '\t', i, len(fasta))
        if dssp[j] not in structures: 
            i, j, s = i+1, j+1, s+1        
        else:
            if dssp[j] == 'H':
                helix[s:] += profile[i:]
                ss_count[0] += 1
            elif dssp[j] == 'E':
                strand[s:] += profile[i:]
                ss_count[1] += 1
            else:
                coil[s:] += profile[i:]
                ss_count[2] += 1
            
            tot[s:] += profile[i:]
            i, j, s = i+1, j+1, s+1

    return(helix, strand, coil, tot, ss_count)


if __name__ == '__main__':
    try:
        PWD = argv[1] # /Users/davidebuzzao/Bioinformatica_Bologna/Second_Year/LB2/ModuleB/Project/Data/TrainingSet/
        ID_FILE = argv[2] # /Users/davidebuzzao/Bioinformatica_Bologna/Second_Year/LB2/ModuleB/Project/Data/TrainingSet/jpred4.domain.list
        OUTDIR = argv[3] # /Users/davidebuzzao/Bioinformatica_Bologna/Second_Year/LB2/ModuleB/Project/Data/TrainingSet/Matrix/
    except:
        print('Program Usage: Script.py <PWD> <ID FILE>')
        raise SystemExit
    else:
        RESIDUES = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        SS = ['H', 'E', '-'] 
        SIZE = 17

        ## Initialization of 17x20 matrices
        HELIX, STRAND = np.zeros((SIZE,len(RESIDUES))), np.zeros((SIZE,len(RESIDUES)))
        COIL, TOT = np.zeros((SIZE,len(RESIDUES))), np.zeros((SIZE,len(RESIDUES))) 
        
        SS_COUNT = np.zeros((len(SS)))

        with open(ID_FILE) as ID:
            ID_LIST = ID.read().splitlines()

        for domain in ID_LIST:
            FASTA_FILE = PWD + 'FASTA/' + domain + '.fasta'
            DSSP_FILE = PWD + 'DSSP/' + domain + '.dssp'
            PROFILE_FILE = PWD + 'PSIBLAST/BINprofiles/' + domain + '.npy'
            with open(FASTA_FILE) as FASTA, open(DSSP_FILE) as DSSP:
                SEQUENCE = FASTA.read().splitlines()[1]
                STRUCTURE = DSSP.read().splitlines()[1]
                PROFILE = np.load(PROFILE_FILE) 
                
                ### Frequencies MATRICES based on evolutionary information ###
                HELIX, STRAND, COIL, TOT, SS_COUNT = gor_training(RESIDUES, SS, SEQUENCE, STRUCTURE, PROFILE, SIZE, HELIX, STRAND, COIL, TOT, SS_COUNT)

        ### NORMALIZATION ###
        SS_TOT = np.sum(SS_COUNT)
        FREQ_H_norm = SS_COUNT[0]/SS_TOT
        FREQ_E_norm = SS_COUNT[1]/SS_TOT
        FREQ_C_norm = SS_COUNT[2]/SS_TOT

        HELIX_norm = HELIX/SS_TOT
        STRAND_norm = STRAND/SS_TOT
        COIL_norm = COIL/SS_TOT
        TOT_norm = TOT/SS_TOT
        SS_norm = SS_COUNT/SS_TOT

        ### INFORMATION MATRICES archive ###
        INFO_H = np.log((HELIX_norm)/(TOT_norm * FREQ_H_norm))
        INFO_E = np.log((STRAND_norm)/(TOT_norm * FREQ_E_norm))
        INFO_C = np.log((COIL_norm)/(TOT_norm * FREQ_C_norm))

        np.save(OUTDIR + 'H.info', INFO_H)
        np.save(OUTDIR + 'E.info', INFO_E)
        np.save(OUTDIR + 'C.info', INFO_C)

        ### NORMALIZED MATRICES archive ###
        np.save(OUTDIR + 'H.norm_matrix', HELIX_norm)
        np.save(OUTDIR + 'E.norm_matrix', STRAND_norm)
        np.save(OUTDIR + 'C.norm_matrix', COIL_norm)
        np.save(OUTDIR + 'TOT.norm_matrix', TOT_norm)
        np.save(OUTDIR + 'SS.norm_matrix', SS_norm)