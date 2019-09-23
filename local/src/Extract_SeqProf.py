#!/usr/local/bin/python3

from sys import argv 
import numpy as np

def pretty_matrix(M):
    for sublist in M:
        print(sublist)

def Extract_SeqProfile(PSSM_FILE):
    PROFILE = []
    with open(PSSM_FILE) as pssm_file:
        for line in pssm_file:
            line = line.rstrip().split()[22:-2]
            PROFILE.append(line)

    for sublist in range(len(PROFILE)): 
        for num in range(20):
            PROFILE[sublist][num] = int(PROFILE[sublist][num])/100
    
    return(PROFILE)

if __name__ == '__main__':
    try:
        PSSM_FILE = argv[1]
        ID = argv[2]
        OUTDIR = argv[3]
    except:
        print('Program Usage: Script Name.py <PSSM FILE> <NAME> <OUTPUT DIRECTORY>')
        raise SystemExit
    else:
        RESIDUE = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        PROFILE = Extract_SeqProfile(PSSM_FILE)
        
        PROFILE_NPY = np.array(PROFILE)
        NAME_FILE = OUTDIR + ID
        np.save(NAME_FILE, PROFILE_NPY)
        #ciao = np.load('/Users/davidebuzzao/Bioinformatica_Bologna/Second_Year/LB2/ModuleB/Project/Results/d1k7ca.npy')
