#!/usr/local/bin/python3
import numpy as np
import sys
import os

class GOR:
    num_of_updates = 0
    residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    ss = ['H', 'E', '-']

    def __init__(self, window=17):
        self.window = window

        tensor = np.zeros((len(GOR.ss), window, len(GOR.residues)), dtype=np.float64)
        self.dictionary = {'H': tensor[0], 'E': tensor[1], '-': tensor[2]}
        self.ss_count = {'H': 0, 'E': 0, '-': 0}
        self.overall = np.zeros((17,20))

    def __str__(self):
        return '\n{0}\n{1}\n{2}\n'.format(self.dictionary[0],self.dictionary[1],self.dictionary[2])

    def fit(self, profile, dssp, padding=True):
        if padding == True: 
            padding = np.zeros((self.window//2,20))
            profile = np.vstack((padding, profile, padding))
            dssp = ' '*8 + dssp + ' '*8

            i, j, k = 0, self.window//2, self.window
            while k <= len(dssp):
                self.dictionary[dssp[j]] += profile[i:k]
                self.ss_count[dssp[j]] += 1
                i, j, k = i+1, j+1, k+1
        
        GOR.num_of_updates += 1

    def normalize(self):
        '''Normalizr matrices'''
        normalizer = 0
        for index in range(3):
            self.overall += self.dictionary[GOR.ss[index]]
            normalizer += self.ss_count[GOR.ss[index]]
        
        for index in range(3):
            self.dictionary[GOR.ss[index]] = np.log(self.dictionary[GOR.ss[index]]/normalizer)
            self.ss_count[GOR.ss[index]]/normalizer
        self.overall/normalizer

    def information(self):
        '''Information matrices'''
        self.normalize()
        for index in range(3):
            self.dictionary[GOR.ss[index]] = np.log(self.dictionary[GOR.ss[index]] / ((self.overall)*(self.ss_count[GOR.ss[index]])))

    def predict(self, profile, padding=True):
        predicted_sequence = ''
        probabilities = [0,0,0]

        if padding == True: 
            padding = np.zeros((self.window//2,20))
            profile = np.vstack((padding, profile, padding))

        i, j, k = 0, 8, 17
        while k <= len(profile):
            for index in range(3):
                probabilities[index] = np.sum(self.dictionary[GOR.ss[index]] * profile[i:k])
            
            predicted_sequence += GOR.ss[probabilities.index([max(probabilities)])]
            
            i, j, k = i+1, j+1, k+1
        return(predicted_sequence)
    
    def save(self, matrices_dir):
        output_names = ['helix_info_matrix', 'strand_info_matrix', 'coil_info_matrix']
        for index in range(3):
            np.save(matrices_dir + output_names[index], self.dictionary[GOR.ss[index]])

    def load(self, matrices_dir):
        input_names = ['helix_info_matrix', 'strand_info_matrix', 'coil_info_matrix']
        for index in range(3):
            self.dictionary[GOR.ss[index]] = np.load(matrices_dir + input_names[index] + '.npy')


class PSSM:

    num_of_instances = 0

    def __init__(self, path_profile, type_profile):
        import numpy as np
        init_profile = []
        
        if type_profile == 'psiblast':
            with open(path_profile) as pssm_file:
                for row in pssm_file:
                    row = row.rstrip().split()[22:-2]

                    init_profile.append(row)
                self.prof = np.array(init_profile, dtype=np.float64)
        
        PSSM.num_of_instances += 1

    def normalize(self):
        residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        
        for row in range(len(self.prof)): 
            for col in range(len(residues)):
                self.prof[row][col] = float(self.prof[row][col])/100

    def shape(self):
        row = 1
        column = len(self.prof)
        return((row, column))

    # def __str__(self):
    #     return str(self.prof)
    
    def __repr__(self):
        return str(self.prof)

def prof_parse(path):
    residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    init_profile=[]
    with open(path) as pssm_file:
        for row in pssm_file:
            row = row.rstrip().split()[22:-2]
            init_profile.append(row)
        prof = np.array(init_profile, dtype=np.float64)
        
    for row in range(len(prof)): 
        for col in range(len(residues)):
            prof[row][col] = float(prof[row][col])/100
    return(prof)


if __name__ == '__main__':
    try:
        filein = sys.argv[1]
        profiles_dir = sys.argv[2]
        dssps_dir = sys.argv[3]
        matrices_dir = sys.argv[4]
    except:
        raise SystemExit
    else:
        padding = np.zeros((17//2,20))
        model = GOR(17)

        with open(filein) as f:
            for id in f:
                id = id.rstrip()
                with open(dssps_dir + id + '.dssp') as dssp_file:
                    path = profiles_dir + id + '.clean.pssm'
                    dssp = dssp_file.read().splitlines()[1]
                    profile = prof_parse(path)
                    out = np.vstack((padding, profile, padding))
                    model.fit(profile, dssp, True)

        model.information()
        #np.set_printoptions(threshold=np.inf)
        model.save(matrices_dir)