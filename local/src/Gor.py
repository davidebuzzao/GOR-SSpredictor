#!/usr/local/bin/python3
import numpy as np, sys
from Tools import Dataset, Pssm, Dssp

class Gor:

    def __init__(self, window=17):
        self.window = window
        self.num_of_updates = 0
        self.residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        self.ss = ['-','H','E']

        tensor = np.zeros((len(self.ss), window, len(self.residues)), dtype=np.float64)
        self.dictionary = {'-': tensor[0], 'H': tensor[1], 'E': tensor[2]}
        self.ss_count = {'-': 0, 'H': 0, 'E': 0}
        self.overall = np.zeros((17,20))
        self.pad = np.zeros((self.window//2,20))

    def train(self, dataset=False, padding=True):
        try: dataset != False
        except: 
            print('Method usage: obj.fit(dataset=X)')
            raise SystemExit
        else:
            for key,val in dataset.items():
                prof = val['profile']
                dssp = val['dssp']

                if padding: 
                    prof = np.vstack((self.pad, prof, self.pad))
                    dssp = ' '*8 + dssp + ' '*8

                i, j, k = 0, self.window//2, self.window
                while k <= len(dssp):
                    if np.sum(prof[i:k]) == 0:
                        i,j,k = i+1, j+1, k+1
                        continue
                    self.dictionary[dssp[j]] += prof[i:k]
                    self.ss_count[dssp[j]] += 1
                    i,j,k = i+1, j+1, k+1

            return(self)

    def normalize(self):
        '''Normalizer of matrices'''
        normalizer = 0
        for index in range(len(self.ss)):
            self.overall += self.dictionary[self.ss[index]]
            normalizer += self.ss_count[self.ss[index]]

        for index in range(len(self.ss)):
            self.dictionary[self.ss[index]] /= normalizer
            self.ss_count[self.ss[index]] /= normalizer

        self.overall /= normalizer
        return(self)

    def information(self):
        '''Information matrices'''
        self.normalize()

        for index in range(len(self.ss)):
            self.dictionary[self.ss[index]] = np.log(self.dictionary[self.ss[index]] / ((self.overall)*(self.ss_count[self.ss[index]])))
        return(self)

    def predict(self, dataset=False, padding=True):
        try: dataset != False
        except: 
            print('Method usage: obj.fit(dataset=X)')
            raise SystemExit
        else:
            zeros = []
            for key,val in dataset.items():
                prof = val['profile']
                if np.sum(prof) == 0:
                    zeros.append(key)

                probabilities = [0,0,0]

                if padding: 
                    prof = np.vstack((self.pad, prof, self.pad))

                seq_pred = ''
                i,j,k = 0, self.window//2, self.window
                while k <= len(prof):
                    for index in range(3):
                        probabilities[index] = np.sum(self.dictionary[self.ss[index]] * prof[i:k])
                    seq_pred += self.ss[probabilities.index(max(probabilities))]
                    i,j,k = i+1, j+1, k+1
                
                val['gor_pred'] = seq_pred
            
            if len(zeros) > 0:
                print('This profiles are all zeros\n', zeros)

            return(dataset)
    
    def save(self, matrices_dir):
        output_names = ['coil_info_matrix', 'helix_info_matrix', 'strand_info_matrix']
        for index in range(len(self.ss)):
            np.save(matrices_dir + output_names[index], self.dictionary[self.ss[index]])
        return(self)

    def load(self, matrices_dir):
        input_names = ['coil_info_matrix', 'helix_info_matrix', 'strand_info_matrix']
        for index in range(len(self.ss)):
            self.dictionary[self.ss[index]] = np.load(matrices_dir + input_names[index] + '.npy')
        return(self)


if __name__ == '__main__':
    try:
        filein = sys.argv[1]
        profiles_dir = sys.argv[2]
        dssps_dir = sys.argv[3]
        matrices_dir = sys.argv[4]
    except:
        print('Program Usage: ')
        raise SystemExit
    else:
        padding = np.zeros((17//2,20))
        model = Gor(17)

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
        model.save(matrices_dir)