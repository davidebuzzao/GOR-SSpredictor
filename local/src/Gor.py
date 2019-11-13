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

    def __str__(self):
        return '\n{0}\n{1}\n{2}\n'.format(self.dictionary[0],self.dictionary[1],self.dictionary[2])

    def train(self, dataset=False, padding=True):
        try: dataset != False
        except: 
            print('Method usage: obj.fit(dataset=X)')
            raise SystemExit
        else:
            for key,val in dataset.items():
                prof = val['profile']
                dssp = val['dssp']

                i, j, k = 0, self.window//2, self.window
                while k <= len(dssp):
                    self.dictionary[dssp[j]] += prof[i:k]
                    self.ss_count[dssp[j]] += 1
                    i,j,k = i+1, j+1, k+1

            return(self)

    def normalize(self):
        '''Normalizer of matrices'''
        normalizer = 0
        for index in range(3):
            self.overall += self.dictionary[self.ss[index]]
            normalizer += self.ss_count[self.ss[index]]
        
        for index in range(3):
            self.dictionary[self.ss[index]] /= normalizer
            self.ss_count[self.ss[index]] /= normalizer
        
        self.overall /= normalizer
        return(self)

    def information(self):
        '''Information matrices'''
        self.normalize()

        for index in range(3):
            self.dictionary[self.ss[index]] = np.log(self.dictionary[self.ss[index]] / ((self.overall)*(self.ss_count[self.ss[index]])))
        return(self)

    def predict(self, dataset=False, padding=True):
        try: dataset != False
        except: 
            print('Method usage: obj.fit(dataset=X)')
            raise SystemExit
        else:
            for key,val in dataset.items():
                prof = val['profile']
                probabilities = [0,0,0]

                if padding: 
                    pad = np.zeros((self.window//2,20))
                    prof = np.vstack((pad, prof, pad))

                seq_pred = ''
                i,j,k = 0, self.window//2, self.window
                while k <= len(prof):
                    for index in range(3):
                        probabilities[index] = np.sum(self.dictionary[self.ss[index]] * prof[i:k])
                    seq_pred += self.ss[probabilities.index(max(probabilities))]
                    i,j,k = i+1, j+1, k+1
                
                val['gor_pred'] = seq_pred

            return(dataset)
    
    def save(self, matrices_dir):
        output_names = ['coil_info_matrix', 'helix_info_matrix', 'strand_info_matrix']
        for index in range(3):
            np.save(matrices_dir + output_names[index], self.dictionary[self.ss[index]])
        return(self)

    def load(self, matrices_dir):
        input_names = ['coil_info_matrix', 'helix_info_matrix', 'strand_info_matrix']
        for index in range(3):
            self.dictionary[self.ss[index]] = np.load(matrices_dir + input_names[index] + '.npy')
        return(self)

# class PSSM:

#     num_of_instances = 0

#     def __init__(self, path_profile, type_profile):
#         import numpy as np
#         init_profile = []
        
#         if type_profile == 'psiblast':
#             with open(path_profile) as pssm_file:
#                 for row in pssm_file:
#                     row = row.rstrip().split()[22:-2]

#                     init_profile.append(row)
#                 self.prof = np.array(init_profile, dtype=np.float64)
        
#         PSSM.num_of_instances += 1

#     def normalize(self):
#         residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        
#         for row in range(len(self.prof)): 
#             for col in range(len(residues)):
#                 self.prof[row][col] = float(self.prof[row][col])/100

#     def shape(self):
#         row = 1
#         column = len(self.prof)
#         return((row, column))

#     # def __str__(self):
#     #     return str(self.prof)
    
#     def __repr__(self):
#         return str(self.prof)

# def prof_parse(path):
#     residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
#     init_profile=[]
#     with open(path) as pssm_file:
#         for row in pssm_file:
#             row = row.rstrip().split()[22:-2]
#             init_profile.append(row)
#         prof = np.array(init_profile, dtype=np.float64)
        
#     for row in range(len(prof)): 
#         for col in range(len(residues)):
#             prof[row][col] = float(prof[row][col])/100
#     return(prof)

# if __name__ == '__main__':
#     try:
#         filein = sys.argv[1]
#         profiles_dir = sys.argv[2]
#         dssps_dir = sys.argv[3]
#         matrices_dir = sys.argv[4]
#     except:
#         print('Program Usage: ')
#         raise SystemExit
#     else:
#         padding = np.zeros((17//2,20))
#         model = Gor(17)

#         with open(filein) as f:
#             for id in f:
#                 id = id.rstrip()
#                 with open(dssps_dir + id + '.dssp') as dssp_file:
#                     path = profiles_dir + id + '.clean.pssm'
#                     dssp = dssp_file.read().splitlines()[1]
#                     profile = prof_parse(path)
#                     out = np.vstack((padding, profile, padding))
#                     model.fit(profile, dssp, True)

#         model.information()
#         #np.set_printoptions(threshold=np.inf)
#         model.save(matrices_dir)