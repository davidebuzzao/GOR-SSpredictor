#!/usr/local/bin/python3

class GOR:
    
    def __init__(self, window=17):
        import numpy as np
        
        residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        self.window = window

        self.helix = np.zeros((window, len(residues)))
        self.strand = np.zeros((window, len(residues)))
        self.coil = np.zeros((window, len(residues)))
        self.ss_count = np.zeros(3)
    
    def fit(self, profile_set, dssp, padding=True):
        import numpy as np
        ss = ['H', 'E', '-']

        if padding == True: 
            padding = np.zeros((self.window//2,20))
            for instance in range(len(profile_set)):
                #print(instance)
                profile_set[instance] = np.vstack((padding, profile_set[instance], padding))
                dssp[instance] = ' '*8 + dssp[instance] + ' '*8

                i, j, k = 0, self.window//2, self.window
                while k <= len(instance):
                    if dssp[instance][j] not in ss: 
                        i, j, k = i+1, j+1, k+1
                    else:
                        if dssp[instance][j] == 'H':
                            self.helix += profile_set[instance][i:k]
                            self.ss_count[0] += 1
                        elif dssp[instance][j] == 'E':
                            self.strand += profile_set[instance][i:k]
                            self.ss_count[1] += 1
                        else:
                            self.coil += profile_set[instance][i:k]
                            self.ss_count[2] += 1
                        
                        tot += profile_set[i:k]
                        i, j, k = i+1, j+1, k+1


    def __str__(self):
        return '\n{0}\n{1}\n{2}\n'.format(self.helix, self.strand, self.coil)


class PSSM:

    def __init__(self, path_profile, type_profile):
        import numpy as np
        init_profile = []
        
        if type_profile == 'psiblast':
            with open(path_profile) as pssm_file:
                for row in pssm_file:
                    row = row.rstrip().split()[22:-2]
                    init_profile.append(row)
                self.profile = np.array(init_profile)

    def normalize(self):
        residues = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        
        for row in range(len(self.profile)): 
            for col in range(len(residues)):
                self.profile[row][col] = int(self.profile[row][col])/100

    def __str__(self):
        return str(self.profile)
    
    def __repr__(self):
        return str(self)


if __name__ == '__main__':
    import sys
    import os
    import numpy as np
    FILEIN = sys.argv[1]
    PWD = sys.argv[2]

    PROFILES = []
    dssps = []
    with open(FILEIN) as filein:
        for id in filein:
            id = id.rstrip()
            with open(PWD + '/dataset/TrainingSet/DSSP/' + id + '.dssp') as dssp:
                ss = dssp.read().splitlines()[1]
                path = PWD + '/dataset/TrainingSet/PSIBLAST/clean_pssm/' + id + '.clean.pssm'
                PROFILE = PSSM(path, 'psiblast')
                PROFILE.normalize()
                #print(PROFILE)
                PROFILES.append(PROFILE)
                dssps.append(ss)
                #print(np.shape(PROFILE))
    '''
    model = GOR(17)
    model.fit(PROFILES, dssps, True)
    print(model)

    profile = PSSM('/Users/davidebuzzao/Projects/GOR/profile.txt', 'psiblast')
    print(profile)
    profile.normalize()
    '''