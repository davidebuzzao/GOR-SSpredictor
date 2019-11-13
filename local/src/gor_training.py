#!/usr/local/bin/python3
import numpy as np, sys
from Tools import Dataset, Pssm, Dssp
from Gor import Gor

'''
    Input file:
    
        Fold                Dataset IDs                     Dataset.pkl
dataset/cv/fold1/     dataset/cv/fold1/cv1234.id     dataset/cv/fold1/cv1234.pkl
'''

if __name__ == '__main__':
    try:
        filein = sys.argv[1]
    except:
        print('Program Usage: gor_training <cv_train.txt>')
        raise SystemExit
    else:
        with open(filein) as f:
            for line in f:
                items = line.split()
                print(items)
                
                wd = items[0]
                data_id = items[1]
                data_pkl = items[2]


        # Build the dataset from scratch
                prof = Pssm(data_id, setype='trainingset', raw_file=False).parse()
                dict_prof = prof.fetch_dict()
                dssp = Dssp(data_id, setype='trainingset', raw_file=False).parse()
                dict_dssp = dssp.fetch_dict()
                dataset = Dataset(data_id, setype='trainingset').build(profile=dict_prof, dssp=dict_dssp).fetch_dict()

                model = Gor(window=17)\
                        .train(dataset=dataset, padding=True)\
                        .information()\
                        .save(matrices_dir=wd)