#!/usr/local/bin/python3
import numpy as np, sys, pickle
from Tools import Dataset, Pssm, Dssp
from Gor import Gor

'''
    Input file:
    
         Fold                Dataset IDs                Dataset.pkl
    dataset/cv/fold1/   dataset/cv/fold1/cv5.id     dataset/cv/fold1/cv5.pkl
'''

if __name__ == '__main__':
    try:
        filein = sys.argv[1]
    except:
        print('Program Usage: python3 svm_encode.py <file.txt>')
        raise SystemExit
    else:
        with open(filein) as f:
            for line in f:
                items = line.split()
                print(items)
                
                wd = items[0]
                data_id = items[1]
                pred_file = items[2]

                prof = Pssm(data_id, setype='trainingset', raw_file=False).parse()
                dict_prof = prof.fetch_dict()
                dssp = Dssp(data_id, setype='trainingset', raw_file=False).parse()
                dict_dssp = dssp.fetch_dict()
                dataset = Dataset(data_id, setype='trainingset').build(profile=dict_prof, dssp=dict_dssp).fetch_dict()

                prediction = Gor(window=17)\
                                .load(matrices_dir=wd)\
                                .predict(dataset=dataset, padding=True)
                
                with open(pred_file, 'wb') as fileout:
                    pickle.dump(prediction, fileout)
                    print('%s: done!' %pred_file)