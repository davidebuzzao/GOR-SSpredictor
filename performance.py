#!/usr/local/bin/python3
import sys
from sklearn.metrics import multilabel_confusion_matrix
import numpy as np

def compute_similarity(expected, predicted):
    # ss = np.array([[1,0,0], [0,1,0], [0,0,1]])
    # dictionary = {'H':ss[0], 'E':ss[1], '-':ss[2]}
    true = []
    pred = []
    with open(expected) as dssp, open(predicted) as gor:
        for i,j in zip(dssp, gor):
            i = i.rstrip()
            j = j.rstrip()
            for k in range(len(i)):
                true.append(i[k])
                pred.append(j[k])
                # np.append(true, dictionary[i[k]])
                # np.append(pred, dictionary[j[k]])
                #print(dictionary[i[k]], dictionary[j[k]])
    return(true, pred)

if __name__ == '__main__':
    try:
        dssp_file = sys.argv[1]
        predicted_file = sys.argv[2]
    except:
        raise SystemExit
    else:
        x,y = compute_similarity(dssp_file, predicted_file)
        multiclass_conf_mat = multilabel_confusion_matrix(x,y, labels=['H','E','-'])
        print(multiclass_conf_mat)