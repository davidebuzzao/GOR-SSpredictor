#!/usr/local/bin/python3
import numpy as np, pandas as pd, sys, pickle
from Tools import Dataset, Pssm, Dssp
from Gor import Gor
from sov import sov, sov_parser, sov_scorer
from performance import compute_similarity, print_performance
from sklearn.metrics import multilabel_confusion_matrix

'''
    Input file:
    
    Fold          Dataset IDs                  Dataset.pkl
    fold1     dataset/cv/fold1/cv5.id    dataset/cv/fold1/cv5.pkl     
'''

if __name__ == '__main__':
    try:
        filein = sys.argv[1]
    except:
        print('Program Usage: gor_performance <cv_pred.txt>')
        raise SystemExit
    else:
        ss = ['-', 'H', 'E']
        cv_list = ['fold1', 'fold2', 'fold3', 'fold4', 'fold5']
        scores = ['MCC', 'ACC', 'TPR', 'PPV', 'SOV']
        cv_dictionary = dict([key1, dict([key2, np.zeros(len(ss))] for key2 in scores)] for key1 in cv_list)
        q3_scores = np.zeros(len(cv_list))

        with open(filein) as f:
            for line in f:
                items = line.split()
                print(items)
                
                fold = items[0]
                data_id = items[1]
                with open(items[2], 'rb') as pred_file:
                    dictionary = pickle.load(pred_file)

                '''
                Performance Steps:
                1. Confusion Matrix and related scores
                2. Segment OVerlapping score
                '''
                prediction = []
                expectation = []
                
                for key in dictionary:
                    expectation.append(dictionary[key]['dssp'])
                    prediction.append(dictionary[key]['gor_pred'])

                x,y = compute_similarity(expectation, prediction)
                multiclass_conf_mat = multilabel_confusion_matrix(x,y, labels=['-','H','E'])
                total = np.sum(multiclass_conf_mat[0])
                #print(multiclass_conf_mat)

                true_predicted = 0

                for index in range(3):
                    true_predicted += multiclass_conf_mat[index][1][1]
                    MCC, ACC, TPR, PPV = print_performance(multiclass_conf_mat[index])
                    score = [MCC, ACC, TPR, PPV]
                    
                    scr = 0
                    while scr < len(score):
                        cv_dictionary[fold][scores[scr]][index] = round(score[scr]*100,2)
                        scr += 1

                q3 = true_predicted/total
                q3_scores[cv_list.index(fold)] = round(q3 * 100, 2)

                sov_scores = sov(expectation, prediction)
                final_sov = np.zeros(len(ss))
                for i in range(len(ss)):
                    final_sov[i] = sov_scores[ss[i]][2]
                
                cv_dictionary[fold]['SOV'] += np.round(final_sov,2)
        
        d = pd.DataFrame(cv_dictionary)
        print('\n',d)

        cv_average_ss = dict([key, np.zeros(len(ss))] for key in scores)
        cv_stderr_ss = dict([key, np.zeros(len(ss))] for key in scores)
        cv_average = dict([key, np.zeros(1)] for key in scores)

        divisor1 = len(cv_list)

        for i in range(len(ss)):
            scr = 0
            while scr < len(scores):
                for j in cv_list:
                    cv_average_ss[scores[scr]][i] += cv_dictionary[j][scores[scr]][i]

                cv_average_ss[scores[scr]][i] /= divisor1
                scr += 1
        
        for i in range(len(ss)):
            scr = 0
            while scr < len(scores):
                for j in cv_list:
                    cv_stderr_ss[scores[scr]][i] += np.power(cv_dictionary[j][scores[scr]][i] - cv_average_ss[scores[scr]][i], 2)

                cv_stderr_ss[scores[scr]][i] = (np.sqrt(cv_stderr_ss[scores[scr]][i]/(len(cv_list)-1)))/np.sqrt(len(cv_list))
                scr += 1

        d1 = pd.DataFrame(cv_average_ss, index=['C','H','E'])
        d2 = pd.DataFrame(cv_stderr_ss, index=['C','H','E'])
        print('\n', d1, '\n')
        print('\n', d2, '\n')

        divisor2 = len(ss)
        for i in range(len(scores)):
            cv_average[scores[i]] += np.sum(cv_average_ss[scores[i]])
            cv_average[scores[i]] /= divisor2

        d3 = pd.DataFrame(cv_average, index=['Average values'])
        print(d3, '\n')

        print('Q3 scores:', q3_scores)
        