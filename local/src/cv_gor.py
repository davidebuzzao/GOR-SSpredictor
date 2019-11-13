#!/usr/local/bin/python3

from gor import Gor
from Database import Database, Pssm, Dssp
from performance import compute_similarity, print_performance
from sov import sov, sov_parser, sov_scorer
from sklearn.metrics import multilabel_confusion_matrix
import argparse, numpy as np, pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('CV', help = 'This is the Cross Validation .txt file with <relative_path_to_dir> <train_set.id> <test_set.id>')
    args = parser.parse_args()
    ss = ['H', 'E', '-']

    file_input = []
    with open(args.CV) as filein:
        for line in filein:
            file_input.append(line.rstrip().split())

    cv_list = ['f1', 'f2', 'f3', 'f4', 'f5']
    scores = ['MCC', 'ACC', 'TPR', 'PPV', 'FPR', 'NPV', 'SOV']
    cv_dictionary = dict([key1, dict([key2, np.zeros(len(ss))] for key2 in scores)] for key1 in cv_list)
    q3_scores = np.zeros(len(cv_list))

    w = 17
    fold = 0
    for cv in file_input:
        training_set = cv[0] + cv[1]
        validation_set = cv[0] + cv[2]

        with open(training_set) as ts, open(validation_set) as vs:
            
            '''
            Training Step
            '''
            train_prof = Database(datatype='psiblast', raw_file=False, window=w, path='./psiblast/bin/')
            train_dssp = Database(datatype='dssp', raw_file=False, window=w, path='./dssp/ss/')
            train_prof.build_dataset(training_set); train_prof_dtset = train_prof.dataset
            train_dssp.build_dataset(training_set); train_dssp_dtset = train_dssp.dataset
            model = Gor(window=w) 
            model.fit(train_prof_dtset, train_dssp_dtset)
            model.information()
            
            '''
            Testing Step
            '''
            test_prof = Database(datatype='psiblast', raw_file=False, window=w, path='./psiblast/bin/')
            test_prof.build_dataset(validation_set); test_prof_dtset = test_prof.dataset
            predictions = model.predict(test_prof_dtset, padd=True)

            expectations = []
            for id in vs:
                id = id.rstrip()
                name = open('./dssp/ss/' + id + '.dssp')
                expectations.append(name.read().splitlines()[1])

            '''
            Performance Step
            1. Confusion Matrix and related scores
            2. Segment OVerlapping score
            '''
            x,y = compute_similarity(expectations, predictions)
            multiclass_conf_mat = multilabel_confusion_matrix(x,y, labels=['H','E','-'])
            total = np.sum(multiclass_conf_mat[0])

            true_predicted = 0

            for index in range(3):
                true_predicted += multiclass_conf_mat[index][1][1]
                MCC, ACC, TPR, PPV, FPR, NPV = print_performance(multiclass_conf_mat[index])
                score = [MCC, ACC, TPR, PPV, FPR, NPV]
                
                scr = 0
                while scr < len(score):
                    cv_dictionary[cv_list[fold]][scores[scr]][index] = round(score[scr]*100,2)
                    scr += 1

            q3 = true_predicted/total
            q3_scores[fold] = q3

            sov_scores = sov(expectations, predictions, file=False)
            final_sov = np.zeros(len(ss))
            for i in range(len(ss)):
                final_sov[i] = sov_scores[ss[i]][2]

            cv_dictionary[cv_list[fold]]['SOV'] += np.round(final_sov,2)
            fold += 1

    ### Find the way to print better c !!! ###
    d = pd.DataFrame(cv_dictionary)

    print(d)
    cv_average_ss = dict([key, np.zeros(len(ss))] for key in scores)
    cv_average = dict([key, np.zeros(1)] for key in scores)

    divisor1 = len(cv_list)
    
    for i in range(len(ss)):
        scr = 0
        while scr < len(scores):
            for j in cv_list:
                cv_average_ss[scores[scr]][i] += cv_dictionary[j][scores[scr]][i]
            cv_average_ss[scores[scr]][i] /= divisor1
            scr += 1
    d1 = pd.DataFrame(cv_average_ss, index=['H','E','C'])
    print('\n', d1, '\n')
    
    divisor2 = len(ss)
    for i in range(len(scores)):
        cv_average[scores[i]] += np.sum(cv_average_ss[scores[i]])
        cv_average[scores[i]] /= divisor2
    d2 = pd.DataFrame(cv_average, index=['Average values'])
    print(d2, '\n')

    d.to_pickle('./dataset/TrainingSet/cv/folds_ss_performance.pkl')
    d1.to_pickle('./dataset/TrainingSet/cv/average_ss_performance.pkl')
    d2.to_pickle('./dataset/TrainingSet/cv/average_performance.pkl')