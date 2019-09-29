#!/usr/local/bin/python3

'''
Residue-level indexes computed from confusion matrices are not sufficient to properly score a SS prediction method.
Assessing the biological significance of predictions requires comparing predicted and observed segments of residues in the same conformation
• α-Helix segments shorter than 4 residues are meaningless in biology;
• Strands usually involves 2 or more residues (except for isolated β-bridges)
Specific segment-level indexes are then required to analyse these aspects and assess the biological significance of SS predictions.
'''

from sklearn.metrics import multilabel_confusion_matrix
import pandas as pd
import numpy as np
import argparse

def pretty_dictionary(D):
    for key,val in D.items():
        print(key, " = ", val)
    print()

def ss_compo(expected_file):
    with open(expected_file) as filein1:
        ss = ['H', 'E', '-']
        dictionary={'H':0, 'E':0, '-':0}
        for line in filein1:
            line = str(line.rstrip())
            
            for ch in line:
                dictionary[ch] += 1
    
    return (dictionary)


def sov(expected_file, predicted_file):
    ss_composition = ss_compo(expected_file)
    secondary_structure = ['H', 'E', '-']
    ss_dictionary = {'H':0, 'E':0, '-':0}
    
    for ss in secondary_structure:
        multi_overlaps, length_seq = sov_parser(expected_file, predicted_file, secondary_structure, ss)
        ss_dictionary[ss]= sov_scorer(multi_overlaps, length_seq, ss_composition[ss])
    
    return(ss_dictionary)


def sov_parser(expected_file, predicted_file, secondary_structure, ss):
    '''
    This function is intended to compute the Segment OVerlap index.
    A segment-level index which evaluated the percent average overlap 
    between predicted and observed segments in the three SS conformations.
    '''

    with open(expected_file) as filein1, open(predicted_file) as filein2:

        for dssp,pred in zip(filein1, filein2):
            dssp = dssp.rstrip()
            pred = pred.rstrip()
            sequences = [dssp, pred]
            try:
                len(sequences)%2 == 0
            except:
                print("Different lenght between expected and predicted sequence: \n%s\n%s" %(dssp,pred))
                raise SystemExit
            else:
                multi_overlaps = []
                skip_number = 0
                double_overlapping_index = None
                i = 0
                while i <= len(dssp):
                    single_overlap = [0,0,0]
                    indeces = [[0,0],[0,0]] ## start,end
                    
                    if skip_number != 0:
                        if indeces[double_overlapping_index][0] == 0: 
                            indeces[double_overlapping_index][0] = i
                        skip_number -= 1
                        i += 1
                        continue

                    else:
                        ### counting length of S1 and S2 ###
                        for it in range(2):
                            if sequences[it][i] == ss: 
                                if indeces[it][0] == 0: indeces[it][0] = i
                                for k in range(i, len(sequences[it])):
                                    if sequences[it][k] != ss: break
                                    else: indeces[it][1] = k

                                single_overlap[it] = (indeces[it][1] - indeces[it][0]) + 1
                        
                        if single_overlap[0] == 0 or single_overlap[1] == 0: 
                            '''
                            Update i iterable by checking if there could be overalp when 
                            1 of 2 fragments has not been detected
                            '''
                            ### counting overlaps btw S1 and S2 ###
                            for j in range(i, i + max(single_overlap[0], single_overlap[1])):
                                if sequences[0][j] == sequences[1][j]: single_overlap[2] += 1

                            ### if there's overlap ###
                            if single_overlap[2] != 0:
                                seq = sequences[single_overlap.index(0)]
                                seq_index = indeces[single_overlap.index(0)]
                                
                                if single_overlap.index(0) == 0:
                                    comp_len = single_overlap[1]
                                else:
                                    comp_len = single_overlap[0]
                                
                                for k in range(i+1, i + comp_len - 1):
                                    if seq[k] == ss:
                                        seq_index[0] = k
                                        for z in range(k, len(seq)):
                                            if seq[z] != ss: break
                                            else: seq_index[1] = z
                                        single_overlap[single_overlap.index(0)] = (seq_index[1] - seq_index[0]) + 1
                                        break
                                multi_overlaps.append(single_overlap)
                                skip_number = single_overlap[2]
                                double_overalpping_index = indeces.index(max(indeces))
                                i += double_overlapping_index[0]    ## In order to go ahead, take the start index of 
                                                                    ## the furthest fragment from the starting point.
                                                                    ## Update the skip_number
                                                                    ## ----HHHHH---
                                                                    ## --HHHH--HHHH
                                                                    ## maybe by using lists of precomputed available fragments and pop them every time?!

                            ### if there's no overlap ###
                            else:
                                i += max(single_overlap)


                        elif single_overlap[0] != 0 and single_overlap[1] != 0:
                            '''
                            Update i iterable when for sure there's overalp since both of 
                            the fragments have been detected
                            '''
                            for j in range(i, i + max(single_overlap[0], single_overlap[1])):
                                if sequences[0][j] == sequences[1][j]: single_overlap[2] += 1
                            multi_overlaps.append(single_overlap)

                            i += max(single_overlap[0], single_overlap[1])
                        
                        else:
                            '''
                            Update i iterable when for sure there's no overalp since both of 
                            the fragments haven't been detected
                            '''
                            i += 1
        
        return multi_overlaps, len(dssp)


def sov_scorer(multi_overlaps, length_seq, number_of_ss):
    summatory = 0
    
    for sublist in multi_overlaps:
        observed = sublist[0]
        predicted = sublist[1]
        overlap = sublist[2]
        
        denominator = observed + predicted - overlap
        delta = max([denominator, overlap, length_seq/2])
        
        summatory += (overlap + delta) / (denominator * overlap[0])
    
    sov_score = (100 * summatory) / number_of_ss
    
    return(sov_score)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('DSSP', \
        help='This is the dssp file with concatenating BlindSet secondary structure')
    parser.add_argument('SSpredictions', \
        help='This is the file with concatenating BlindSet secondary structure predicted via GOR')
    args = parser.parse_args()    
    
    sov_scores = sov(args.DSSP, args.SSpredictions)
    pretty_dictionary(sov_scores)