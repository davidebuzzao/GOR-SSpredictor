#!/usr/local/bin/python3

'''
Residue-level indexes computed from confusion matrices are not sufficient to properly score a SS prediction method.
Assessing the biological significance of predictions requires comparing predicted and observed segments of residues in the same conformation
• α-Helix segments shorter than 4 residues are meaningless in biology;
• Strands usually involves 2 or more residues (except for isolated β-bridges)
Specific segment-level indexes are then required to analyse these aspects and assess the biological significance of SS predictions.
'''


import pandas as pd
import numpy as np
import argparse

def sov(expected_file, predicted_file):
    secondary_structure = ['E']#, 'E', '-']
    ss_dictionary = {'H':[0,0,0], 'E':[0,0,0], '-':[0,0,0]}
    with open(expected_file) as filein1, open(predicted_file) as filein2:
        for dssp,pred in zip(filein1, filein2):

            for ss in secondary_structure:

                multi_overlaps, length_seq, num_of_ss = sov_parser(dssp, pred, secondary_structure, ss)
                if num_of_ss != 0: 
                    ss_dictionary[ss][0] = sov_scorer(multi_overlaps, length_seq, num_of_ss)
                    ss_dictionary[ss][1] += 1
    for val in ss_dictionary.values():
        val[2] = val[1] / val[2]

    return(ss_dictionary)


# def sov_parser(dssp, pred, secondary_structure, ss):
#     '''
#     This function is intended to compute the Segment OVerlap index.
#     A segment-level index which evaluated the percent average overlap 
#     between predicted and observed segments in the three SS conformations.
#     '''
#     dssp = dssp.rstrip()
#     pred = pred.rstrip()
#     sequences = [dssp, pred]
#     # try:
#     #     len(sequences)%2 == 0
#     # except:
#     #     print("Different lenght between expected and predicted sequence: \n%s\n%s" %(dssp,pred))
#     #     raise SystemExit
#     # else:
#     ss_dictionary = {}
#     num_of_ss = int(dssp.count(ss))
#     multi_overlaps = []
#     skip_number = 0
#     dbl_overl_index = 0
    
#     single_overlap = [0,0,0]
#     indeces = [[0,0],[0,0]] ## start,end
#     i = 0
#     while i <= len(dssp):
#         print(i, skip_number)

#         if skip_number != 0:
#             if indeces[dbl_overl_index][0] == 0: 
#                 indeces[dbl_overl_index][0] = i
#             skip_number -= 1
#             i += 1
#             continue
#         else:
#         # if skip_number != 0:
#         #     #i += indeces[dbl_overl_index][1]
#         #     #print(indeces[dbl_overl_index][1])
#         #     skip_number = 0

#             single_overlap = [0,0,0]
#             indeces = [[0,0],[0,0]] ## start,end

#             ### counting length of S1 and S2 ###
#             for it in range(2):
#                 if sequences[it][i] == ss: 
#                     print(sequences[it][i]) 

#                     if indeces[it][0] == 0: 
#                         indeces[it][0] = i
#                         print(indeces)
#                     for k in range(i, len(dssp)):

#                         if sequences[it][k] != ss: break
#                         else: 
#                             indeces[it][1] = k

#                     single_overlap[it] = (indeces[it][1] - indeces[it][0]) + 1
#                 break
#                 # else: i += 1; print(i); break

#             if (single_overlap[0] == 0 and single_overlap[1] != 0) or (single_overlap[0] != 0 and single_overlap[1] == 0): 
#                 '''
#                 Update i iterable by checking if there could be overalp when 
#                 1 of 2 fragments has not been detected
#                 '''
#                 ### counting overlaps btw S1 and S2 ###
#                 print(indeces)
#                 for j in range(i, i + max(single_overlap[0], single_overlap[1])):
#                     # print(j)
#                     if sequences[0][j] == sequences[1][j]: single_overlap[2] += 1; #print('ciao')

#                 ### if there's overlap ###
#                 if single_overlap[2] != 0:

#                     seq = sequences[single_overlap.index(0)]
#                     seq_index = indeces[single_overlap.index(0)]

#                     if single_overlap.index(0) == 0:
#                         comp_len = single_overlap[1]
#                     else:
#                         comp_len = single_overlap[0]

#                     for k in range(i+1, i + comp_len - 1):
#                         if seq[k] == ss: 
#                             seq_index[0] = k

#                             for z in range(k, len(seq)):
#                                 if seq[z] != ss: break
#                                 else: seq_index[1] = z

#                             single_overlap[single_overlap.index(0)] = (seq_index[1] - seq_index[0]) + 1

#                             break
#                         pass

#                     if abs(indeces[0][1] - indeces[1][1]) > 1:
#                         multi_overlaps.append(single_overlap)
#                         skip_number = single_overlap[2]
#                         dbl_overl_index = indeces.index(max(indeces))
#                         i += int(indeces[dbl_overl_index][0]) 

#                     else:
#                         print(i)
#                         print(single_overlap)
#                         multi_overlaps.append(single_overlap)
#                         i += int(indeces[dbl_overl_index][1]) +1
                    

#                 ### if there's no overlap ###
#                 else:
#                     i += max(single_overlap)

#             elif single_overlap[0] != 0 and single_overlap[1] != 0:
#                 '''
#                 Update i iterable when for sure there's overalp since both of 
#                 the fragments have been detected
#                 '''
#                 for j in range(i, i + max(single_overlap[0], single_overlap[1])):
#                     if sequences[0][j] == sequences[1][j]: single_overlap[2] += 1
#                 multi_overlaps.append(single_overlap)

#                 #i += max(single_overlap[0], single_overlap[1])
            
#             else:
#                 print(i)
#                 print(single_overlap)
#                 '''
#                 Update i iterable when for sure there's no overalp since both of 
#                 the fragments haven't been detected
#                 '''
#                 i += 1
#             if i == 30: break
#     print(multi_overlaps)   
#     return multi_overlaps, len(dssp), num_of_ss


def sov_scorer(multi_overlaps, length_seq, number_of_ss):
    summatory = 0

    for sublist in multi_overlaps:
        observed = sublist[0]
        predicted = sublist[1]
        overlap = sublist[2]
        
        denominator = observed + predicted - overlap
        delta = max([denominator, overlap, observed/2, predicted/2])
        
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
    dataframe = pd.DataFrame(sov_scores)
    print(dataframe)
    #pretty_dictionary(sov_scores)