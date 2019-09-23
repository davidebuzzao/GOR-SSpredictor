#!/usr/local/bin/python3

import sys

Norm_Acc={
            "A" :106.0,  "B" :160.0, \
            "C" :135.0,  "D" :163.0,  "E" :194.0, \
            "F" :197.0,  "G" : 84.0,  "H" :184.0, \
            "I" :169.0,  "K" :205.0,  "L" :164.0, \
            "M" :188.0,  "N" :157.0,  "P" :136.0, \
            "Q" :198.0,  "R" :248.0,  "S" :130.0, \
            "T" :142.0,  "V" :142.0,  "W" :227.0, \
            "X" :180.0,  "Y" :222.0,  "Z" :196.0
            }

def pretty_matrix(output_file, M):
    for element in M:
        output_file.write('\t'.join([str(i) for i in element]) + '\n')


def parse_dssp(dssp_file, ch, surface):   
    
    dssp = []
    dssp2 = []
    str_length = 0
    c = 0

    for line in dssp_file:
        if line.find('  #  RESIDUE') == 0: ## when --# Residue appears in line
            c = 1
            continue
        if c == 0: continue 
        if line[13] == '!': continue
        if line[11] == ch or ch == '_':
            r = line[13].upper()
            c = line[11].upper()
            str_length += 1
            pos = line[5:10].strip()
            ss = line[16]
            if ss == ' ': ss = 'C' ## for coils there's empty space
            acc = float(line[35:38])
            phi = float(line[103:109])
            psi = float(line[109:115])
            racc = round(min(acc/Norm_Acc[r], 1.0), 4) #relative accessibility
            v = [r, pos, c, ss, acc, racc, phi, psi]
            v2 = [r, pos, c, acc]
            dssp2.append(v2)
            dssp.append(v)
    
    if surface == 'NO':
        return(dssp, str_length)

    if surface == 'YES':
        return(dssp2)

def get_SS(fileout, dssp_id_ch, id):
    H = ['H','G','I']
    E = ['B','E']
    
    SEQ = ''
    for sublist in DSSP_ID_CH:
        SS = sublist[3]
        if SS in H: 
            SEQ += 'H'
        
        elif SS in E:
            SEQ += 'E'

        else:
            SEQ += '-' 

    fileout.write('>' + id + '\n' + SEQ)


def find_surface(dssp_monomer, dssp_complex, chain, surface):
    DSSP_MONOMER = parse_dssp(dssp_monomer, chain, surface)
    DSSP_COMPLEX = parse_dssp(dssp_complex, chain, surface)
    INTERACT = []
    for i, j in zip(DSSP_MONOMER, DSSP_COMPLEX):
        if float(i[3]) == float(j[3]): continue
        else: 
            i.append(float(i[3]) - float(j[3]))
            INTERACT.append(i)
    
    return(INTERACT)

if __name__ == '__main__':
    try:
        ID_CHAIN = sys.argv[1]
        FASTA = sys.argv[2]
        PWD = sys.argv[3]
        FIND_SURFACE = sys.argv[4]
    except:
        print('Program Usage: text.py <ID_CHAIN> <ID_FASTA> <CHAIN> <YES/NO FIND_SURFACE>')
        raise SystemExit
    else:
        ID = ID_CHAIN.split('_')[0]
        CHAIN = ID_CHAIN.split('_')[1]
        FILEIN = PWD + 'raw/' + ID + '.dssp'

        with open(FILEIN) as DSSP:
        
            if FIND_SURFACE == 'YES':
                try:
                    DSSP_COMPLEX = sys.argv[5]
                except:
                    print('Program Usage: text.py <ID_CHAIN> <CHAIN> <YES for FIND_SURFACE> <DSSP_COMPLEX>')
                    INTERACT = find_surface(ID, DSSP_COMPLEX, CHAIN, FIND_SURFACE)
                    pretty_matrix(INTERACT)
            else:
                OUTDIR1 = PWD + 'clean/'
                OUTDIR2 = PWD + 'ss/'
                FILEOUT1 = OUTDIR1 + ID_CHAIN + '.clean.dssp'
                FILEOUT2 = OUTDIR2 + ID_CHAIN + '.ss'
                
                DSSP_ID_CH, STRUCTURE_LENGTH = parse_dssp(DSSP, CHAIN, FIND_SURFACE)
                with open(FASTA) as SEQUENCE_FILEIN:
                    SEQUENCE = SEQUENCE_FILEIN.readlines()[1].rstrip()
                    if len(SEQUENCE) == STRUCTURE_LENGTH:
                        with open(FILEOUT1, 'w') as CLEAN_DSSP_OUT, open(FILEOUT2, 'w') as SS_OUT:
                            pretty_matrix(CLEAN_DSSP_OUT, DSSP_ID_CH)
                            get_SS(SS_OUT, DSSP_ID_CH, ID_CHAIN)
                    else:
                        DIFFERENCE = len(SEQUENCE)-STRUCTURE_LENGTH
                        FILEOUT3 =  PWD + 'NotDefined_DSSP.list'
                        with open(FILEOUT1, 'w') as CLEAN_DSSP_OUT, open(FILEOUT3, 'a') as DIFF_LENGTH:
                            pretty_matrix(CLEAN_DSSP_OUT, DSSP_ID_CH)
                            DIFF_LENGTH.write(ID_CHAIN + '\t' + str(DIFFERENCE) + '\n')