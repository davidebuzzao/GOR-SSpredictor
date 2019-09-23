#!/usr/local/bin/python3

import sys 
import matplotlib.pyplot as plt

def ss_compo(dssp_file):
    H = E = C = 0
    for line in dssp_file:
        clean_line = str(line.rstrip())
        for ch in clean_line:
            if ch == '-':
                C += 1
            else:
                if ch == 'H':
                    H += 1
                else: 
                    E += 1
    return H, E, C

if __name__ == '__main__':
    try:
        DSSP_FILE = sys.argv[1] #
        OUT_DIR = sys.argv[2]
    except:
        print('Program Usage: text.py <DSSP_FILE> <OUT_DIR>')
        raise SystemExit
    else:
        with open(DSSP_FILE) as FILE_IN:
            H, E, C = ss_compo(FILE_IN)
            TOT = H + E + C
            print(TOT)
        '''
        # Pie chart, where the slices will be ordered and plotted counter-clockwise:
        labels = 'Helix', 'Strand', 'Coil'
        sizes = [H/TOT, E/TOT, C/TOT]

        #colors
        colors = ['#ff9999','#66b3ff','#99ff99']
        
        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, colors = colors, labels=labels, autopct='%1.1f%%', startangle=90)
        
        #draw circle
        centre_circle = plt.Circle((0,0),0.70,fc='white')
        fig = plt.gcf()
        fig.gca().add_artist(centre_circle)

        # Equal aspect ratio ensures that pie is drawn as a circle
        ax1.axis('equal')
        plt.title('SS composition', fontname="Times New Roman", fontweight="bold")
        plt.savefig(OUT_DIR + 'SS_composition', box_inches='tight')
        '''