#!/usr/local/bin/python3

import sys, matplotlib.pyplot as plt, matplotlib.ticker as mtick, numpy as np, pandas as pd
from Tools import Dataset, Dssp, Fasta

residues = ['G','A','V','P','L','I','M','F','W','Y','S','T','C','N','Q','H','D','E','K','R']
structures = ['-','H','E']

## res propensity to fold in SS
def rSS_compo(dataset):
   
    res_dict = dict([(res, [0.0, 0.0, 0.0, 0.0]) for res in residues]) ## Dictionary format = {res: [-, H, E, total_res]}
    str_dict = dict([(struct, [0.0, 0.0]) for struct in structures]) ## Dictionary format = {SS: [num_SS, freq_SS]}

    for val in dataset.values():
        fasta = val['fasta']
        dssp = val['dssp']
        try: (len(fasta) + len(dssp)) % 2 == 0
        except:
            print("Different lenght between fasta and dssp: \n%s\n%s" %(fasta, dssp))
            raise SystemExit
        else:
            for i in range(len(structures)):
                str_dict[structures[i]][0] += dssp.count(structures[i])
            for ch in range(len(fasta)):
                if fasta[ch] not in res_dict: continue 
                res_dict[fasta[ch]][structures.index(dssp[ch])] += 1
                res_dict[fasta[ch]][-1] += 1
            
    return(res_dict, str_dict)

def rSS_stat(res_dict,str_dict):
    percentages = dict([(res, [0.0, 0.0, 0.0, 0.0]) for res in residues])

    ss_total = 0
    for val in str_dict.values():
        ss_total += val[0]
    for key,val in str_dict.items():
        str_dict[key][1] = round((str_dict[key][0]/ss_total)*100,2)

    for key,val in res_dict.items():
        for i in range(len(structures)):
            percentages[key][i] = round((val[i]/str_dict[structures[i]][0])*100,2)
        percentages[key][-1] = round((val[-1]/ss_total)*100,2)
    
    dataframe = pd.DataFrame(percentages).T
    dataframe.columns = ['Coil', 'Helix', 'Strand', 'Residue']

    return(dataframe)

def rSS_barplot(dataframe, setype, outdir=False, RGB='Greys_r'):
    color = RGB
    df = dataframe.plot(kind='bar', rot=0, colormap=color, edgecolor='black') # use a method for dataframe
    df.set(xlabel="Residues", ylabel="Residue Frequency (%)", title='Residue Composition')
    df.yaxis.set_major_formatter(mtick.PercentFormatter())

    if outdir:
        plt.savefig(outdir + '/Rcompo_' + setype + '.png', box_inches='tight')
    else:
        plt.show()

def ss_piechart(str_dict, setype, outdir=False, RGB='Greys_r'):
    color = RGB
    sizes = [str_dict['-'][1], str_dict['H'][1], str_dict['E'][1]]
    labels = ['Coil', 'Helix', 'Strand']
    
    fig1, ax1 = plt.subplots()
    theme = plt.get_cmap(color)
    ax1.set_prop_cycle("color", [theme(1. * i / len(sizes))
                                    for i in range(len(sizes))])

    ax1.pie(sizes,\
            #explode=(0.05,0.05,0.05),
            labeldistance=1.2,
            #autopct='%1.1f%%',
            startangle=70,
            wedgeprops={"edgecolor":"k", 'linewidth': 1, 'linestyle': 'solid', 'antialiased': True})
    
    plt.legend(loc='upper left',
            labels=['%s, %1.1f%%' % (
                    l, s) for l, s in zip(labels, sizes)],
            prop={'size': 9},
            bbox_to_anchor=(0.0, 1),
            bbox_transform=fig1.transFigure)

    #draw circle
    centre_circle = plt.Circle((0,0),0.50,fc='white', edgecolor='black')
    fig = plt.gcf()
    fig.gca().add_artist(centre_circle)

    # Equal aspect ratio ensures that pie is drawn as a circle
    ax1.axis('equal')
    plt.title('SS Composition', fontname="Times New Roman", fontweight="bold")
    
    if outdir:
        plt.savefig(outdir + '/SScompo_' + setype + '.png', box_inches='tight')
    else:
        plt.show()


if __name__ == '__main__':
    try:
        setype = sys.argv[1]
    except:
        print('Program Usage: res_compo <setype (trainingset, blindset)>')
        raise SystemExit
    else:
        wd = 'dataset/' + setype + '/'
        if setype == 'trainingset':
            data_id = wd + 'ts.id'
        elif setype == 'blindset':
            data_id = wd + 'bs.id'

        # Build the dataset from scratch
        fasta = Fasta(data_id, setype=setype).parse(fastatype='singleline').fetch_dict()
        dssp = Dssp(data_id, setype=setype, raw_file=False).parse().fetch_dict()
        dataset = Dataset(data_id, setype=setype).build(fasta=fasta, dssp=dssp).fetch_dict()

        res_dict, str_dict = rSS_compo(dataset)
        dataframe = rSS_stat(res_dict,str_dict)

        rSS_barplot(dataframe, setype)
        ss_piechart(str_dict, setype)