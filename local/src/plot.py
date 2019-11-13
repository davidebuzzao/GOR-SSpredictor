#!/usr/local/bin/python3

import argparse 
import pandas as pd
import seaborn as sns
import pylab as plt

def change_width(ax, new_value) :
    for patch in ax.patches :
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)


if __name__ == '__main__':

    scores = ['MCC', 'ACC', 'TPR', 'PPV', 'FPR', 'NPV', 'SOV']
    cv_list = ['f1', 'f2', 'f3', 'f4', 'f5']
    ss = ['H', 'E', 'C']

    unpickled_df = pd.read_pickle("./dataset/TrainingSet/cv/folds_ss_performance.pkl")
    unpickled_df1 = pd.read_pickle("./dataset/TrainingSet/cv/average_ss_performance.pkl")
    unpickled_df2 = pd.read_pickle("./dataset/TrainingSet/cv/average_performance.pkl")

    # print(unpickled_df1)

    data = unpickled_df.to_csv(index=ss, index_label='Score Type')
    data1 = unpickled_df1.to_csv(index=ss, index_label='SS')
    data2 = unpickled_df2.to_csv(index_label='MEAN')
    with open('./dataframe.csv', 'w') as df, open('./dataframe1.csv', 'w') as df1, open('./dataframe2.csv', 'w') as df2:
        df.write(data)
        df1.write(data1)
        df2.write(data2)

    data = pd.read_csv("dataframe.csv") 
    data1 = pd.read_csv("dataframe1.csv")
    data2 = pd.read_csv("dataframe2.csv")
    # print(data)#,'\n\n',data2)
    # df_melt = pd.melt(data, id_vars="Score Type", var_name="Fold", value_name="%score")
    # print(df_melt)
    # df_plot = sns.catplot(x="Fold", y="%score", hue="Score Type", kind="bar", palette="ch:.25", data=df_melt)

    fig, ax = plt.subplots()
    # df1_melt = pd.melt(data1, id_vars="SS", var_name="Score Type", value_name="%score")
    # df1_plot = sns.catplot(ax=ax, x="SS", y="%score", hue="Score Type", kind="bar", palette="ch:.25", data=df1_melt)
    # change_width(ax, .25)
    # df1_plot.savefig('./dataset/TrainingSet/cv/average_ss_performance.png')

    fig1, ax1 = plt.subplots()
    df2_melt = pd.melt(data2, id_vars="MEAN", var_name="Score Type", value_name="%score")
    df2_plot = sns.catplot(ax=ax1, x="Score Type", y="%score", hue="Score Type", kind="bar", height=4, aspect=2, palette="ch:.25", data=df2_melt)
    change_width(ax, .15)
    plt.show(df2_plot)
    df2_plot.savefig('./dataset/TrainingSet/cv/average_performance.png')





    