import numpy as np
import pandas as pd
import operator
import math
import matplotlib.pyplot as plt 
from sklearn.decomposition import PCA


'''
Author: YC3

Run with:
python PCA_scatter.py
'''

def plot_pca(data, targets, colors):
     
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(data.T)
    pcs_df = pd.DataFrame(data=pcs, columns=['PC1', 'PC2'])

    finalDf = pd.concat([pcs_df, pd.DataFrame({'target':data.columns.values})], axis = 1)

    fig = plt.figure(figsize = (6,5))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('PC1', fontsize = 15)
    ax.set_ylabel('PC2', fontsize = 15)

    for target, color in zip(targets,colors):
        indicesToKeep = finalDf['target'] == target
        ax.scatter(finalDf.loc[indicesToKeep, 'PC1'], finalDf.loc[indicesToKeep, 'PC2'], c = color, s = 50)

    ax.legend(targets, loc='center left', bbox_to_anchor=(1, 0.5))
    ax.grid()
  


if __name__ == '__main__':
 
    if_plot = True
    # read in data
    exp = pd.read_csv("GSE139242_normalised_abundance_fpkm_all_samples.csv")

    # preprocessing
    exp2 = exp.groupby('Geneid').mean()  
    exp3 = np.log(exp2 + 1)
    exp3.to_csv('exp_processed.csv')

    # PCA plot
    colors = ['red']*4 + ['blue']*4 + ['lightblue'] + ['green']*5 + ['orange']*3 + ['gold']
    plot = plot_pca(data = exp3, targets = exp3.columns.values, colors = colors)
    plt.show()
    print('Result: it is reasonalble to combine "bloodinfant" and "bloodadult"')
