#%% miscellaneous functions
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

def heatmap_sns(data, title:str= '', figsize:tuple=(15,8)):
    ''' This function allows larger displays for the heatmap compared to the imshow.

    As a drawback it depends of seaborn.
    '''
    fig, ax = plt.subplots(1,1, figsize=figsize)
    ax = sns.heatmap(data,ax =ax,
                cbar_kws={'label': 'Power [MW]'})
    ax.set_title ('')
    ax.set_xlabel ('Day of year')
    ax.set_ylabel ('Hour')
    ax.set_title(title)

def heatmap2d(arr: np.ndarray):
    plt.imshow(arr, cmap='hot', interpolation='gaussian')
    plt.xlabel('Day')
    plt.ylabel('Hour')
    #plt.savefig('maps.png') # <<<<<<<<<<<check operation
    ax = plt.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cax=cax)
    plt.ylabel('MW')#W/m${^2}$
    plt.show()
# %%
