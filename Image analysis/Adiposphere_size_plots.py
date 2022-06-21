import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator 


def input_xlsx_data (input_dir):
    os.chdir(input_dir)   # Change the working directory
    file_names = os.listdir() # Check the files in the directory
    metadata = pd.DataFrame()
    for filename in file_names:
        excel_tabs = pd.ExcelFile(filename).sheet_names 
        for tab in excel_tabs:
           data = pd.read_excel(filename, sheet_name=tab, index_col=1, header = None)
           # Remove all nans
           data = data.iloc[:, 1:(data.shape[1]-1)] # data.shape[1] gives number of columns
           while data.iloc[0].isnull().values.any():
               data = data.iloc[1:(data.shape[0]),:] # data.shape[1] gives number of rows
           while data.iloc[-1:].isnull().values.any():
               data = data.iloc[0:(data.shape[0]-1),:]
           # Remove false colnames
           data = data.iloc[1:(data.shape[0]),:]
           # Set colnames
           data.rename(columns = {2:'Area', 3:'Perim.', 4:'Major', 5:'Minor',
                                  6:'Angle', 7:'Circ.', 8:'Feret', 9:'FeretX',
                                  10:'FeretY', 11:'FeretAngle', 12:'MinFeret',
                                  13:'AR', 14:'Round', 15:'Solidity', 16:'Volume',
                                  17:'oblate ellipse'}, inplace = True)
           for day in ['J0','J4','J8','J14','J17','J24','J30','J44',]:
               if day in filename:
                   data ['Day'] = day 
           for hSVF in ['1', '2', '3', '4']:
               if hSVF in tab:
                   data ['Donor'] = "hSVF" + hSVF
           media = 'EGM'
           for lipid in ['HG', 'LG']:
               if lipid in tab:
                   media = lipid
           data ['Media type'] = media
           metadata = metadata.append(data)
    return metadata 

def plot_oblatElipse_perP (metadata, hSVF, exceptions, output_dir, plot, verbose = True):
    metadata_toplot = metadata[metadata['Donor'] == hSVF]
    for exception in exceptions : 
        metadata_toplot = metadata_toplot[metadata_toplot['Day'] != exception]    
    # Create a figure
    plt.rcParams["figure.figsize"] = (12, 12)
    sns.set(font_scale=3) # Originally 1.5
    plotting_parameters = {'data': metadata_toplot,
                            'x': 'Day',
                            'y': 'oblate ellipse',
                            'hue': 'Media type',
                            'palette' : "YlOrBr",
                            'dodge': True}
    # Plot a barplot or a boxplot 
    if plot == 0:
        ax = sns.boxplot (**plotting_parameters) 
    elif plot == 1:
        ax = sns.barplot (**plotting_parameters)
    # Tune some visual parameters
    ax.set(xlabel=None)
    ax.legend(fontsize = 30, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
    ax.set_title(hSVF)
    plt.grid(axis = 'y', color = 'gray', linestyle = '--', linewidth = 0.5)
    # Add the statistics 
    if 'J17' in metadata_toplot['Day'].values:
        stat_pairs = [[('J0', 'EGM'),('J4', 'EGM')],
                  [('J8', 'LG'),('J8', 'HG')],
                  [('J14', 'LG'),('J14', 'HG')],
                  [('J17', 'LG'),('J17', 'HG')],
                  [('J24', 'LG'),('J24', 'HG')],
                  [('J30', 'LG'),('J30', 'HG')],
                  [('J44', 'LG'),('J44', 'HG')],
                  [('J4', 'EGM'),('J14', 'LG')]]
    else:
        stat_pairs = [[('J0', 'EGM'),('J4', 'EGM')],
                  [('J8', 'LG'),('J8', 'HG')],
                  [('J14', 'LG'),('J14', 'HG')],
                  [('J24', 'LG'),('J24', 'HG')],
                  [('J30', 'LG'),('J30', 'HG')],
                  [('J44', 'LG'),('J44', 'HG')],
                  [('J4', 'EGM'),('J14', 'LG')]]
    annotator = Annotator(ax, stat_pairs, data=metadata_toplot, x='Day', y='oblate ellipse', hue = 'Media type')
    annotator.configure(test='Mann-Whitney', text_format='star', comparisons_correction="bonferroni")
    annotator.apply_and_annotate()
    # Save the plots
    os.chdir(output_dir)   # Change the working directory
    if plot == 0:
        plt.savefig('Boxplot '+ hSVF +'.png', bbox_inches='tight')
    elif plot == 1:
        plt.savefig('Barplot '+ hSVF +'.png', bbox_inches='tight')
    # Optional visualisation 
    if verbose:
        plt.show()
    plt.close('all') # Close any open plots to avoid a mess
    return metadata_toplot

if __name__ == "__main__":
    input_dir = r'C:\Input folder directory'
    output_dir = r'C:\Output folder directory'
    metadata = input_xlsx_data (input_dir)
    # Convert um^3 to mm^3
    metadata ['oblate ellipse'] = (metadata ['oblate ellipse']*0.000000001).astype('float64')
    # Additional parameters
    hSVF = ['hSVF1','hSVF2','hSVF3','hSVF4']
    exceptions = [['J17'],['17'],[],[]]
    for k in range (4):
        for i in [0,1]:
            metadata_toplot = plot_oblatElipse_perP (metadata, hSVF[k], exceptions[k], 
                                   output_dir, plot = i)