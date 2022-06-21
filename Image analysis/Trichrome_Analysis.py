import pandas as pd
import numpy as np
from PIL import Image, ImageOps
import matplotlib.pyplot as plt
import seaborn as sns

# Read return metadata from a file name (e.g. media type or microscope magnification)
def photo_metadata_read (file_name):
    #  Normal or obesogenic?
    if ('norm' in file_name) or ('Norm' in file_name):
        lipid = 'Normogenic'
    elif 'obes' in file_name:
        lipid = 'Obesogenic'
    # Ascorbate?
    if 'ASCp' in file_name:
        ASC = '+'
    elif 'ASC-' in file_name:
        ASC = '-' 
    # Magnificiation?
    if 'x10' in file_name:
        magnification = 'x10'
    elif 'x5' in file_name:
        magnification = 'x5'
    # Replicate?
    replicate = file_name [-5]
    while replicate == ' ':
        print ('Error!')
        import time
        time.sleep(1) # Sleep for 3 seconds 
    # Donor?
    donor = file_name[5]
    # Collate and return the data
    df = pd.DataFrame ([[file_name, lipid, ASC, magnification, 
                         replicate, donor]],
                       columns=('File Name','Media type', 'Ascorbic acid',
                                'Magnificiation', 'Replicate', 
                                'Donor'))
    return df

#♠ Collate data from multiple files across several input folders
def photo_metadata_collate (dir, folder_input):
    #♠ Gather and order metadata of files across all folders
    metadata = pd.DataFrame()
    for folder in folder_input:
        os.chdir(dir + folder)   # Change the working directory
        file_names = os.listdir()    # Check the files in the directory
        for file in file_names:
            metadata = metadata.append(photo_metadata_read(file))
    return metadata 

# Analyse collagen content in an image aka blue channel
def photo_binary_analysis (image_name, verbose = False):
    plt.rcParams["figure.figsize"] = (20, 20)
    # Load image and prep a dataset
    input_image = Image.open(image_name)
    width = input_image.width
    height = input_image.height
    monocolor = np.empty((width, height))
    # Invert the input
    inverted = ImageOps.invert(input_image)
    input_pixels = inverted.load()
    # Extract one channel info
    non_zero_counter = 0
    for x in range(width):
        for y in range(height):
            pixel = input_pixels[x, y]
            monocolor[x, y] = pixel[2] / 3
            #•Threshold the values
            if monocolor[x, y]  < 5:
                monocolor[x, y]  = 0
            else:
                non_zero_counter = non_zero_counter + 1
    # Optional visualisation                 
    if verbose:
        plt.imshow(monocolor, cmap='gray')
        plt.show()
    return (100-non_zero_counter/width/height*100)

# Analyse several photos in a batch using functions above
def bulk_photo_analysis (dir, folder_input, verbose = True):
    metadata = photo_metadata_collate (dir, folder_input)
    means = [] 
    for index in range(len (metadata)):
        os.chdir(dir + folder_input[int(metadata.iloc[index]['Donor'])-1])
        means.append (photo_binary_analysis(metadata.iloc[index]['File Name'], verbose = verbose))
    metadata ['Percentage'] = means
    return metadata

# Export results as an Excel file for quicker access later
def export_microscopy_data (metadata, dir, folder_output):
    os.chdir(dir + folder_output)   # Change the working directory
    metadata.loc[:, metadata.columns.isin(['File Name','Media type','Ascorbic acid','Magnificiation','Replicate','Donor','Percentage'])].to_excel("[Microscopy] Data Summary.xlsx")
    return

# Import previously exported results
def import_microscopy_data (dir, folder_output,filename = "[Microscopy] Data Summary.xlsx"):
    os.chdir(dir + folder_output)   # Change the working directory
    return pd.read_excel(filename, index_col=None)

# Plotting the results with pooled ASC+ and ASC-
def barplot_it_condensed_ASC (metadata, plot, dir, folder_output, tag = "", verbose = True):
    # Annotations: https://pythonlang.dev/repo/trevismd-statannotations/
    from statannotations.Annotator import Annotator 
    metadata ['hSVF'] = 'hSVF#'+ metadata['Donor'].apply(str) 
    # Sort the datasets 
    metadata = metadata.sort_values(by=['hSVF','Media type'])
    # Create a figure
    plt.rcParams["figure.figsize"] = (12, 12)
    sns.set(font_scale=3)
    sns.set_style("ticks")
    plotting_parameters = {'data': metadata,
                            'x': 'hSVF',
                            'y': 'Percentage',
                            'hue': 'Media type',
                            'palette': "winter_r",
                            'dodge': True}
    # Plot a barplot, a boxplot or a violin plot
    if plot == 0:
        ax = sns.boxplot (**plotting_parameters) 
    elif plot == 1:
        ax = sns.barplot (**plotting_parameters)
    elif plot == 2:
        ax = sns.violinplot (**plotting_parameters)        
    # Add stats annotations
    pairs=[[('hSVF#1', 'Obesogenic'),('hSVF#1', 'Normogenic')],
           [('hSVF#2', 'Obesogenic'),('hSVF#2', 'Normogenic')],
           [('hSVF#3', 'Obesogenic'),('hSVF#3', 'Normogenic')],
           [('hSVF#4', 'Obesogenic'),('hSVF#4', 'Normogenic')]]
    annotator = Annotator(ax, pairs, data=metadata, x='hSVF', y='Percentage', hue = 'Media type')
    annotator.configure(test='Mann-Whitney', text_format='star', comparisons_correction="bonferroni")
    annotator.apply_and_annotate()
    # Tune some visual parameters
    ax.set(xlabel=None)
    ax.legend(fontsize = 30, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
    plt.ylabel("Surface area labelled as collagen (%)", fontsize = 35)
    ax.set(ylim=(0, 100))
    plt.grid(axis = 'y', color = 'gray', linestyle = '--', linewidth = 0.5)
    # Save the figure
    os.chdir(dir + folder_output)   # Change the working directory
    if plot == 0:
        plt.savefig(tag + 'Condensed boxplot (no ASC).png', bbox_inches='tight')
    elif plot == 1:
        plt.savefig(tag + 'Condensed barplot (no ASC).png', bbox_inches='tight')
    elif plot == 2:
        plt.savefig(tag + 'Condensed violin plot (no ASC).png', bbox_inches='tight')
    if verbose:
        plt.show()
    plt.close('all') # Close any open plots to avoid a mess
    return 

# Plotting the results with pooled HG and LG
def barplot_it_condensed_lipids (metadata, plot, dir, folder_output, tag = "", verbose = True):
    # Annotations: https://pythonlang.dev/repo/trevismd-statannotations/
    from statannotations.Annotator import Annotator 
    metadata ['hSVF'] = 'hSVF#'+ metadata['Donor'].apply(str) 
    # Sort the datasets 
    metadata = metadata.sort_values(by=['hSVF','Media type'])
    # Create a figure
    plt.rcParams["figure.figsize"] = (12, 12)
    sns.set(font_scale=3)
    sns.set_style("ticks")
    plotting_parameters = {'data': metadata,
                            'x': 'hSVF',
                            'y': 'Percentage',
                            'hue': 'Ascorbic acid',
                            'palette': "winter_r",
                            'dodge': True}
    # Plot a barplot, a boxplot or a violin plot
    # https://seaborn.pydata.org/generated/seaborn.swarmplot.html
    if plot == 0:
        ax = sns.boxplot (**plotting_parameters) 
    elif plot == 1:
        ax = sns.barplot (**plotting_parameters)
    elif plot == 2:
        ax = sns.violinplot (**plotting_parameters)        
    # Add stats annotations
    pairs=[[('hSVF#1', '+'),('hSVF#1', '-')],
           [('hSVF#2', '+'),('hSVF#2', '-')],
           [('hSVF#3', '+'),('hSVF#3', '-')],
           [('hSVF#4', '+'),('hSVF#4', '-')]]
    annotator = Annotator(ax, pairs, data=metadata, x='hSVF', y='Percentage', hue = 'Ascorbic acid')
    annotator.configure(test='Mann-Whitney', text_format='star', comparisons_correction="bonferroni")
    annotator.apply_and_annotate()
    # Tune some visual parameters
    ax.set(xlabel=None)
    ax.legend(fontsize = 30, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
    plt.ylabel("Surface area labelled as collagen (%)", fontsize = 35)
    ax.set(ylim=(0, 100))
    plt.grid(axis = 'y', color = 'gray', linestyle = '--', linewidth = 0.5)
    # Save the figure
    os.chdir(dir + folder_output)   # Change the working directory
    if plot == 0:
        plt.savefig(tag + 'Condensed boxplot (no lipids).png', bbox_inches='tight')
    elif plot == 1:
        plt.savefig(tag + 'Condensed barplot (no lipids).png', bbox_inches='tight')
    elif plot == 2:
        plt.savefig(tag + 'Condensed violin plot (no lipids).png', bbox_inches='tight')
    if verbose:
        plt.show()
    plt.close('all') # Close any open plots to avoid a mess
    return 

# Plotting the results with separated conditions
def barplot_it_detailed (metadata, plot, dir, folder_output, tag = "", verbose = True):
    # Annotations: https://pythonlang.dev/repo/trevismd-statannotations/
    from statannotations.Annotator import Annotator 
    # Prepare additional functions to do statistics and plots
    metadata ['hSVF'] = "hSVF#"+ metadata['Donor'].apply(str) 
    metadata['Growth media'] = metadata['Media type'] + ' (ASC' + metadata['Ascorbic acid'] + ')'
    metadata = metadata.sort_values(by=['hSVF','Media type'])
    # Create a figure
    plt.rcParams["figure.figsize"] = (12, 12)
    sns.set(font_scale=3)
    sns.set_style("ticks")
    plotting_parameters = {'data': metadata,
                            'x': 'hSVF',
                            'y': 'Percentage',
                            'hue': 'Growth media',
                            'palette': "winter_r",
                            'dodge': True}
    # Plot a barplot or a boxplot
    if plot == 0:
        ax = sns.boxplot (**plotting_parameters) 
    elif plot == 1:
        ax = sns.barplot (**plotting_parameters)
    elif plot == 2:
        ax = sns.violinplot (**plotting_parameters)   
    # Add stats annotations
    pairs=[[('hSVF#1', 'Obesogenic (ASC-)'),('hSVF#1', 'Obesogenic (ASC+)')],
           [('hSVF#1', 'Normogenic (ASC-)'),('hSVF#1', 'Normogenic (ASC+)')],
           [('hSVF#2', 'Obesogenic (ASC-)'),('hSVF#2', 'Obesogenic (ASC+)')],
           [('hSVF#2', 'Normogenic (ASC-)'),('hSVF#2', 'Normogenic (ASC+)')],
           [('hSVF#3', 'Obesogenic (ASC-)'),('hSVF#3', 'Obesogenic (ASC+)')],
           [('hSVF#3', 'Normogenic (ASC-)'),('hSVF#3', 'Normogenic (ASC+)')],
           [('hSVF#4', 'Obesogenic (ASC-)'),('hSVF#4', 'Obesogenic (ASC+)')],
           [('hSVF#4', 'Normogenic (ASC-)'),('hSVF#4', 'Normogenic (ASC+)')],]
    annotator = Annotator(ax, pairs, data=metadata, x='hSVF', y='Percentage', hue = 'Growth media')
    annotator.configure(test='Mann-Whitney', text_format='star', comparisons_correction="bonferroni")
    annotator.apply_and_annotate()
    # Tune some visualparameters
    ax.set(xlabel=None)
    ax.legend(fontsize = 30, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
    ax.set(ylim=(0, 100))
    plt.ylabel("Surface area labelled as collagen (%)", fontsize = 35)
    plt.grid(axis = 'y', color = 'gray', linestyle = '--', linewidth = 0.5)
    # Save the figure
    os.chdir(dir + folder_output)   # Change the working directory
    if plot == 0:
        plt.savefig(tag + 'Detailed boxplot.png', bbox_inches='tight')
    elif plot == 1:
        plt.savefig(tag + 'Detailed barplot.png', bbox_inches='tight')
    elif plot == 2:
        plt.savefig(tag + 'Detailed violin plot.png', bbox_inches='tight')
    if verbose:
        plt.show()
    plt.close('all') # Close any open plots to avoid a mess
    return 

    
if __name__ == "__main__":
    import os
    dir = r'C:\Input folder directory'
    folder_input = ['\Input folder directory #1',
               '\Input folder directory #2',
               '\Input folder directory #3',
               '\Input folder directory #4']
    folder_output = r'\Output folder directory'
    # Update (uncomment) the section below to use pre-exported data
    tag = "[Microscopy] "
    metadata_microscopy = bulk_photo_analysis (dir, folder_input, verbose = False)
    export_microscopy_data(metadata_microscopy, dir, folder_output)
    # metadata_microscopy = import_microscopy_data(dir, folder_output)
    for plot in range(0,3):
        barplot_it_condensed_ASC(metadata_microscopy, plot, dir, folder_output, tag=tag)
        barplot_it_condensed_lipids(metadata_microscopy, plot, dir, folder_output, tag=tag)
        barplot_it_detailed(metadata_microscopy, plot, dir, folder_output, tag=tag)
    