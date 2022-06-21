import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator 

# References to data that we need to plot - reduced compared to JLat gating script
data_references = ['% Live cells','% GFP+ cells','MFI GFP+']

def plot_export (metadata, pairs, hue, figure_size=8, rotation=45, name="x", verbose = False, palette = "viridis"):
    # hue = [X-axis label, Legend labels]
    # Start by preparing the plot template etc.
    plt.close('all') # Close any open plots to avoid a mess
    plt.rcParams["figure.figsize"] = (figure_size*4,figure_size*1.25)
    # Prepare 4 subplot axes & appropriate subplot titles
    fig, ax = plt.subplots(1, 3, constrained_layout = True)  
    subplots = data_references.copy()
    # different axis spans
    subplot_titles = ['Live cells', 'GFP+ cells', 'Live cell MFI']
    # Sort by stimulus (and by timepoint if there are applicable)
    if ('Timepoint' in hue):  
        plt.rcParams["figure.figsize"] = (figure_size*12,figure_size)
        timepoints = ['0h', '1h', '2h', '4h', '8h', '16h', '24h', '30h', '48h']
        metadata ['Timepoint_tag'] = metadata ['Timepoint']
        metadata ['Timepoint_tag'] = metadata['Timepoint_tag'].replace(to_replace=timepoints, 
                                                                       value=range(len(timepoints)))
        metadata = metadata.sort_values(by=['Timepoint_tag','Stimulus'])
    else:
        metadata = metadata.sort_values(by=['iCasp'], ascending=False)
        metadata = metadata.sort_values(by=['Stimulus'])
    # Cycle through 4 gates and draw the subplots
    for i in range(len(ax)):
        # Index the current plot
        ax[i]= plt.subplot(1, 3, i+1)
        # Refine the plot -> appropriate legends and plots 
        if hue[1] != '':
            ax[i]= sns.barplot (x=hue[0], y=subplots[i],  hue=hue[1], data=metadata, 
                            dodge=True, palette = palette)
        else:
            ax[i]= sns.barplot (x=hue[0], y=subplots[i],  data=metadata, 
                            dodge=False, palette = palette)
        ax[i].tick_params(axis='x', rotation=rotation)
        if i < 2:
            ax[i].set_ylim([0,105])
        plt.grid(axis = 'y', color = 'gray', linestyle = '--', linewidth = 0.5)
        plt.title(subplot_titles[i])
        # Clean up the plot by removing excess labels 
        plt.xlabel("")
        if (i != 2):
            if hue[1] != '': ax[i].get_legend().remove()
        else:
            ax[i].legend(fontsize = 30, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
        if i <= 1:
            plt.ylabel("Fraction of events (%)")
        else:
            plt.ylabel("Mean Fluorescence Intensity")
        if pairs != []:
            annotator = Annotator(ax[i], pairs, data=metadata, x= hue[0], y=subplots[i], hue=hue[1])
            annotator.configure(test='t-test_welch', text_format='star', comparisons_correction="bonferroni")
            annotator.apply_and_annotate()
        # This has to come after data is plotted
        if i == 2:
            ax[i].set_ylim(200, None)
        # Scale the plot
        sns.set_context("talk", font_scale=2)

    # Change the working directory
    os.chdir(r"C:\Output folder directory")
    # Save the figure
    plt.savefig(name, bbox_inches='tight')
    # Optional visualisation
    if verbose:
        plt.show()
    return metadata


if __name__ == "__main__":
    # Pick working directory
    dir = r'C:\Input folder directory'
    os.chdir(dir)   # Change the working directory
    # Import the Excel summaries 
    metadata_kinetics = pd.read_excel("PMA Kinetics Experiment Summary.xlsx", index_col=None)
    metadata_qvdoph = pd.read_excel("Q-VD-OPh Experiment Summary.xlsx", index_col=None)
    metadata_potentiators =  pd.read_excel("PMA Enhancers Experiment Summary.xlsx", index_col=None)
    # Exclude Jurkat and 7AAD
    metadata_kinetics = metadata_kinetics [metadata_kinetics['Cell type'] != "Jurkat"]
    metadata_kinetics = metadata_kinetics [metadata_kinetics['7AAD'] != "-"]
    metadata_qvdoph = metadata_qvdoph [metadata_qvdoph['Cell type'] != "Jurkat"]
    metadata_qvdoph = metadata_qvdoph [metadata_qvdoph['7AAD'] != "-"]
    metadata_potentiators = metadata_potentiators [metadata_potentiators['Cell type'] != "Jurkat"]
    metadata_potentiators = metadata_potentiators [metadata_potentiators['7AAD'] != "-"]
    # Split and reduce the datasets
    metadata_kinetics_AS = metadata_kinetics[metadata_kinetics['Folder tag'] == 0].copy()
    metadata_kinetics_D2 = metadata_kinetics[metadata_kinetics['Folder tag'] == 2].copy()
    metadata_kinetics_mixed = metadata_kinetics[metadata_kinetics['Stimulus'] == 'PMA Iono '].copy()
    metadata_qvdoph_D2 = metadata_qvdoph[metadata_qvdoph['Folder tag'] == 2]
    metadata_potentiators_D3 = metadata_potentiators[metadata_potentiators['Folder tag'] == 2]
    # More tailored edits
    metadata_kinetics_AS = metadata_kinetics_AS[metadata_kinetics_AS['Timepoint'] != '1h']
    metadata_kinetics_AS = metadata_kinetics_AS[metadata_kinetics_AS['Timepoint'] != '30h']
    metadata_kinetics_AS = metadata_kinetics_AS[metadata_kinetics_AS['Timepoint'] != '48h']
    metadata_kinetics_D2 = metadata_kinetics_D2[metadata_kinetics_D2['Timepoint'] != '1h']
    metadata_kinetics_D2 = metadata_kinetics_D2[metadata_kinetics_D2['Timepoint'] != '30h']
    metadata_kinetics_mixed = metadata_kinetics_mixed[metadata_kinetics_mixed['Timepoint'] != '1h']
    metadata_kinetics_mixed = metadata_kinetics_mixed[metadata_kinetics_mixed['Timepoint'] != '30h']
    metadata_kinetics_mixed = metadata_kinetics_mixed[metadata_kinetics_mixed['Timepoint'] != '48h']
    metadata_kinetics_mixed = metadata_kinetics_mixed[metadata_kinetics_mixed['Folder tag'] != 3]
    potentiator_categories = ["CTL", "Iono", "JQ1", "RVX", "OXA", "PMA", "PMA Iono ", "PMA+JQ1", "PMA+RVX", "PMA+OXA"]
    metadata_potentiators_D3["Stimulus"] = pd.Categorical(metadata_potentiators_D3["Stimulus"], categories = potentiator_categories)
    metadata_potentiators_D3.sort_values(by = "Stimulus")
    # Pairs for the stats
    pairs_kinetics_AS = [[('2h','PMA Iono '),('2h','CTL')],
                         [('4h','PMA Iono '),('4h','CTL')],
                         [('8h','PMA Iono '),('8h','CTL')],
                         [('16h','PMA Iono '),('16h','CTL')],
                         [('24h','PMA Iono '),('24h','CTL')]
                         ]
    pairs_kinetics_D2 = [[('2h','PMA Iono '),('4h','PMA Iono ')],
                         [('8h','PMA Iono '),('16h','PMA Iono ')],
                         [('4h','PMA Iono '),('8h','PMA Iono ')],
                         [('16h','PMA Iono '),('24h','PMA Iono ')],
                         ]
    pairs_qvdoph_D2 = [[('CTL','no iCasp'),('CTL','iCasp')],
                       [('PMA','no iCasp'),('PMA','iCasp')],
                       [('PMA-Iono','no iCasp'),('PMA-Iono','iCasp')],
                       [('PMA-JQ1','no iCasp'),('PMA-JQ1','iCasp')],
                       [('PMA-OXA','no iCasp'),('PMA-OXA','iCasp')],
                       [('PMA','no iCasp'),('PMA-Iono','no iCasp')],
                       [('PMA-Iono','iCasp'),('PMA-JQ1','iCasp')],
                       [('PMA','no iCasp'),('PMA-JQ1','no iCasp')],
                       ]
    # Final plots
    metadata_kinetics_AS = plot_export(metadata_kinetics_AS, pairs_kinetics_AS, hue=['Timepoint','Stimulus'], name = "JLat_kinetics_D0", palette = "magma", rotation = 0)
    metadata_kinetics_D2 = plot_export(metadata_kinetics_D2, pairs_kinetics_D2, hue=['Timepoint','Stimulus'], name = "JLat_kinetics_D2", palette = "magma", rotation = 0)
    metadata_qvdoph_D2 = plot_export(metadata_qvdoph_D2, pairs_qvdoph_D2, hue=['Stimulus','iCasp'], name = "JLat_qvdoph_D2",palette = "rocket_r")
    metadata_potentiators_D3 = plot_export(metadata_potentiators_D3, [], hue=['Stimulus',''], name = "JLat_enhancers_D3", rotation = 90)
    metadata_kinetics_mixed = plot_export(metadata_kinetics_mixed, [], hue=['Timepoint','Folder tag'], name = "JLat_kinetics_mixed", palette = "flare", rotation = 0)
