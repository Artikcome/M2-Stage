import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator 

# References to data that we need to plot - reduced compared to JLat gating script
data_references = ['Total GFP+ cells','% Live cells','% GFP+ cells','Total Live Cells',]

def plot_characteristics_detailed (metadata_coculture, figure_size=4, verbose=False, palette="cubehelix_r"):
    metadata = metadata_coculture.copy ()
    metadata.loc[metadata['Donor'].ne('-'), 'Donor'] = "[co-culture]"
    metadata.loc[metadata['Donor'].eq('-'), 'Donor'] = "[control]"
    # Add legend, remove out-of-scope datapoints and sort
    metadata['Legend'] = metadata['Stimulus'] +' in '+ metadata['Media'] + " " + metadata['Donor']
    temp = metadata [metadata['Folder tag'] == 0]
    metadata = metadata [(metadata['Folder tag'] != 4) & (metadata['Media'] != 'RPMI')]
    metadata = metadata.sort_values(by=['Stimulus','Donor', 'Media']) 
    metadata = temp.append (metadata)
    # Make the plot !
    for i in range(len(data_references)):
        hue = ['Folder tag', data_references[i], 'Legend']
        plt.rcParams["figure.figsize"] = (figure_size*2.5,figure_size)
        ax = sns.barplot (x=hue[0], y=hue[1],  hue=hue[2], data=metadata, 
                        dodge=True, palette = palette)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
        plt.xlabel("Co-culture day")
        # Optional log scale
        if i in [0,3]:
            plt.yscale('symlog')
        # Change the working directory
        os.chdir(r"C:\Output folder directory")
        # Save the figure
        plt.savefig("Co-culture_characteristicsDetailed_" + hue[1], bbox_inches='tight')
        # Optional visualisation
        if verbose:
            plt.show()
        ax.clear()
    return metadata

def plot_characteristics_condensed (metadata_coculture, figure_size=5, verbose=False, palette="cubehelix_r"):
    metadata = metadata_coculture.copy ()
    metadata.loc[metadata['Donor'].ne('-'), 'Donor'] = "Co-culture"
    metadata.loc[metadata['Donor'].eq('-'), 'Donor'] = "Monoculture"
    metadata.loc[metadata['Media'].eq('HG'), 'Media'] = "Obesogenic"
    metadata.loc[metadata['Media'].eq('LG'), 'Media'] = "Normogenic"
    # Add legend, remove out-of-scope datapoints and sort
    metadata = metadata [(metadata['Folder tag'] == 3)&(metadata['Media'] != 'RPMI')]
    metadata = metadata.sort_values(by=['Media']) 
    metadata = metadata.sort_values(by=['Donor'], ascending = False) 
    metadata = metadata.sort_values(by=['Stimulus']) 
    # Make the plot !
    for i in range(len(data_references)):
        if i == 3:
            metadata['Legend'] = metadata['Stimulus'] + " " + metadata['Media']
            hue = ['Donor', data_references[i], 'Legend']
        else:
            metadata['Legend'] = metadata['Stimulus'] + " " + metadata['Donor']
            hue = ['Media', data_references[i], 'Legend']
        plt.rcParams["figure.figsize"] = (figure_size,figure_size)
        ax = sns.barplot (x=hue[0], y=hue[1],  hue=hue[2], data=metadata, 
                        dodge=True, palette = palette)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
        # Optional log scale
        if i in [0,3]:
            plt.yscale('log')
        # Statistics! 
        pairs = []
        if i == 2:
            pairs = [[('Obesogenic','PMA-JQ1 Monoculture'),('Obesogenic','PMA-JQ1 Co-culture')],
                     [('Normogenic','PMA-JQ1 Monoculture'),('Normogenic','PMA-JQ1 Co-culture')],
                     [('Obesogenic','PMA-JQ1 Monoculture'),('Normogenic','PMA-JQ1 Monoculture')],
                     [('Obesogenic','PMA-JQ1 Co-culture'),('Normogenic','PMA-JQ1 Co-culture')]]
        elif i == 3:
            pairs = [[('Monoculture','PMA-JQ1 Normogenic'),('Co-culture','PMA-JQ1 Normogenic')],
                     [('Monoculture','CTL Obesogenic'),('Monoculture','CTL Normogenic')],
                     [('Co-culture','CTL Obesogenic'),('Co-culture','CTL Normogenic')],
                     [('Monoculture','PMA-JQ1 Obesogenic'),('Monoculture','PMA-JQ1 Normogenic')],
                     [('Co-culture','PMA-JQ1 Obesogenic'),('Co-culture','PMA-JQ1 Normogenic')]]
        if pairs != []:
            annotator = Annotator(ax, pairs, data=metadata, x= hue[0], y=hue[1],  hue=hue[2])
            annotator.configure(test='t-test_welch', text_format='star', comparisons_correction="bonferroni")
            annotator.apply_and_annotate()
        # Scale the plot
        sns.set_context("talk")
        # Change the working directory
        os.chdir(r"C:\Output folder directory")
        # Save the figure
        plt.savefig("Co-culture_characteristicsCondensed_" + hue[1], bbox_inches='tight')
        # Optional visualisation
        if verbose:
            plt.show()
        ax.clear()
    return metadata

# Alternative visualisation as percent stacked barplot 
# https://www.python-graph-gallery.com/stacked-and-percent-stacked-barplo
# https://stackoverflow.com/questions/35692781/python-plotting-percentage-in-seaborn-bar-plot
def plot_migration (metadata_coculture, figure_size=5, verbose=False, palette="cubehelix_r"):
    metadata = metadata_coculture.copy ()
    # Add legend, remove out-of-scope datapoints and sort
    metadata = metadata [((metadata['Folder tag'] == 3)|(metadata['Folder tag'] == 4)) & (metadata['Media'] != 'RPMI')]
    metadata = metadata.sort_values(by=['Media'], ascending = False) 
    metadata = metadata.sort_values(by=['Stimulus']) 
    # Fix some metadata
    metadata.loc[metadata['Media'].eq('HG'), 'Media'] = "Obesogenic"
    metadata.loc[metadata['Media'].eq('LG'), 'Media'] = "Normogenic"
    metadata.loc[metadata['Folder tag'].eq(3), 'Folder tag'] = 'Insert'
    metadata.loc[metadata['Folder tag'].eq(4), 'Folder tag'] = 'Well'
    # Add legend
    metadata['Legend'] = metadata['Stimulus'] +' in '+ metadata['Folder tag'] 
    # Make the plot !
    for i in range(len(data_references)):
        # Optional log scale
        if i in [0,3]:
            plt.semilogy()
            # Scale the insert data: Well volume is 20 times the insert sample! 
            metadata.loc[metadata['Folder tag'].eq('Well'), data_references[i]] = metadata[data_references[i]]/20
        hue = ['Media', data_references[i], 'Legend']
        plt.rcParams["figure.figsize"] = (figure_size,figure_size)
        ax = sns.barplot (x=hue[0], y=hue[1],  hue=hue[2], data=metadata, 
                        dodge=True, palette = palette)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
       
        # Statistics! 
        pairs = [[('Normogenic','CTL in Insert'),('Normogenic','CTL in Well')],
                 [('Normogenic','PMA-JQ1 in Insert'),('Normogenic','PMA-JQ1 in Well')],
                 [('Obesogenic','CTL in Insert'),('Obesogenic','CTL in Well')],
                 [('Obesogenic','PMA-JQ1 in Insert'),('Obesogenic','PMA-JQ1 in Well')]]
        annotator = Annotator(ax, pairs, data=metadata, x= hue[0], y=hue[1],  hue=hue[2])
        annotator.configure(test='t-test_welch', text_format='star', comparisons_correction="bonferroni")
        annotator.apply_and_annotate()
        # Scale the plot & other corrections
        #sns.set_context("talk")
        # Change the working directory
        os.chdir(r"C:\Output folder directory")
        # Save the figure
        plt.savefig("Co-culture_migration_" + hue[1], bbox_inches='tight')
        # Optional visualisation
        if verbose:
            plt.show()
        ax.clear()
    return metadata

def plot_patientVariability (metadata_coculture, figure_size=5, verbose=False, palette="cubehelix_r"):
    metadata = metadata_coculture.copy ()
    # Add legend, remove out-of-scope datapoints and sort
    metadata = metadata [(metadata['Folder tag'] != 0) & (metadata['Donor'] != '-')]
    metadata['Legend'] = metadata['Stimulus'] +' with hSVF#'+ metadata['Donor'].str[1:]
    metadata = metadata.sort_values(by=['Stimulus','Donor', 'Media']) 
    # Make the plot !
    for i in range(len(data_references)):
        hue = ['Folder tag', data_references[i], 'Legend']
        plt.rcParams["figure.figsize"] = (figure_size*2.5,figure_size)
        ax = sns.barplot (x=hue[0], y=hue[1],  hue=hue[2], data=metadata, 
                        dodge=True, palette = palette)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
        plt.xlabel("Co-culture day")
        # Optional log scale
        if i in [0,3]:
            plt.yscale('symlog')
        # Statistics! 
        pairs = [[(1,'CTL with hSVF#1'),(1,'CTL with hSVF#2')],
                 [(1,'CTL with hSVF#3'),(1,'CTL with hSVF#4')],
                 [(1,'CTL with hSVF#2'),(1,'CTL with hSVF#3')],
                 [(1,'PMA-JQ1 with hSVF#1'),(1,'PMA-JQ1 with hSVF#2')],
                 [(1,'PMA-JQ1 with hSVF#2'),(1,'PMA-JQ1 with hSVF#4')],
                 [(2,'CTL with hSVF#1'),(2,'CTL with hSVF#2')],
                 [(2,'CTL with hSVF#3'),(2,'CTL with hSVF#4')],
                 [(2,'CTL with hSVF#2'),(2,'CTL with hSVF#3')],
                 [(2,'PMA-JQ1 with hSVF#1'),(2,'PMA-JQ1 with hSVF#2')],
                 [(2,'PMA-JQ1 with hSVF#2'),(2,'PMA-JQ1 with hSVF#4')],
                 [(3,'CTL with hSVF#1'),(3,'CTL with hSVF#2')],
                 [(3,'CTL with hSVF#3'),(3,'CTL with hSVF#4')],
                 [(3,'CTL with hSVF#2'),(3,'CTL with hSVF#3')],
                 [(3,'PMA-JQ1 with hSVF#1'),(3,'PMA-JQ1 with hSVF#2')],
                 [(3,'PMA-JQ1 with hSVF#2'),(3,'PMA-JQ1 with hSVF#4')],
                 [(4,'CTL with hSVF#1'),(4,'CTL with hSVF#2')],
                 [(4,'CTL with hSVF#3'),(4,'CTL with hSVF#4')],
                 [(4,'CTL with hSVF#2'),(4,'CTL with hSVF#3')],
                 [(4,'PMA-JQ1 with hSVF#1'),(4,'PMA-JQ1 with hSVF#2')],
                 [(4,'PMA-JQ1 with hSVF#2'),(4,'PMA-JQ1 with hSVF#4')]
                 ]
        annotator = Annotator(ax, pairs, data=metadata, x= hue[0], y=hue[1],  hue=hue[2])
        annotator.configure(test='t-test_welch', text_format='star', comparisons_correction="bonferroni")
        annotator.apply_and_annotate()
        # Scale the plot
        #sns.set_context("talk")
        # Change the working directory
        os.chdir(r"C:\Output folder directory")
        # Save the figure
        plt.savefig("Co-culture_patientVariability_" + hue[1], bbox_inches='tight')
        # Optional visualisation
        if verbose:
            plt.show()
        ax.clear()
    return metadata

def plot_LTculture (metadata_LTculture, metadata_coculture, figure_size=5, verbose=False, palette="dark:green_r"):
    # Merge two datasets
    metadata = metadata_coculture.copy ()
    metadata = metadata [metadata['Media'] != "RPMI"]
    metadata = metadata [metadata['Folder tag'] != 4] # Keep insert data only
    metadata = metadata.append(metadata_LTculture)
    metadata = metadata [metadata['Folder tag'] != 0] # remove day 0
    #Sort
    metadata = metadata.sort_values(by=['Donor']) 
    # Merge different hSVF
    metadata.loc[metadata['Donor'].ne('-'), 'Donor'] = "Co-culture"
    metadata.loc[metadata['Donor'].eq('-'), 'Donor'] = "Monoculture"
    # Keep LG only
    metadata = metadata [metadata['Media'] != "HG"]
    # Add legend, remove out-of-scope datapoints and sort
    metadata['Legend'] = metadata['Stimulus'] + " " + metadata['Donor']
    for i in range(len(data_references)):
        hue = ['Folder tag', data_references[i], 'Legend']
        plt.rcParams["figure.figsize"] = (figure_size,figure_size)

        # remove CTL if looking at GFP
        if i in [0, 2]:
            metadata_toplot = metadata [metadata['Stimulus'] != "CTL"]
        else:
            metadata_toplot = metadata
        ax = sns.barplot (x=hue[0], y=hue[1],  hue=hue[2], data=metadata_toplot, 
                        dodge=True, palette = palette)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
        plt.xlabel("Co-culture day")
        # Optional log scale
        if i in [0,3]:
            plt.semilogy()
        # Statistics! 
        pairs = [[(1,'PMA-JQ1 Co-culture'),(1,'PMA-JQ1 Monoculture')],
                 [(2,'PMA-JQ1 Co-culture'),(2,'PMA-JQ1 Monoculture')],
                 [(3,'PMA-JQ1 Co-culture'),(3,'PMA-JQ1 Monoculture')]]
        annotator = Annotator(ax, pairs, data=metadata, x= hue[0], y=hue[1],  hue=hue[2])
        annotator.configure(test='t-test_welch', text_format='star', comparisons_correction="bonferroni")
        annotator.apply_and_annotate()
        # Change the working directory
        os.chdir(r"C:\Output folder directory")
        # Save the figure
        plt.savefig("LTculture_detailed_" + hue[1], bbox_inches='tight')
        # Optional visualisation
        if verbose:
            plt.show()
        ax.clear()
    return metadata

def plot_restimulation (metadata_restimulation, figure_size=5, verbose=False, palette="dark:green"):
    metadata = metadata_restimulation.copy ()
    # Correct some metadata
    metadata.loc[metadata['Stimulus'].eq('PMA'), 'Stimulus'] = 'PMA-JQ1'
    metadata.loc[metadata['Stimulus'].eq('PMA-CTL'), 'Stimulus'] = 'PMA-JQ1'
    metadata.loc[metadata['Stimulus'].eq('CTL-PMA'), 'Stimulus'] = 'CTL restimulated'
    metadata.loc[metadata['Stimulus'].eq('PMA-PMA'), 'Stimulus'] = 'PMA-JQ1 restimulated'
    metadata.loc[metadata['Folder tag'].eq(0), 'Folder tag'] = 16
    metadata.loc[metadata['Folder tag'].eq(1), 'Folder tag'] = 17
    for i in range(len(data_references)):
        hue = ['Folder tag', data_references[i], 'Stimulus']
        plt.rcParams["figure.figsize"] = (figure_size,figure_size)
        ax = sns.barplot (x=hue[0], y=hue[1],  hue=hue[2], data=metadata, 
                        dodge=True, palette = palette)
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)
        plt.xlabel("Culture day")
        # Optional log scale
        if i in [0,3]:
            plt.yscale('symlog')
            plt.semilogy()
        # Statistics! 
        pairs = [[(16,'PMA-JQ1'),(16,'CTL')],
                 [(17,'PMA-JQ1'),(17,'CTL')],
                 [(17,'PMA-JQ1'),(17,'CTL restimulated')],
                 [(17,'PMA-JQ1 restimulated'),(17,'CTL restimulated')]
                 ]
        annotator = Annotator(ax, pairs, data=metadata, x= hue[0], y=hue[1],  hue=hue[2])
        annotator.configure(test='t-test_welch', text_format='star', comparisons_correction="bonferroni")
        annotator.apply_and_annotate()
        # Change the working directory
        os.chdir(r"C:\Output folder directory")
        # Save the figure
        plt.savefig("Reinduction_" + hue[1], bbox_inches='tight')
        # Optional visualisation
        if verbose:
            plt.show()
        ax.clear()
    return metadata

if __name__ == "__main__":
    # Pick working directory
    dir = r'C:\Input folder directory'
    os.chdir(dir)   # Change the working directory
    # Import and pool Excel summaries 
    metadata_coculture = pd.read_excel("Co-culture FACS Experiment Summary.xlsx", index_col=None)
    metadata_coculture['Experiment tag'] = '0'
    metadata_coculture_CTL = pd.read_excel("Co-culture FACS controls Experiment Summary.xlsx", index_col=None)
    metadata_coculture_CTL['Experiment tag'] = '1'
    metadata_coculture = metadata_coculture.append(metadata_coculture_CTL, ignore_index=True)
    # Exclude Jurkat and 7AAD-
    metadata_coculture = metadata_coculture [metadata_coculture['7AAD'] != "-"]
    metadata_coculture = metadata_coculture [metadata_coculture['Cell type'] != "Jurkat"]
    # Reduce some columns & fill nan's
    metadata_coculture = metadata_coculture.drop(['Unnamed: 0','Timepoint', 'Med+Fil', 'Stimulation', '7AAD'], axis=1)
    metadata_coculture = metadata_coculture.fillna("-")
    # Plots! - careful with the plt.rcParams["figure.figsize"]
    plot_patientVariability(metadata_coculture)
    plot_characteristics_detailed(metadata_coculture)
    plot_migration(metadata_coculture)
    plot_characteristics_condensed(metadata_coculture)
    ## Long term culture X co-culture 
    os.chdir(dir)   # Change the working directory
    metadata_LTculture =pd.read_excel("Post-stimulation culture (RPMIx2) Experiment Summary.xlsx", index_col=None)
    # Exclude Jurkat and 7AAD-
    metadata_LTculture = metadata_LTculture [metadata_LTculture['7AAD'] != "-"]
    metadata_LTculture = metadata_LTculture [metadata_LTculture['Cell type'] != "Jurkat"]
    # Reduce some columns & fill nan's
    metadata_LTculture = metadata_LTculture.drop(['Unnamed: 0','Timepoint', 'Med+Fil', 'Stimulation', '7AAD', 'Insert type'], axis=1)
    metadata_LTculture = metadata_LTculture.fillna("-")
    # Fix the days and media data
    metadata_LTculture.loc[metadata_LTculture['Media'].eq('RPMI'), 'Media'] = "LG"
    metadata_LTculture.loc[metadata_LTculture['Folder tag'].eq(2), 'Folder tag'] = 7
    metadata_LTculture.loc[metadata_LTculture['Folder tag'].eq(3), 'Folder tag'] = 9
    metadata_LTculture.loc[metadata_LTculture['Folder tag'].eq(4), 'Folder tag'] = 14
    metadata_LTculture.loc[metadata_LTculture['Folder tag'].eq(5), 'Folder tag'] = 16
    metadata_LTculture.loc[metadata_LTculture['Folder tag'].eq(1), 'Folder tag'] = 4
    # Plotting - careful with the plt.rcParams["figure.figsize"]
    plot_LTculture(metadata_LTculture, metadata_coculture)
    ## Restimulation data
    os.chdir(dir)   # Change the working directory
    metadata_restimulation =pd.read_excel("Reinduction Experiment Summary.xlsx", index_col=None)
    # Exclude Jurkat and 7AAD-
    metadata_restimulation = metadata_restimulation [metadata_restimulation['7AAD'] != "-"]
    metadata_restimulation = metadata_restimulation [metadata_restimulation['Cell type'] != "Jurkat"]
    # Reduce some columns & fill nan's
    metadata_restimulation = metadata_restimulation.drop(['Unnamed: 0','Timepoint', 'Stimulation', '7AAD'], axis=1)
    metadata_restimulation = metadata_restimulation.fillna("-")
    # Plot! Careful with the plt.rcParams["figure.figsize"]
    metadata = plot_restimulation (metadata_restimulation)