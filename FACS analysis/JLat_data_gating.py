import os
import numpy as np
import pandas as pd

# Several ways metadata was added to the file name (e.g. cell type or stimulus) 
# The first two naming conventions (NC) are identical if there are no replicates
# Separate NC3 and NC4 were made for co-culture and reinduction experiments 

def metadata_read_NC0 (file_name):
    # What cell type?
    if 'Jurkat' in file_name:
        cell_type = 'Jurkat'
    if 'J-LAT' in file_name:
        cell_type = 'J-LAT'
    # 7-AAD or PBS?
    if 'sans' in file_name:
        AAD = '-'
    else:
        AAD = '+'
    # What stimuli?
    costimuli = ['Iono', 'OXA', 'JQ1', 'RVX']
    stimulation = 'N/A'
    if 'PMA' in file_name:
        stimulation = 'PMA'
        for costimulus in costimuli:
            if costimulus in file_name:
                stimulation = stimulation + '+' + costimulus
    else:
        for costimulus in costimuli:
            if costimulus in file_name:
                stimulation = costimulus
    # What timepoint?
    timepoints = ['0h', '1h', '2h', '4h', '8h', '16h', '24h', '30h', '48h']
    timepoint = ''
    for time in timepoints:
        if time in file_name:
            timepoint = time
    # What replicate?
    if stimulation == 'PMA+Iono': 
        stimulation = 'PMA Iono '
    if stimulation == 'N/A': 
        stimulation = 'CTL'
    replicates = [stimulation+'1', stimulation+'2', stimulation+'3']
    replicate = 1
    for i in range (3):
        if replicates [i] in file_name:
            replicate = i+1
    # Restore the Stimulation title
    if stimulation == 'PMA Iono': 
        stimulation = 'PMA+Iono '
    # Caspase inhibitor?
    iCasp = 'no iCasp'
    if 'Caspi' in file_name:
        iCasp = 'iCasp'
    # Collate and return the data
    df = pd.DataFrame ([[cell_type, stimulation, AAD, timepoint, replicate, iCasp]],
                       columns=('Cell type', 'Stimulus', '7AAD', 'Timepoint', 
                                'Replicate', 'iCasp'))
    return df

def metadata_read_NC1 (file_name):
    # What cell type?
    if 'Jurkat' in file_name:
        cell_type = 'Jurkat'
    if 'J-LAT' in file_name:
        cell_type = 'J-LAT'   
    # 7-AAD or PBS?
    if 'sans' in file_name:
        AAD = '-'
    else:
        AAD = '+'      
    # What stimuli?
    costimuli = ['Iono', 'OXA', 'JQ1', 'RVX']
    stimulation = 'N/A'
    if 'PMA' in file_name:
        stimulation = 'PMA'+ '-' 
        for costimulus in costimuli:
            if costimulus in file_name:
                stimulation = stimulation + costimulus  + '-' 
    else:
        for costimulus in costimuli:
            if costimulus in file_name:
                stimulation = costimulus + '-' 
    # What timepoint?
    timepoint = ''
    # What replicate?
    if stimulation == 'N/A': 
        stimulation = 'CTL-'
    replicates = [stimulation+'1', stimulation+'2', stimulation+'3']
    replicate = 1
    for i in range (3):
        if replicates [i] in file_name:
            replicate = i+1
    # Remove the tailing '-'
    if  '-' in stimulation: 
        stimulation = stimulation[0:len(stimulation)-1]
    # Caspase inhibitor?
    iCasp = 'no iCasp'
    if 'Caspi' in file_name:
        iCasp = 'iCasp'
    # Collate and return the data
    df = pd.DataFrame ([[cell_type, stimulation, AAD, timepoint, replicate, iCasp]],
                       columns=('Cell type', 'Stimulus', '7AAD', 'Timepoint', 
                                'Replicate', 'iCasp'))
    return df

def metadata_read_NC2 (file_name):
    # Specifically for co-culture experiments
    cell_type = 'N/A'
    timepoint = ''
    #  What cell type?
    if 'Jurkat' in file_name:
        cell_type = 'Jurkat'
    if 'J-LAT' in file_name:
        cell_type = 'J-LAT'   
    # 7-AAD or PBS?
    if 'sans' in file_name:
        AAD = '-'
    else:
        AAD = '+'      
    # What stimuli?
    stimulation = 'Error'
    if 'PMA' in file_name:
        stimulation = 'PMA-JQ1-' 
    elif 'CTL' in file_name:
        stimulation = 'CTL-'
    # What replicate?
    replicates = [stimulation+'1', stimulation+'2', stimulation+'3']
    replicate = 0
    for i in range (3):
        if replicates [i] in file_name:
            replicate = i+1
    # Remove the tailing '-'
    if  '-' in stimulation: 
        stimulation = stimulation[0:len(stimulation)-1]
    # Filter size
    filter_size = ''
    if ('insert5' in file_name) or (' 5' in file_name):
        filter_size = '5 um'
    elif '04 ' in file_name:
        filter_size = '0.4 um'
    # Normogenic or adipogenic 
    lipid = 'RPMI'
    if 'LG' in file_name:
        lipid = 'LG'
    elif 'HG' in file_name:
        lipid = 'HG'
    # what patient? 
    patient = ''
    if 'P1' in file_name:
        patient = 'P1'
    elif 'P2' in file_name:
        patient = 'P2'
    elif 'P3' in file_name:
         patient = 'P3'
    elif 'P4' in file_name:
          patient = 'P4'
    # Some pooled metadata for graphing/
    meta_pool = lipid + ' '+ patient +'(' + filter_size + ')'
    # Collate and return the data
    df = pd.DataFrame ([[cell_type, stimulation, AAD, timepoint, patient, replicate, 
                         lipid, filter_size, meta_pool]],
                       columns=('Cell type', 'Stimulus', '7AAD','Timepoint', 
                                'Donor','Replicate', 'Media', 'Insert type', 'Med+Fil'))
    return df

def metadata_read_NC3 (file_name):
    # Specifically for reinduction experiments
    cell_type = 'N/A'
    timepoint = ''
    #  What cell type?
    if 'Jurkat' in file_name:
        cell_type = 'Jurkat'
    if 'J-LAT' in file_name:
        cell_type = 'J-LAT'   
    # 7-AAD or PBS?
    if 'sans' in file_name:
        AAD = '-'
    else:
        AAD = '+'      
    # What stimuli?
    stimulation = "Error"
    if 'CTL-PMA' in file_name:
        stimulation = 'CTL-PMA-'
    elif 'PMA-CTL' in file_name:
        stimulation = 'PMA-CTL-'
    elif 'PMA-PMA' in file_name:
        stimulation = 'PMA-PMA-'
    else: 
        if 'PMA' in file_name:
            stimulation = 'PMA-'
        if 'CTL' in file_name:
            stimulation = 'CTL-'
    # What replicate?
    replicates = [stimulation+'1', stimulation+'2', stimulation+'3']
    replicate = 1
    for i in range (3):
        if replicates [i] in file_name:
            replicate = i+1
    # Remove the tailing '-'
    if  '-' in stimulation: 
        stimulation = stimulation[0:len(stimulation)-1]
    # Collate and return the data
    df = pd.DataFrame ([[cell_type, stimulation, AAD, timepoint, replicate]],
                       columns=('Cell type', 'Stimulus', '7AAD','Timepoint', 
                                'Replicate'))
    return df


# Function that gates a single sample (.fcs file)
def JLAT_gating (file_name, verbose=False, export=False, figure_size=5,
                     point_size=1, alpha=0.2, iteration=0, naming_convention=2):
    import FlowCytometryTools # https://eyurtsev.github.io/FlowCytometryTools/tutorial.html
    import matplotlib.pyplot as plot
    # Initialise 4 subplots and set figure size
    plot.rcParams["figure.figsize"] = (4*figure_size,figure_size)
    fig, (ax1, ax2, ax3, ax4) = plot.subplots(1, 4, constrained_layout = True)    
    # Generate the major plot title 
    if naming_convention==0:
        df = metadata_read_NC0 (file_name)
    elif naming_convention==1:
        df = metadata_read_NC1 (file_name)
    elif naming_convention==2:
        df = metadata_read_NC2 (file_name)
    elif naming_convention==3:
        df = metadata_read_NC3 (file_name)
    # Setup some plot aparameters 
    plot_title = df['Cell type'] + ' & ' + df['Stimulus']  + ' '
    plot_title = plot_title + df['Timepoint'] + ' (7AAD' + df['7AAD'] + ')'
    plot_title = plot_title.to_string (index=False)
    plot.suptitle(plot_title, fontsize=30, fontweight='roman')
    # On-graph text parameters
    font_gates = {'family': 'sans-serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 20}
    # Plotting parameters
    channel = ['FSC-A','FSC-H','SSC-A','7AAD-A','GFP-A']
    channel_pairs = [['FSC-A','FSC-H'],['FSC-A','SSC-A'],
                     ['GFP-A', '7AAD-A'],['GFP-A', 'FSC-A']]
    # Load the sample datafile
    sample = FlowCytometryTools.FCMeasurement(ID='Test Sample', datafile=file_name)
    sample_0 = sample.transform('hlog', channels=channel, b=500.0)
    # FSC-A vs. FSC-H: gating single cells pt. 1
    x = 8190
    dx = 1800
    y = 8000-50
    dy = 1750
    ddy = 125
    gate_xy = [(x, y),(x+dx, y+dy),(x+dx, y+dy+ddy),(x, y+ddy)]
    single_cell_1 = FlowCytometryTools.PolyGate (vert = gate_xy,
                                               channels=channel_pairs[0], 
                                               name='Single cells')
    # Gate the sample and record the % of gated events
    sample_1 = sample_0.gate(single_cell_1)
    fraction_1 = sample_1.counts/sample_0.counts*100
    if verbose or export:
        # Subplot index key: 1 roz, 4 columns, graph index = 1
        ax1 = plot.subplot(141)
        # Fix axis across all samples
        ax1.set_xlim([7900,10000])
        ax1.set_ylim([7900,10000])
        # Plot the dataset
        sample_0.plot(channel_pairs[0], kind='scatter',
                      s=point_size, alpha=alpha, color='black')
        # Plot the gate and add a description
        single_cell_1.plot(color='red') 
        plot.text(x-90, y+700,'Singlets\n'+str(round(fraction_1, 2))+ 
                  '%', fontdict=font_gates)
        # Visual parameter tweaking
        plot.xticks(np.arange(8000, 10000+1, 500))
        ax1.tick_params(axis="x", labelsize=14)
        ax1.xaxis.label.set_size(20)
        plot.yticks(np.arange(8000, 10000+1, 500))
        ax1.tick_params(axis="y", labelsize=14)
        ax1.yaxis.label.set_size(20)
    
    # FSC-A vs. SSC-A: gating single cells pt.2
    import math
    x = 9200
    y = 8050
    r = 600*1.3 # Radius of the pseudo-eliptic gate
    f = 1.3 # "Elipsisity" factor (horizontal)
    gate_xy = [(x+r/f, y),(x+r*math.sqrt(1/2)/f, y+r*math.sqrt(1/2)),
               (x, y+r),(x-r*math.sqrt(1/2)/f, y+r*math.sqrt(1/2)),
               (x-r/f, y),(x-r*math.sqrt(1/2)/f, y-r*math.sqrt(1/2)),
               (x, y-r),(x+r*math.sqrt(1/2)/f, y-r*math.sqrt(1/2))]
    single_cell_2 = FlowCytometryTools.PolyGate (vert = gate_xy,
                                               channels=channel_pairs[1], 
                                               name='Single cells')
    sample_2 = sample_1.gate(single_cell_2)
    fraction_2 = sample_2.counts/sample_1.counts*100
    if verbose or export:
        ax2 = plot.subplot(142)
        ax2.set_xlim([7900,10000])
        ax2.set_ylim([6350,10000])
        sample_1.plot(channel_pairs[1], kind='scatter',
                    s=point_size, alpha=alpha, color='orange')
        single_cell_2.plot(color='red')
        plot.text(x-200, y-1300,'Granularity\n'+str(round(fraction_2, 2))+ 
                  '%', fontdict=font_gates)
        # Visual parameter tweaking
        plot.xticks(np.arange(8000, 10000+1, 500))
        ax2.tick_params(axis="x", labelsize=14)
        ax2.xaxis.label.set_size(20)
        plot.yticks(np.arange(6500, 10000+1, 700))
        ax2.tick_params(axis="y", labelsize=14)
        ax2.yaxis.label.set_size(20)
        
    # GFP-A vs. 7AAD-A: gating live cells
    x = 3500
    y = 2000
    live_cell = FlowCytometryTools.ThresholdGate(y, ['7AAD-A'], region='below')
    sample_3 = sample_2.gate(live_cell)
    fraction_3 = sample_3.counts/sample_2.counts*100
    if verbose or export:
        ax3 = plot.subplot(143)  
        ax3.set_xlim([-500,10000])
        ax3.set_ylim([-2300,10000])
        sample_2.plot(channel_pairs[2], kind='scatter', gates=[live_cell],
                    s=point_size, alpha=alpha, color='blue')
        plot.text(x+2200, y+600,'Live cells\n'+str(round(fraction_3, 2))+ 
                  '%', fontdict=font_gates)
        # Visual parameter tweaking
        plot.xticks(np.arange(0, 10000+1, 2000))
        ax3.tick_params(axis="x", labelsize=14)
        ax3.xaxis.label.set_size(20)
        plot.yticks(np.arange(-2000, 10000+1, 2000))
        ax3.tick_params(axis="y", labelsize=14)
        ax3.yaxis.label.set_size(20)

    # GFP-A vs. FSC-A: gating GFP+ cells
    x = 1000
    y = 9210
    GFP_cell = FlowCytometryTools.ThresholdGate(1000.0, ['GFP-A'], region='above')
    sample_4 = sample_3.gate(GFP_cell)
    fraction_4 = sample_4.counts/sample_3.counts*100
    if verbose or export:
        ax4 = plot.subplot(144)
        ax4.set_xlim([-500,10000])
        ax4.set_ylim([7900,10000])
        sample_3.plot(channel_pairs[3], kind='scatter', gates=[GFP_cell],
                      s=point_size, alpha=alpha, color='green')
        plot.text(x+700, y-1100,'GFP+ cells\n'+str(round(fraction_4, 2))+ 
                  '%', fontdict=font_gates)
        # Visual parameter tweaking
        plot.xticks(np.arange(0, 10000+1, 2000))
        ax4.tick_params(axis="x", labelsize=14)
        ax4.xaxis.label.set_size(20)
        plot.yticks(np.arange(8000, 10000+1, 500))
        ax4.tick_params(axis="y", labelsize=14)
        ax4.yaxis.label.set_size(20)
    
    # Export the analysis results
    if export:
        figure_name = file_name [0:(len(file_name)-4)]
        plot.savefig(str(figure_name) + '.png') # Save the final figure
    if verbose:
        plot.show()
    plot.close('all') #https://stackoverflow.com/questions/24500065/closing-matplotlib-figures
    
    # Return metadata 
    df.insert (df.shape[1], '% Single cells [1]', round(fraction_1, 2))
    df.insert (df.shape[1], '% Single cells [2]', round (fraction_2, 2))
    df.insert (df.shape[1], '% Live cells', round (fraction_3, 2))
    df.insert (df.shape[1], '% GFP+ cells', round (fraction_4, 2))
    # Absolute event numbers:
    df.insert (df.shape[1], 'Total Single cells [1]', sample_1.counts)
    df.insert (df.shape[1], 'Total Single cells [2]', sample_2.counts)
    df.insert (df.shape[1], 'Total Live Cells', sample_3.counts)
    df.insert (df.shape[1], 'Total GFP+ cells', sample_4.counts)
    # MFI:
    df.insert (df.shape[1], 'MFI GFP+', sample_3['GFP-A'].median())
    return df

def excel_export (metadata):
    # Re-arrange the data to keep replicates near one another
    metadata = metadata.sort_values(by=['7AAD', 'Cell type', 'Stimulus', 'Timepoint', 'Replicate'])
    # Drop extra columns if empty
    if 'Timepoint' in metadata:
        if metadata['Timepoint'].sum()=="":
            metadata = metadata.loc[:, metadata.columns!='Timepoint']
    if 'Donor' in metadata:
        if metadata['Donor'].sum()=="":
            metadata = metadata.loc[:, metadata.columns!='Donor']
    if 'Insert type' in metadata:
        if metadata['Insert type'].sum()=="":
            metadata = metadata.loc[:, metadata.columns!='Insert type']
    # Skip all analysis if there are no replicates 
    n = metadata.nlargest(1, ['Replicate'])['Replicate']
    n.reset_index(drop=True, inplace=True)
    n = n[0]
    if n == 1:
        metadata = metadata.loc[:, metadata.columns!='Replicate']
        metadata.to_excel(r'Analysis Results.xlsx')
    metadata.to_excel(r'Analysis Results.xlsx')
    return metadata

# References to data that we need to plot
data_references = ['% Single cells [1]','% Single cells [2]','% Live cells',
                   '% GFP+ cells','Total Single cells [1]','Total Single cells [2]',
                   'Total Live Cells','Total GFP+ cells', 'MFI GFP+']

def plot_export (metadata, figure_size=4, Jurkat=True, rotation=90, hue=['Stimulation',''], folder=''):
    # hue = [X-axis label, Legend labels]
    # Start by preparing the plot template etc.
    import seaborn as sns
    import matplotlib.pyplot as plt
    plt.close('all') # Close any open plots to avoid a mess
    plt.rcParams["figure.figsize"] = (figure_size*5,figure_size*2.21)
    # Prepare 4 subplot axes & appropriate subplo titles
    fig, ax = plt.subplots(1, 10, constrained_layout = True)
    
    subplots = data_references.copy()
    subplots.insert(4, '% GFP+ cells')
    # '% GFP+ cells' appear twice because there are two graphs with 
    # different axis spans
    subplot_titles = ['Singlets', 'Granularity', 'Live cells', 'GFP+ cells', 
                      'GFP+ cells',
                      'Singlets', 'Granularity', 'Live cells', 'GFP+ cells', 
                      'Live cells']
    
    # Sort dataframe to ensure proper order on the barplot
    metadata ['Stimulation'] = metadata.loc [:, "Cell type"] + ' * ' + metadata.loc [:, "Stimulus"]
    
    # Sort by stimulus (and by timepoint if there are applicable)
    if ('Timepoint' in hue):  
        plt.rcParams["figure.figsize"] = (figure_size*7,figure_size*2.21)
        rotation = 0
        timepoints = ['0h', '1h', '2h', '4h', '8h', '16h', '24h', '30h', '48h']
        metadata ['Timepoint_tag'] = metadata ['Timepoint']
        metadata ['Timepoint_tag'] = metadata['Timepoint_tag'].replace(to_replace=timepoints, 
                                                                       value=range(len(timepoints)))
        metadata = metadata.sort_values(by=['Stimulation', 'Timepoint_tag'])
    else:
        metadata = metadata.sort_values(by=['Stimulation']) 
    
    # Drop 7AAD(-) datapoints from the plot
    metadata = metadata [metadata['7AAD'] != "-"]
    # Optionally exclude Jurkat
    if Jurkat != True:
        metadata = metadata [metadata['Cell type'] != "Jurkat"]
        metadata ['Stimulation'] = metadata.loc [:, "Stimulus"] # No need for 'J-Lat * ' if no Jurkat
    
    # Adjust rotation if there are few x ticks
    n = len(pd.unique(metadata['Stimulation']))
    if (n <= 6) and ('Timepoint' not in hue):
        rotation = 45

    # Cycle through 4 gates and draw the subplots
    for i in range(len(ax)):
       
        # Index the current plot
        ax[i]= plt.subplot(2, 5, i+1)
        
        # Refine the plot -> appropriate legends and plots 
        if len(hue[1]) != 0:
            ax[i]= sns.barplot (x = hue[0], y=subplots[i],  hue=hue[1], data=metadata)
        else:
            ax[i]= sns.barplot (x = hue[0], y=subplots[i],  data=metadata)
        ax[i].tick_params(axis='x', rotation=rotation)
        if i < 4:
            ax[i].set_ylim([0,105])
        elif i == 4:
            ax[i].set_ylim([0,52])
        plt.grid(axis = 'y', color = 'gray', linestyle = '--', linewidth = 0.5)
        plt.title(subplot_titles[i])
        
        # Clean up the plot by removing excess labels 
        plt.xlabel("")
        if (i != 4) and (hue[1] != ''): # if hue[1] == '', there is no legend!
            ax[i].get_legend().remove()
        if i == 0:
            plt.ylabel("Fraction of events (%)")
            #ax[i].axes.yaxis.set_ticklabels([])
        elif i == 5:
            plt.ylabel("Number of events")
        elif i == 9:
            plt.ylabel("Mean Fluorescence Intensity")
            ax[i].set_ylim(bottom = 250)
        else:
            plt.ylabel("")
        # Introduce logarithmic y-axes for absolute cell count graphs 
        # symlog also?  Logit doesn't work
        if i in range (5, 9):
            plt.yscale('symlog')

    # Export the resulting p^lot
    if '\\' not in folder[1:]:
        fig_name = folder[1:] + ' Analysis Results '
    else:
        fig_name = folder[24:] + ' Analysis Results '
    if Jurkat == True:
        fig_name = fig_name + '(with Jurkat).png'
        plt.savefig(fig_name, bbox_inches='tight')
    else:
        fig_name = fig_name + '(without Jurkat).png'
        plt.savefig(fig_name, bbox_inches='tight')



def folder_analysis (metadata, hue=['Stimulation',''], folder = ''):
    metadata_means = excel_export(metadata)
    metadata_means = metadata
    plot_export (metadata, Jurkat=True, hue=hue, folder=folder)
    plot_export (metadata, Jurkat=False, hue=hue, folder = folder)
    return metadata_means 
    
# Analysing content within multiple input folders    
def multiple_folders (dir='', folders=[''], verbose=False, export = True, 
                      hue=['Stimulation',''], naming_convention = 0):
    folder_tag = range(len(folders))
    
    # Initiate several empty datasets for downstream manipulation 
    metadata = pd.DataFrame ()
    metadata_temp = pd.DataFrame ()
    metadata_means = pd.DataFrame ()
    
    for k in range(0, len(folders)):
        os.chdir(dir + folders[k])   # Change the working directory
        file_names = os.listdir() # Check the files in the directory
        metadata_temp = pd.DataFrame ()
        for i in range(0, len(file_names)):
            if '.fcs' in file_names [i]:
                df = JLAT_gating(file_names [i], iteration=i, alpha=0.5, 
                                 point_size=2, verbose=verbose, export=export,
                                 naming_convention=naming_convention)
                # The mess below to correctly insert the day tag
                df.insert(loc = 5, column = 'Folder tag', value = folder_tag[k])
                metadata_temp = metadata_temp.append(df, ignore_index=True)
        
        # Collate raw metadata from each folder to return
        metadata = metadata.append(metadata_temp, ignore_index=True)
        # Re-use metadata_temp to a treated data chunk (replicate means and STDs)
        metadata_temp = folder_analysis (metadata_temp.loc[:, metadata_temp.columns!='Folder tag'], 
                       hue=hue, folder = folders[k])
        # Folder tag gets lost along the way - bring it back!
        metadata_temp.insert(loc = 5, column = 'Folder tag', value = folder_tag[k])
        # Collate treated metadata from each folder to return
        metadata_means=metadata_means.append(metadata_temp, ignore_index=True)

    return metadata, metadata_means, df
        
# Metadata for experiments within multiple folders
experiments = ['Experiment name #1', 
               'Experiment name #2', 
               'Experiment name #3',
               'Experiment name #4']

folders = [[r'\Folder directory #1',r'\Folder directory #2',r'\Folder directory #3'],
           [r'\Folder directory #4',r'\Folder directory #5',r'\Folder directory #6'],
           [r'\Folder directory #7',r'\Folder directory #8',r'\Folder directory #9'],
           [r'\Folder directory #10',r'\Folder directory #11',r'\Folder directory #12']]

naming_convnetion_0 =  [n for n in range(0, 1)]     # Relevant to metadata_read_NC0
naming_convnetion_1 =  [n for n in range(1, 2)]     # Relevant to metadata_read_NC1
naming_convnetion_2 =  [n for n in range(2, 3)]    # Relevant to metadata_read_NC2
naming_convnetion_3 = [4]

skips = [n for n in range(0, 11)] 

def analyse_everything_thus_far_v1 (dir = r'C:\FACS folder directory',
                                    experiments=experiments,folders=folders,
                                    export = True, verbose=False):
    for i in range(len(experiments)):
        if i in naming_convnetion_0:
            naming_convention = 0
        elif i in naming_convnetion_1:
            naming_convention = 1
        elif i in naming_convnetion_2:
            naming_convention = 2
        elif i in naming_convnetion_3:
            naming_convention = 3
        if i not in skips:
            folders_subset = folders[i]
            hue=['Stimulation',''] # hue = [X-axis label, Legend labels]
            if 'PMA Kinetics' in experiments[i]:
                hue=['Timepoint','Stimulation']
            if 'Q-VD-OPh' in experiments[i]:
                hue=['Stimulation','iCasp']
            if 'Co-culture' in experiments[i]:
                hue=['Stimulation','Med+Fil']
            if 'RPMIx2' in experiments[i]:
                hue=['Stimulation','Media']
            metadata, metadata_means, df = multiple_folders (dir, folders_subset, verbose=verbose, 
                                                         export = export, hue = hue,
                                                         naming_convention = naming_convention)
            os.chdir(dir)
            # Export the resulting p^lot
            if '\\' not in folders[i][0][1:]:
                spreadsheet_name = r'(' + folders[i][0][1:] + r') '
            else:
                spreadsheet_name =  r'(' + folders[i][0][24:] + r') '
            filename = spreadsheet_name + experiments[i] + ' Experiment Summary.xlsx'
            metadata_means.to_excel(filename)
    return metadata

if __name__ == "__main__":
    # Export = plotting data qnd gates for each .fcs
    metadata = analyse_everything_thus_far_v1 (export = False) 
    
