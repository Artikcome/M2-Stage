"""
This section was not used - the .fcs are compensated already
"""

def custom_compensate(original_sample, f1=0.15, f2=0.32):
    # Copy the original sample
    new_sample = original_sample.copy()
    new_data = new_sample.data
    original_data = original_sample.data

    # Our transformation goes here
    new_data['GFP-A'] = original_data['GFP-A'] - f1 * original_data['7AAD-A']
    new_data['7AAD-A'] = original_data['7AAD-A'] - f2 * original_data['GFP-A']
    new_data = new_data.dropna()  # Removes all NaN entries
    new_sample.data = new_data
    return new_sample

def data_load (files):
    import FlowCytometryTools
    # Read/import the four 'control' data sets
    JLAT_CTL_sans7AAD = FlowCytometryTools.FCMeasurement(ID='Test Sample', datafile=files[1])
    JLAT_CTL_7AAD = FlowCytometryTools.FCMeasurement(ID='Test Sample', datafile=files[0])
    JLAT_PMA_sans7AAD = FlowCytometryTools.FCMeasurement(ID='Test Sample', datafile=files[2])
    JLAT_PMA_7AAD = FlowCytometryTools.FCMeasurement(ID='Test Sample', datafile=files[3])
    # Pack in a single array for script simplicity
    samples = [JLAT_CTL_sans7AAD,JLAT_CTL_7AAD,JLAT_PMA_sans7AAD,JLAT_PMA_7AAD]
    # Transform the data
    for sample in samples:
        sample = sample.transform('hlog') 
    return samples

def plot_4 (sample, channel, celltype):
    import matplotlib.pyplot as plt
    # General plot parameters
    plt.rcParams["figure.figsize"] = (10,10)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4) 
    plt.suptitle(celltype, fontsize=30, fontweight='roman') # Major title settigns
    # Subplot 1: CTL sans 7AAD
    ax1 = plt.subplot(221) # Establish sublot position: 221 = 2 rows, 2 columns, graph index = 1
    sample[0].plot(channel, color='green', alpha=0.7, bins=100) # Plot data from one channel
    plt.title('sans 7AAD') # Minor (subplot) title 
    plt.ylabel("Control")
    # Subplot 2: CTL 7AAD
    ax2 = plt.subplot(222)
    sample[1].plot(channel, color='green', alpha=0.7, bins=100)
    plt.title('7AAD')
    plt.ylabel("")
    # Subplot 1: PMA-Iono sans 7AAD
    ax3 = plt.subplot(223)
    sample[2].plot(channel, color='green', alpha=0.7, bins=100)
    plt.ylabel("PMA + Ionolycin")
    # Subplot 3: PMA-Iono 7AAD
    ax4 = plt.subplot(224)
    sample[3].plot(channel, color='green', alpha=0.7, bins=100)
    plt.ylabel("")
    plt.savefig(celltype + ' - ' + channel +'.png') #Export the figure
    

if __name__ == "__main__":
    import os
    dir = os.getcwd() + r'\Compensation'
    os.chdir(dir)   # Change the working directory
    samples = data_load (os.listdir()) # Load all files in the directory
    # Visualise GFP and 7AAD channel data across 4 relevant conditions 
    # to tune the compensation factors (f1, f2) 
    plot_4(samples, channel='GFP-A', celltype='J-LAT') 
    plot_4(samples, channel='7AAD-A', celltype='J-LAT')
    