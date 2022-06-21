import pandas as pd
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

def merge_three_channels (image_names, verbose = False):
    plt.rcParams["figure.figsize"] = (20, 20)
    merger = []
    for i in range(len(image_names)):
        input_image = Image.open(image_names[i])
        input_pixels = input_image.load()
        width = input_image.width
        height = input_image.height
        # Split the image into three channels
        monocolor = np.empty((width, height))
        summa = 0
        counter = 0
        for x in range(width):
            for y in range(height):
                pixel = input_pixels[x, y]
                monocolor[x, y] = pixel[i] / 3
                if pixel [i] > 0:
                    summa = summa + pixel [i]
                    counter = counter + 1
        # Modify intensities based on average values
        average = summa / counter 
        intensity_scales = [200, 200, 100]
        monocolor = monocolor*(intensity_scales[i]/average)
        # Stretch if max too low
        max_pixel = [175, 175, 150]
        if np.max(monocolor) < max_pixel[i]:
            print (np.max(monocolor))
            monocolor = monocolor*(max_pixel[i]/np.max(monocolor))
        # Remove oversaturated values (max<255)
        for x in range(width):
            for y in range(height):
                if monocolor[x, y] > 225:
                    monocolor[x, y] = 225
        # Optional visualisation                 
        if verbose:
            plt.imshow(monocolor, cmap='gray')
            plt.title([i])
            plt.show()
        # Preserve the result
        merger.append(monocolor)
    # Stack RGB and export
    rgb = (np.dstack((merger[0],merger[1],merger[2]))).astype(np.uint8)
    plt.imshow(rgb)
    plt.axis('off')  
    if verbose:
        plt.show()
    return

# Read return metadata from a file name (e.g. media type or microscope magnification)
def metadata_read (file_name):
    #  Normal or obesogenic?
    if 'norm' in file_name:
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
    elif 'x20' in file_name:
        magnification = 'x20'
    # Channel name?
    if 'Bodi' in file_name:
        channel = 'Bodipy'
    elif 'dapi' in file_name:
        channel = 'Hoechst'
    elif 'Mito' in file_name:
        channel = 'MitoTracker'
    # Replicate?
    replicate = file_name [-5]
    while replicate == ' ':
        print ('Error!')
        import time
        time.sleep(1) # Sleep for 3 seconds 
    # Donor?
    donor = file_name[5]
    # Collate and return the data
    df = pd.DataFrame ([[file_name, lipid, ASC, channel, magnification, 
                         replicate, donor]],
                       columns=('File Name','Media type', 'Ascorbic acid',
                                'Channel', 'Magnificiation', 'Replicate', 
                                'Donor'))
    return df


def bulk_analysis (dir, folder_input, folder_output = r'\Output folder directory', test = False):
    #♠ Gather and order metadata of files across all folders
    metadata = pd.DataFrame()
    for folder in folder_input:
        os.chdir(dir + folder)   # Change the working directory
        file_names = os.listdir()    # Check the files in the directory
        for file in file_names:
            metadata = metadata.append(metadata_read(file))
    # Re-arrange the data to keep replicates near one another
    metadata = metadata.sort_values(by=['Donor','Media type', 'Ascorbic acid','Magnificiation',
                                        'Replicate','Channel'])
    # Merge it all!
    n = 3
    for i in range ( len (metadata)//n):
        # Subset metadata
        k = i*n
        metadata_temp = (metadata [k:k+n])
        # Move to the correct folder
        if test :
            os.chdir(dir + folder_input[0])
            folder_output = folder_input[0]
        else:
            indexer = 1
            os.chdir(dir + folder_input[int(metadata_temp.iloc[0]['Donor'])-indexer])      
        # Subset file names
        temp_1 = metadata_temp [metadata_temp['Channel'] == 'MitoTracker']
        temp_2 = metadata_temp [metadata_temp['Channel'] == 'Bodipy']
        temp_3 = metadata_temp [metadata_temp['Channel'] == 'Hoechst']
        image_names = [temp_1.iloc[0]['File Name'],
                       temp_2.iloc[0]['File Name'],
                       temp_3.iloc[0]['File Name']]
        plt.close('all') # Close any open plots to avoid a mess
        merge_three_channels(image_names)
        os.chdir(dir + folder_output) # Where to save the figure
        fig_name = metadata_temp.iloc[0]['File Name'][0:10] + metadata_temp.iloc[0]['Ascorbic acid'] + ' ' + metadata_temp.iloc[0]['Media type'] + ' ' + metadata_temp.iloc[0]['Magnificiation'] + ' '  + metadata_temp.iloc[0]['Replicate']
        if not test:
            plt.savefig(fig_name, transparent=True, bbox_inches='tight')
        else:
            plt.show()
            

if __name__ == "__main__":
    import os
    dir = r'C:\Input folder directory'
    folder_input = ['\Input folder directory #1',
               '\Input folder directory #2',
               '\Input folder directory #3',
               '\Input folder directory #4']
    bulk_analysis (dir, folder_input)
    
    # Subsection to run tests on a smaller isolated image batch
    dir_test = r'C:\Test folder directory'
    folder_input_test = ['\Test\Fluorescence']
    bulk_analysis(dir_test, folder_input_test, test = True)

        
        
    
