import os
import math
import numpy as np
import matplotlib.pyplot as plt
from skimage import io
from skimage.filters import threshold_otsu
from scipy.ndimage.filters import median_filter
from statistics import mean

def bin(array, bin_size):
    processed = []
    z, i, j = array.shape
    for image in array:
        binned_image = image.reshape(i//bin_size, bin_size, j//bin_size, bin_size).mean(-1).mean(1)
        processed.append(binned_image)
    return np.array(processed)

def filter(array, filter_size):
    processed = []
    z, i, j = array.shape
    for image in array:
        filtered_image = median_filter(image, size = filter_size)
        processed.append(filtered_image)
    return np.array(processed)

def plot(array):
    #convenience function for plotting a full Z-stack
    z, i, j = array.shape
    #arrange subplots with one row for each channel
    fig, ax = plt.subplots(
        nrows = 4, ncols = int(z/4),
        sharex='col', sharey='row',
        gridspec_kw={'hspace': 0, 'wspace': 0}
        )
    xi_max = z/4
    xi = 0
    for i, image in enumerate(array):
        #set which row of the sublot
        yi = math.floor(i/xi_max)
        ax[yi, xi].imshow(image)
        #iterate through each column of the subplot
        xi += 1
        if xi == xi_max:
            xi = 0
    plt.show()

def extract_membrane_mask(array):
    #retrieve the cell membrane mask channel only
    z, i, j = array.shape
    channel_length = int(z/4)
    membrane_channel = array[3*channel_length:z,:,:]
    #set the threshold value with otsu's algorithm
    threshold = threshold_otsu(membrane_channel)
    thresholded_membrane_channel = []
    for i, image in enumerate(membrane_channel):
        #only keep intensity values above the threshold
        thresholded_image = image > threshold
        thresholded_membrane_channel.append(thresholded_image)
    return(np.tile(thresholded_membrane_channel, (4,1,1)))

def pcc(array1, array2):
    #calculate Pearsons correlation coefficient
    array1_mean = np.mean(array1)
    array2_mean = np.mean(array2)
    r = np.sum((array1 - array1_mean) * (array2 - array2_mean)) / np.sqrt(np.sum((array1 - array1_mean) ** 2) * (np.sum((array2 - array2_mean) ** 2)))
    return(r)

def plot_intensities(array):
    z, i, j = array.shape
    channel_length = int(z/4)
    gfp_images = array[channel_length:2*channel_length,:,:]
    anap_images = array[2*channel_length:3*channel_length,:,:]
    gfp_intensities = gfp_images[gfp_images>0].flatten()
    anap_intensities = anap_images[anap_images>0].flatten()
    plt.plot(gfp_intensities, anap_intensities, ',')
    plt.show()
    return(pcc(gfp_intensities, anap_intensities))


#read and bin the image stacks
a1 = filter(bin(io.imread("data/confocal_images/WT-GFP1.tif"), bin_size=2), filter_size=10)
a2 = filter(bin(io.imread("data/confocal_images/WT-GFP2.tif"), bin_size=2), filter_size=10)
a3 = filter(bin(io.imread("data/confocal_images/WT-GFP3.tif"), bin_size=2), filter_size=10)

b1 = filter(bin(io.imread("data/confocal_images/W311-GFP1.tif"), bin_size=2), filter_size=10)
b2 = filter(bin(io.imread("data/confocal_images/W311-GFP2.tif"), bin_size=2), filter_size=10)
b3 = filter(bin(io.imread("data/confocal_images/W311-GFP3.tif"), bin_size=2), filter_size=10)

plot(a1)
plot(b1)

#invert the binary matrix in order to set all values below the threshold (not membrane) to zero
a1[np.invert(extract_membrane_mask(a1))] = 0
a2[np.invert(extract_membrane_mask(a2))] = 0
a3[np.invert(extract_membrane_mask(a3))] = 0
b1[np.invert(extract_membrane_mask(b1))] = 0
b2[np.invert(extract_membrane_mask(b2))] = 0
b3[np.invert(extract_membrane_mask(b3))] = 0

a1_pcc = plot_intensities(a1)
a2_pcc = plot_intensities(a2)
a3_pcc = plot_intensities(a3)
b1_pcc = plot_intensities(b1)
b2_pcc = plot_intensities(b2)
b3_pcc = plot_intensities(b3)

print(mean([a1_pcc, a2_pcc, a3_pcc]), mean([b1_pcc, b2_pcc, b3_pcc]))