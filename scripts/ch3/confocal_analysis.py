import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skimage import io
from skimage.filters import threshold_otsu
from scipy.ndimage.filters import median_filter

def filter(array, filter_size):
    processed = []
    for image in array:
        filtered_image = median_filter(image, size = filter_size)
        processed.append(filtered_image)
    return np.array(processed)

def max_proj(array):
    processed = []
    for stack in np.split(array, 4, axis=0):
        maximum = np.max(stack, axis=0)
        processed.append(maximum)
    return np.stack(processed, axis=0)

def normalise(array):
    processed = []
    for stack in np.split(array, 4, axis=0):
        maximum = np.max(stack.flatten())
        processed.append(stack / maximum)
    return np.concatenate(processed, axis=0)

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
        ax[yi, xi].imshow(image, cmap = "gray")
        ax[yi, xi].axis("off")
        #iterate through each column of the subplot
        xi += 1
        if xi == xi_max:
            xi = 0
    plt.show()

def extract_mask(array, channel):
    #retrieve the cell membrane mask channel only
    z, i, j = array.shape
    channel_length = int(z/4)
    channel_images = array[channel*channel_length:(channel+1)*channel_length,:,:]
    #set the threshold value with otsu's algorithm
    threshold = threshold_otsu(channel_images)
    thresholded_channel_images = []
    for i, image in enumerate(channel_images):
        #create binary mask with values below the threshold as true
        thresholded_image = image < threshold
        thresholded_channel_images.append(thresholded_image)
    return(np.tile(thresholded_channel_images, (4,1,1)))

def pcc(array1, array2):
    #calculate Pearsons correlation coefficient
    array1_mean = np.mean(array1)
    array2_mean = np.mean(array2)
    r = np.sum((array1 - array1_mean) * (array2 - array2_mean)) / np.sqrt(np.sum((array1 - array1_mean) ** 2) * (np.sum((array2 - array2_mean) ** 2)))
    return(r)

def plot_intensities(array, channel_1, channel_2):
    z, i, j = array.shape
    channel_length = int(z/4)
    channel_1_images = array[channel_1*channel_length:(channel_1+1)*channel_length,:,:]
    channel_2_images = array[channel_2*channel_length:(channel_2+1)*channel_length,:,:]
    channel_1_intensities = channel_1_images[channel_1_images!=0].flatten()
    channel_2_intensities = channel_2_images[channel_2_images!=0].flatten()
    plt.plot(channel_1_intensities, channel_2_intensities, ',')
    plt.show()
    return(channel_1_intensities, channel_2_intensities)

def plot_profiles(array, z_position, j1, j2):
    z, i, j = array.shape
    profiles = []
    fig, ax = plt.subplots(ncols = 5)
    for index, channel in enumerate(np.split(array, 4, axis=0)):
        image = channel[z_position,:,:]
        ax[index].imshow(image, cmap = "gray")
        ax[index].axis("off")
        profile = np.mean(image[:,j1:j2], axis=1)
        ax[4].plot((profile - np.min(profile)) / (np.max(profile) - np.min(profile)), label = index)
    ax[4].legend()        
    plt.show()

def read(directory, trim):
    stack = io.imread(directory)
    trimmed = []
    for channel in np.split(stack, 4, axis=0):
        z, i, j = channel.shape
        trimmed_channel = channel[trim:(z-trim),:,:]
        trimmed.append(trimmed_channel)
    return np.concatenate(trimmed, axis=0)

#read the image stacks
a1 = read("data/confocal_images/WT-GFP1.tif", trim=4)[:,:,0:340]
a2 = read("data/confocal_images/WT-GFP2.tif", trim=4)[:,122:450,140:]
a3 = read("data/confocal_images/WT-GFP3.tif", trim=4)[:,40:400,:]
b1 = read("data/confocal_images/W311-GFP1.tif", trim=4)[:,0:380,165:440]
b2 = read("data/confocal_images/W311-GFP2.tif", trim=4)[:,100:,50:]

b2_processed = normalise(filter(b2, filter_size=5))
plot(b2_processed)
plot_profiles(b2_processed, 5, 200, 220)

a3_processed = normalise(filter(a3, filter_size=5))
plot(a3_processed)
plot_profiles(a3_processed, 5, 280, 300)

for stack in (a1, a2, a3, b1, b2):
    filtered_stack = filter(stack, filter_size = 20)
    normalised_stack = normalise(filtered_stack)
    plot(normalised_stack)
    maxproj = max_proj(normalised_stack)
    fig, ax = plt.subplots(4)
    for i, image in enumerate(maxproj):
        ax[i].imshow(image, cmap = "gray")
    plt.show()

#set the fluorescence channels
dic_channel = 0
gfp_channel = 1
anap_channel = 2
membrane_channel = 3

wt = []
for i, stack in enumerate((a1, a2, a3), start=1):
    #extract masks for gfp and the membrane stain such that true should be masked
    gfp_mask = extract_mask(stack, channel=gfp_channel)
    membrane_mask = extract_mask(stack, channel=membrane_channel)
    #values inside either of the masks (i.e. not gfp and membrane stain positive) are set to zero
    stack[np.logical_or(gfp_mask, membrane_mask)] = 0
    #extract per pixel intesities for gfp and anap
    gfp, anap = plot_intensities(stack, channel_1=gfp_channel, channel_2=anap_channel)
    wt.append(pd.DataFrame({
        "construct" : "WT-GFP+SUR1",
        "image_stack" : i,
        "gfp_intensity" : gfp,
        "anap_intensity" : anap
    }))

wt = pd.concat(wt)

w311 = []
for i, stack in enumerate((b1, b2), start=1):
    #extract masks for gfp and the membrane stain such that true should be masked
    gfp_mask = extract_mask(stack, channel=gfp_channel)
    membrane_mask = extract_mask(stack, channel=membrane_channel)
    #values inside either of the masks (i.e. not gfp and membrane stain positive) are set to zero
    stack[np.logical_or(gfp_mask, membrane_mask)] = 0
    #extract per pixel intesities for gfp and anap
    gfp, anap = plot_intensities(stack, channel_1=gfp_channel, channel_2=anap_channel)
    w311.append(pd.DataFrame({
        "construct" : "W311*-GFP+SUR1",
        "image_stack" : i,
        "gfp_intensity" : gfp,
        "anap_intensity" : anap
    }))

w311 = pd.concat(w311)

intensities = wt.append(w311)

# intensities.to_csv("data/confocal_images/per_pixel_intensities.csv")
