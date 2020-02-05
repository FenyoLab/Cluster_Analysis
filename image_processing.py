#!/usr/bin/env python

#    Copyright (C) 2020  Sarah Keegan
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from skimage import io, exposure, measure, color, filter
from scipy import ndimage
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import glob
import mahotas as mh

#contains general functions to help process microscopy images
pixel_saturation = 0
bb_line_width = .5
bb_text_size = 6

pixels_to_nm = 20 # (1 px = 20 nm for reconstructed images)

def SaveFigureAsImage(fileName,fig=None,**kwargs):
    ''' Save a Matplotlib figure as an image without borders or frames.
       Args:
            fileName (str): String that ends in .png etc.

            fig (Matplotlib figure instance): figure you want to save as the image
        Keyword Args:
            orig_size (tuple): width, height of the original image used to maintain 
            aspect ratio.
    '''
    fig.patch.set_alpha(0)
    if kwargs.has_key('orig_size'): # Aspect ratio scaling if required
        w,h = kwargs['orig_size']
        fig_size = fig.get_size_inches()
        w2,h2 = fig_size[0],fig_size[1]
        fig.set_size_inches([(w2/w)*w,(w2/w)*h])
        fig.set_dpi((w2/w)*fig.get_dpi())
    a=fig.gca()
    a.set_frame_on(False)
    a.set_xticks([]); a.set_yticks([])
    plt.axis('off')
    plt.xlim(0,h); plt.ylim(w,0)
    fig.savefig(fileName, transparent=True, bbox_inches='tight', \
        pad_inches=0)
    
def intensity_stretch(image):
    if(pixel_saturation == 0):
        min_in_range = image.min()
        max_in_range = image.max()
    else:
        pixels = image.flatten()
        pixels.sort()
        num_saturated = int((pixel_saturation*(len(pixels)))/2)
        index = len(pixels) - num_saturated - 1
        max_in_range = pixels[index]
        min_in_range = pixels[num_saturated]
    
    return exposure.rescale_intensity(
        image, in_range=(min_in_range,max_in_range), out_range=(0,255))

def filter_area(props, mask, min_area, display_image):
    for prop in props:
        area = prop['Area']*pixels_to_nm**2
        
        if(area < min_area):
            
            for c in prop['Coordinates']:
                #set coordinates of this segment to False in the image mask
                mask[c[0],c[1]] = False 
                
                #set coordinates of this segment to 0 in RGB display image
                display_image[c[0],c[1]] = [0,0,0]

def calculate_overlap(props, overlap_labeled_clusters, overlap_dict):
    #overlap_dict contains, for each cluster, list of labels/areas for each overlapping
    #cluster, and the total area overlapped
    for prop in props: 
        coords = prop['Coordinates']
        label = prop['Label']
        area = prop['Area']
        total_area = 0
        cur_overlap_dict = {}
        for c in coords:
            if overlap_labeled_clusters[c[0],c[1]] > 0: #found overlap at point c[0],c[1]: record the label of the overlapping
                                                        #cluster and start/increment the area counter for that label
                if overlap_labeled_clusters[c[0],c[1]] in cur_overlap_dict:
                    cur_overlap_dict[overlap_labeled_clusters[c[0],c[1]]] += 1
                else: cur_overlap_dict[overlap_labeled_clusters[c[0],c[1]]] = 1
                total_area += 1
        ov_list = []
        for ov_label in cur_overlap_dict.keys():
            ov_list.append([ov_label, cur_overlap_dict[ov_label]])
        overlap_dict[label] = [ov_list, total_area]
        
def sort_by_area(props):
    #return list of sorted indexes into the props array
    #the indexes indicate the order of props sorted by largest to smallest c1 area
    prop_positions = []
    for i in range(0,len(props)):
        prop_positions.append([i, props[i]['Area']])
    return sorted(prop_positions, key=lambda p: p[1], reverse=True)

def get_imagej_measures(image, area): #input image is the image of the cluster ( image = prop['Image'] )

    #surround image with 1x1 blank rectangle so entire perimeter will be found
    new_row = np.zeros(len(image[0]))
    image = np.insert(image, 0, new_row, axis=0)
    image = np.insert(image, len(image), new_row, axis=0)
    image = np.insert(image, 0, 0, axis=1)
    image = np.insert(image, len(image[0]), 0, axis=1)
    
    #find perimeter of the sliced image
    image_perim = mh.labeled.bwperim(image)
    
    #remove fake border
    image_perim = np.delete(image_perim, 0, axis=0)
    image_perim = np.delete(image_perim, len(image_perim)-1, axis=0)
    image_perim = np.delete(image_perim, 0, axis=1)
    image_perim = np.delete(image_perim, len(image_perim[0])-1, axis=1)
    
    image = np.delete(image, 0, axis=0)
    image = np.delete(image, len(image)-1, axis=0)
    image = np.delete(image, 0, axis=1)
    image = np.delete(image, len(image[0])-1, axis=1)
    
    #find perimeter manually, trying to match ImageJ perimater (and therefore circulatiry)
    perim_length = 0
    corner_count = 0
    #count perimeter pixels and count corners
    for i,val1 in enumerate(image_perim):
        for j,val2 in enumerate(val1):
            if(val2):
                num_open_edges1 = 0
                num_open_edges2 = 0
                if(i == (len(image)-1) or not image[i+1,j]):
                    perim_length += 1
                    num_open_edges1 += 1
                if(j == (len(val1)-1) or not image[i,j+1]):
                    perim_length += 1
                    num_open_edges2 += 1
                if(i == 0 or not image[i-1,j]):
                    perim_length += 1
                    num_open_edges1 += 1
                if(j == 0 or not image[i,j-1]):
                    perim_length += 1
                    num_open_edges2 += 1
                num_open_edges = num_open_edges1 + num_open_edges2    
                if(num_open_edges == 3 or num_open_edges == 4): corner_count += 2
                elif(num_open_edges == 2 and num_open_edges1 == 1):
                    corner_count += 1
    manual_perim = perim_length - corner_count*(2-(2**.5))
    manual_perim *= 20
    manual_circ = 4 * math.pi * ((area - .5*manual_perim) / manual_perim**2)
    if(manual_circ > 1.): manual_circ = 1.
    
    return (manual_perim, manual_circ)
     
def save_to_csv(file_name, props, overlap_dict, draw_on_figure):
    try: 
        out_file = open(file_name, 'w') 
    except(IOError):
        print 'Error in opening/writing: ' + '"' + file_name + '"' + '.\n'
        return 0
        
    if(overlap_dict == {}):
        out_file.write(',Area,CentroidX,CentroidY,MajorAxisLength,MinorAxisLength,Perimeter,Eccentricity,MaxIntensity,MeanIntensity,MinIntensity\n')
    else:
        out_file.write(',Area,CentroidX,CentroidY,MajorAxisLength,MinorAxisLength,Perimeter,Eccentricity,OverlapSegmentList,OverlapAreaList,TotalOverlapArea,MaxIntensity,MeanIntensity,MinIntensity\n')
    for prop in props: 
        label = prop['Label']
        
        if(overlap_dict != {}):
            #format overlap entries to print:
            ov_labels = [str(ov[0]) for ov in overlap_dict[label][0]]
            ov_areas = [str(ov[1]*pixels_to_nm**2) for ov in overlap_dict[label][0]]
        
        #print to file
        out_file.write(str(label) + ',' + 
                       str(prop['Area']*pixels_to_nm**2) + ',' + 
                       str(prop['Centroid'][1]) + ',' + 
                       str(prop['Centroid'][0]) + ',' + 
                       str(prop['MajorAxisLength']*pixels_to_nm) + ',' + 
                       str(prop['MinorAxisLength']*pixels_to_nm) + ',' + 
                       str(prop['Perimeter']*pixels_to_nm) + ',' +
                       str(prop['Eccentricity']) + ',')
        
        if(overlap_dict != {}):
            out_file.write(';'.join(ov_labels) + ',' +
                           ';'.join(ov_areas) + ',' +
                           str(overlap_dict[label][1]*pixels_to_nm**2) + ',')
            
        out_file.write(str(prop['MaxIntensity']) + ',' + 
                       str(prop['MeanIntensity']) + ',' + 
                       str(prop['MinIntensity']) + '\n')
        
        if(draw_on_figure):
                minr, minc, maxr, maxc = prop['BoundingBox']
                bx = (minc, maxc, maxc, minc, minc)
                by = (minr, minr, maxr, maxr, minr)
                plt.plot(bx, by, '-b', linewidth=bb_line_width)
                plt.text(minc, minr, str(label), color='blue', size=bb_text_size)

#same as above but also writes to an all file 
def save_to_csv2(file_name, all_file, props, overlap_dict, draw_on_figure):
    try: 
        out_file = open(file_name, 'w') 
    except(IOError):
        print 'Error in opening/writing: ' + '"' + file_name + '"' + '.\n'
        return 0
        
    if(overlap_dict == {}):
        out_file.write(',Area,CentroidX,CentroidY,MajorAxisLength,MinorAxisLength,Perimeter,Eccentricity,MaxIntensity,MeanIntensity,MinIntensity\n')
    else:
        out_file.write(',Area,CentroidX,CentroidY,MajorAxisLength,MinorAxisLength,Perimeter,Eccentricity,OverlapSegmentList,OverlapAreaList,TotalOverlapArea,MaxIntensity,MeanIntensity,MinIntensity\n')
    for prop in props: 
        label = prop['Label']
        
        if(overlap_dict != {}):
            #format overlap entries to print:
            ov_labels = [str(ov[0]) for ov in overlap_dict[label][0]]
            ov_areas = [str(ov[1]*pixels_to_nm**2) for ov in overlap_dict[label][0]]
        
        #print to file
        out_file.write(str(label) + ',' + 
                       str(prop['Area']*pixels_to_nm**2) + ',' + 
                       str(prop['Centroid'][1]) + ',' + 
                       str(prop['Centroid'][0]) + ',' + 
                       str(prop['MajorAxisLength']*pixels_to_nm) + ',' + 
                       str(prop['MinorAxisLength']*pixels_to_nm) + ',' + 
                       str(prop['Perimeter']*pixels_to_nm) + ',' +
                       str(prop['Eccentricity']) + ',')
        
        all_file.write(str(label) + ',' + 
                       str(prop['Area']*pixels_to_nm**2) + ',' + 
                       str(prop['Centroid'][1]) + ',' + 
                       str(prop['Centroid'][0]) + ',' + 
                       str(prop['MajorAxisLength']*pixels_to_nm) + ',' + 
                       str(prop['MinorAxisLength']*pixels_to_nm) + ',' + 
                       str(prop['Perimeter']*pixels_to_nm) + ',' +
                       str(prop['Eccentricity']) + ',')
        
        if(overlap_dict != {}):
            out_file.write(';'.join(ov_labels) + ',' +
                           ';'.join(ov_areas) + ',' +
                           str(overlap_dict[label][1]*pixels_to_nm**2) + ',')
            
            all_file.write(';'.join(ov_labels) + ',' +
                           ';'.join(ov_areas) + ',' +
                           str(overlap_dict[label][1]*pixels_to_nm**2) + ',')
            
        out_file.write(str(prop['MaxIntensity']) + ',' + 
                       str(prop['MeanIntensity']) + ',' + 
                       str(prop['MinIntensity']) + '\n')
        
        all_file.write(str(prop['MaxIntensity']) + ',' + 
                       str(prop['MeanIntensity']) + ',' + 
                       str(prop['MinIntensity']) + '\n')
        
        if(draw_on_figure):
                minr, minc, maxr, maxc = prop['BoundingBox']
                bx = (minc, maxc, maxc, minc, minc)
                by = (minr, minr, maxr, maxr, minr)
                plt.plot(bx, by, '-b', linewidth=bb_line_width)
                plt.text(minc, minr, str(label), color='blue', size=bb_text_size)
                
#same as above but only writes to all file, no figure, writes file name, no circularity
def save_to_csv3(all_file, props, overlap_dict, file_name):
    for prop in props: 
        label = prop['Label']
        
        if(overlap_dict != {}):
            #format overlap entries to print:
            ov_labels = [str(ov[0]) for ov in overlap_dict[label][0]]
            ov_areas = [str(ov[1]*pixels_to_nm**2) for ov in overlap_dict[label][0]]
        
        #print to file
        all_file.write(str(label) + ',' + file_name + ',' +
                       str(prop['Area']*pixels_to_nm**2) + ',' + 
                       str(prop['Centroid'][1]) + ',' + 
                       str(prop['Centroid'][0]) + ',' + 
                       str(prop['MajorAxisLength']*pixels_to_nm) + ',' + 
                       str(prop['MinorAxisLength']*pixels_to_nm) + ',' + 
                       str(prop['Perimeter']*pixels_to_nm) + ',' +
                       str(prop['Eccentricity']) + ',')
        
        if(overlap_dict != {}):
            all_file.write(';'.join(ov_labels) + ',' +
                           ';'.join(ov_areas) + ',' +
                           str(overlap_dict[label][1]*pixels_to_nm**2) + ',')
            
        all_file.write(str(prop['MinIntensity']) + ',' + 
                       str(prop['MeanIntensity']) + ',' + 
                       str(prop['MaxIntensity']) + '\n')