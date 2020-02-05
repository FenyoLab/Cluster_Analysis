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

from skimage import io, measure, color, draw
from scipy import ndimage
import numpy as np
import matplotlib.pyplot as plt
import math
import glob
import mahotas as mh
import os
import sys
import readroi as roi
import image_processing as ip
import matplotlib.cm as cm
import datetime as datetime
import re

plugin = 'pil'
reversed_mask = True
pixels_to_nm = 10 # (1 px = 10 nm for reconstructed images)
min_area = 600 # 1000 (simone) # 2400 (esperanza)  
max_distance_to_ID = 600 #nm
show_image = False
id_suffix = '' #'_Reconstruction' #'_ID' #'_IDline'

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
        #fig_size = fig.get_size_inches()
        #w2,h2 = fig_size[0],fig_size[1]
        #fig.set_size_inches([(w2/w)*w,(w2/w)*h])
        #fig.set_dpi((w2/w)*fig.get_dpi())
    a=fig.gca()
    a.set_frame_on(False)
    a.set_xticks([]); a.set_yticks([])
    plt.axis('off')
    plt.xlim(0,h); plt.ylim(w,0)

    #fig.savefig(fileName, dpi=300, transparent=True, bbox_inches='tight', pad_inches=0)
    fig.set_size_inches([w/(90*2),h/(90*2)])
    fig.set_dpi(900)
    fig.savefig(fileName, dpi=900, transparent=True, bbox_inches='tight', pad_inches=0, format='pdf')

############################################################################################################

def cluster_distances(base_dir):
    #reads ROI file containing poly-line that determines the ID line
    #reads ROI file containing line that determines the direction of cell fibers
    #reads MASK of green clusters in cell
    
    #finds the distance from clusters to the ID line along the direction of the cell fibers
    #also calculates measurements of the clusters (and colocalization)

    ###log file
    stamp=datetime.datetime.now().strftime("%Y-%m-%d-%H.%M.%S")
    logfile = open(base_dir + '/log-' + stamp + '.txt', 'w')
    logfile.write('Processing ' + base_dir + '\n')
    ###
    
    file_list = glob.glob(base_dir + '/' + '*.tif')
    if(len(file_list) == 0):
        logfile.write('No tif files found in directory.\n')
        sys.exit(1)
        
    #open CSV file for writing measurements - c2 (GREEN) - includes distances
    c2_out_file = open(base_dir + '/Green_Mask_Measurements.csv', 'w') 
    c2_out_file.write(",File,Label,DistanceID,DistancePosID,ID_Overlap,Area,CentroidX,CentroidY,MajorAxisLength,MinorAxisLength,Perimeter,Eccentricity,Circularity,Perimeter2,Circularity2,OverlapSegmentList,OverlapAreaList,TotalOverlapArea\n")
    
    c2_out_file2 = open(base_dir + '/Green_Mask_Measurements-extra.csv', 'w') 
    c2_out_file2.write(",File,Label,Area,Cx,Cy,Px,Py,Ix,Iy,x1,y1,x2,y2,MinDP\n")
    
    #open CSV file for writing measurements - c1 (RED)
    c1_out_file = open(base_dir + '/Red_Mask_Measurements.csv', 'w') 
    c1_out_file.write(",File,Label,Area,CentroidX,CentroidY,MajorAxisLength,MinorAxisLength,Perimeter,Eccentricity,Circularity,Perimeter2,Circularity2,OverlapSegmentList,OverlapAreaList,TotalOverlapArea\n")
        
    c1_cluster_i = 1
    c2_cluster_i = 1
    
    for image_file in file_list:
        
        if('Green_Mask' not in image_file): continue
        
        logfile.write('Processing ' + image_file + '...\n')
        
        #get corresponding ID line ROI file
        #image file name looks like: 'spool_0.1.tifRGB_Reconstruction_Green_Mask.tif'
        #roi file name looks like: 'spool_0.1.tifRGB_Reconstruction_IDline.zip'
        (head, tail) = os.path.split(image_file)
        (root, ext) = os.path.splitext(tail)
        
        roi_file = root.replace('_Green_Mask', id_suffix)
        roi_file = head + '/' + roi_file + '.zip'
        
        #add title of immediate directory to the root
        (head1,tail1) = os.path.split(head)
        root = tail1 + '_' + root
        
        
        #read the ROI file
        #return is [x1, y1, x2, y2] for a line
        #and array of x,y points for a poly line
        #x value is distance from left, y value is distance from top
        try:
            roi_points = roi.read_roi_zip(roi_file)
        except(IOError, ValueError):
            logfile.write('Error in reading ROI file: ' + roi_file + '\n')
            continue
        
        #old way to read roi file
        #found_line = False
        #found_poly_line = False
        #for points in roi_points:
        #    if((len(points.shape) == 1 and points.shape[0] == 4) or
        #        (len(points.shape) == 2 and points.shape[0] == 2 and points.shape[1] == 2)): #it is a fiber directional line
        #        if(not found_line):
        #            if(len(points.shape) == 1):
        #                x1 = points[0]
        #                y1 = points[1]
        #                x2 = points[2]
        #                y2 = points[3]
        #            else:
        #                x1 = points[0][0]
        #                y1 = points[0][1]
        #                x2 = points[1][0]
        #                y2 = points[1][1]
        #            found_line = True
        #    elif(len(points.shape) == 2 and points.shape[0] > 2 and points.shape[1] == 2):
        #        id_line_coords = points
        #        found_poly_line = True
        #    if(found_line and found_poly_line): break
        
        #new way to read roi file - first roi is the fiber directional line,
        #remaining roi's are the id-lines (there can be more than one, with gaps between)
        read_error = False
        id_line_coords_list = []
        for i,points in enumerate(roi_points):
            if(i == 0):
                if(len(points.shape) == 1 and points.shape[0] == 4):
                    x1 = points[0]
                    y1 = points[1]
                    x2 = points[2]
                    y2 = points[3]
                elif(len(points.shape) == 2 and points.shape[0] == 2 and points.shape[1] == 2):
                    x1 = points[0][0]
                    y1 = points[0][1]
                    x2 = points[1][0]
                    y2 = points[1][1]
                else:
                    read_error = True
                    logfile.write('Error in reading fiber directional line: ' + roi_file + '\n')
                    break
            else:
                if(len(points.shape) == 1 and points.shape[0] == 4): #shape is 4x1, change to 2x2 shape
                    points_f = []
                    points_f.append(points[:2])
                    points_f.append(points[2:4])
                elif(len(points.shape) == 2 and points.shape[0] >= 2 and points.shape[1] == 2):
                    points_f = points
                else:
                    read_error = True
                    logfile.write('Error in reading ID line: ' + roi_file + '\n')
                    break
                id_line_coords_list.append(points_f)
                
        if(read_error): continue
        
        #read in green cluster mask and red cluster mask
        image_file2 = tail.replace('Green_Mask', 'Red_Mask')
        image_file2 = head + '/' + image_file2
        c2_image = io.imread(image_file, True, plugin) 
        c1_image = io.imread(image_file2, True, plugin) 
        
        if(plugin == 'freeimage'): 
            c2_mask = c2_image > 0
        else:
            c2_mask = c2_image < 255
        s = [[1,1,1],[1,1,1],[1,1,1]]
        
        #filter out clusters with area less than min area
        c2_labeled_clusters, c2_num_labels = ndimage.label(c2_mask, structure=s)
        c2_props = measure.regionprops(c2_labeled_clusters, ['Area', 'Coordinates'])
        c2_mask = filter_mask_for_area(c2_mask, c2_props)
        
        c2_labeled_clusters, c2_num_labels = ndimage.label(c2_mask, structure=s)
        c2_props = measure.regionprops(c2_labeled_clusters, [
                    'Area',
                    'Centroid',
                    'Eccentricity',
                    'MajorAxisLength',
                    'MinorAxisLength',
                    'MaxIntensity',
                    'MeanIntensity',
                    'MinIntensity',
                    'Perimeter',
                    'BoundingBox',
                    'Coordinates',
                    'Image'],c2_image)
        
        if(plugin == 'freeimage'): 
            c1_mask = c1_image > 0
        else:
            c1_mask = c1_image < 255
        s = [[1,1,1],[1,1,1],[1,1,1]]
        
        #filter out clusters with area less than min area
        c1_labeled_clusters, c1_num_labels = ndimage.label(c1_mask, structure=s)
        c1_props = measure.regionprops(c1_labeled_clusters, ['Area', 'Coordinates'])
        c1_mask = filter_mask_for_area(c1_mask, c1_props)
        
        c1_labeled_clusters, c1_num_labels = ndimage.label(c1_mask, structure=s)
        c1_props = measure.regionprops(c1_labeled_clusters, [
                    'Area',
                    'Centroid',
                    'Eccentricity',
                    'MajorAxisLength',
                    'MinorAxisLength',
                    'MaxIntensity',
                    'MeanIntensity',
                    'MinIntensity',
                    'Perimeter',
                    'BoundingBox',
                    'Coordinates',
                    'Image'],c1_image)
        
        #extend the poly line outward at the ends so that any c2 clusters near the end will not be missed
        #don't do this - esperanza will try to draw the line properly
        #x = id_line_coords[0][0] - (id_line_coords[1][0]-id_line_coords[0][0])
        #y = id_line_coords[0][1] - (id_line_coords[1][1]-id_line_coords[0][1])
        #id_line_coords = np.insert(id_line_coords,0,[x,y],axis=0)
        #l = len(id_line_coords) - 1
        #x = id_line_coords[l][0] - (id_line_coords[l-1][0]-id_line_coords[l][0])
        #y = id_line_coords[l][1] - (id_line_coords[l-1][1]-id_line_coords[l][1])
        #id_line_coords = np.append(id_line_coords,[[x,y]],axis=0)
        
        #draw lines on display image
        if(show_image):
            display_image = color.gray2rgb(c2_image)
            #draw fb line on image
            rr,cc = draw.line(int(y1),int(x1),int(y2),int(x2))
            display_image[rr, cc] = [255,0,0]
            #draw id line on image
            for id_line_coords in id_line_coords_list:
                for i, coords in enumerate(id_line_coords):
                    if(i < (len(id_line_coords)-1)): #examining 2 at a time
                        rr,cc = draw.line(int(coords[1]),int(coords[0]),int(id_line_coords[i+1][1]),int(id_line_coords[i+1][0]))
                        #drawing start and end lines in a different color - also checking if they go out of bounds of image
                        #not used since poly line extension feature not used (see above)
                        #if(i == 0 or i == (len(id_line_coords)-2)):
                        #    mask = np.ones(len(rr), dtype=bool)
                        #    for i, r in enumerate(rr):
                        #        if(r < 0 or r >= len(display_image)):
                        #            mask[i] = False
                        #    for i, c in enumerate(cc):
                        #        if(c < 0 or c >= len(display_image[0])):
                        #            mask[i] = False
                        #    rr = rr[mask]
                        #    cc = cc[mask]        
                        #    display_image[rr, cc] = [255,127,0]
                        #else: display_image[rr, cc] = [255,0,0]
                        display_image[rr, cc] = [255,0,0]
            #draw orange line for the gaps in the ID line - seeing if this will work...
            #for i,id_line_coords in enumerate(id_line_coords_list):
            #    if(i < (len(id_line_coords_list)-1)): #examining 2 at a time
            #        j = len(id_line_coords)-1
            #        rr,cc = draw.line(int(id_line_coords[j][1]),int(id_line_coords[j][0]),int(id_line_coords_list[i+1][0][1]),int(id_line_coords_list[i+1][0][0]))
            #        display_image[rr, cc] = [255,127,0]
                    
        #calculate overlap
        #calculate overlap area for each c1 cluster
        c1_overlap_dict = {}
        ip.calculate_overlap(c1_props, c2_labeled_clusters, c1_overlap_dict)
                
        #calculate overlap area for each c2 cluster
        c2_overlap_dict = {}
        ip.calculate_overlap(c2_props, c1_labeled_clusters, c2_overlap_dict)
        
        #slope, of fiber directional line
        if(x2 == x1): fm = None 
        else: fm = (float(y2 - y1))/(x2-x1)
    
        #c2 - measure clusters, find distance to ID line and output data to the CSV file
        labels_array = []
        for prop in c2_props:
            
            #############DISTANCE FROM ID line CALCULATION#########################################################
            #measure distance from centroid of cluster, along fiber directional line, to intersecting ID line:
            
            #cx,cy is coordinates of centroid of cluster
            cy,cx = prop['Centroid']
            label = prop['Label']
            image = prop['Image'] #sliced binary region image same size as bounding box
            r1,c1,r2,c2 = prop['BoundingBox']
            cluster_coords = prop['Coordinates']
            perim = prop['Perimeter'] * pixels_to_nm
            area = prop['Area'] * pixels_to_nm * pixels_to_nm
            
            #we already filter by area, above
            #if(area <= min_area):
            #    continue
            
            circ = 4 * math.pi * ((area - .5*perim) / perim**2)
            if(circ > 1.): circ = 1.
            
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
            
            #color green the c2 cluster
            
            #get perimeter pixels
            perim_coords = []
            for i,val1 in enumerate(image_perim):
                for j,val2 in enumerate(val1):
                    if(val2):
                        perim_coords.append([j + c1,i + r1]) #put back in whole image
                        #outline perimeters of clusters we are measuring on the image
                        if(show_image):
                            display_image[i + r1][j + c1] = [255,0,0]
            
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
            
            touch_line = False
            cross_line_count = 0
            found_min_d = False
            for h, id_line_coords in enumerate(id_line_coords_list):
                for i, coords in enumerate(id_line_coords):
                    if(i < (len(id_line_coords)-1)): #examining 2 at a time
                        
                        id_x1 = float(coords[0])
                        id_y1 = float(coords[1])
                        id_x2 = float(id_line_coords[i+1][0])
                        id_y2 = float(id_line_coords[i+1][1])
                        
                        #slope, of current ID line
                        if(id_x2 == id_x1): id_m = None
                        else: id_m = (float(id_y2-id_y1))/(id_x2-id_x1)
                        
                        if(id_m == fm): continue
                        
                        #check if poly line is touching cluster
                        for cl_y,cl_x in cluster_coords:
                            if(id_m != None): v = abs(cl_y - id_m*cl_x - (-1*id_m*id_x1+id_y1)) / (1+id_m*id_m)**.5 #distance between cluster point and line
                            if((cl_x <= max(id_x2,id_x1) and cl_x >= min(id_x2,id_x1)) and 
                                (cl_y <= max(id_y2,id_y1) and cl_y >= min(id_y2,id_y1))): # cl_x,cl_y between st/end of ID line
                                if((id_m == None and abs(cl_x-id_x1) < 1.) or (v < 1.)): 
                                    #set diff to 0
                                    min_d = 0
                                    found_min_d = True
                                    inter_i = [h,i]
                                    edge_x = 0
                                    edge_y = 0
                                    min_x_inter = 0
                                    min_y_inter = 0
                                    touch_line = True
                                    break
                            
                        #NOTE: this may not work b/c of gaps in id line ---
                        #checking for which side of the poly line the cluster is on
                        #count how many times line from centroid to a point on the 'green' side crosses poly line
                        #even - it is on 'green' side, odd - not on 'green' side
                        if(x2 != cx or id_m != None): #otherwise lines are both vertical and do not cross
                            if(x2 == cx):
                                x_inter = x2
                                y_inter = id_m*(x_inter-id_x1)+id_y1
                            else:
                                cl_m = (float(y2-cy))/(x2-cx)
                                if(id_m == None):
                                    x_inter = id_x1
                                else:
                                    x_inter = (id_y1 - cy + cl_m*cx - id_m*id_x1)/(cl_m-id_m)
                                y_inter = cl_m*(x_inter-cx) + cy
                            if((x_inter <= max(id_x2,id_x1) and x_inter >= min(id_x2,id_x1)) and 
                                (y_inter <= max(id_y2,id_y1) and y_inter >= min(id_y2,id_y1))  and
                                (x_inter < max(cx,x2) and x_inter > min(cx,x2)) and
                                (y_inter < max(cy,y2) and y_inter > min(cy,y2))):
                                #check if it intersection is between the st/end of ID line and cx,cy and x2,y2
                                #NOTE: what if it touches line but does not cross? *CHECK*
                                if(x_inter == id_x2 and y_inter == id_y2 and i < (len(id_line_coords)-2) ):
                                    #if it crosses at x2,y2 of id line section, don't count it b/c it will be counted
                                    #in the next iteration, since the point is part of both line sections
                                    #unless it is the last id line section, then count it since it is only part of one line
                                    pass 
                                else:
                                    cross_line_count += 1
                        
                        if(not touch_line):    
                            for px,py in perim_coords:
                                #find intersection between perimeter point and ID line, along line with slope fm (fiber directional line)
                                if(fm == None):
                                    x_inter = px
                                    y_inter = id_m*(x_inter-id_x1) + id_y1
                                elif(id_m == None):
                                    x_inter = id_x1
                                    y_inter = fm*(x_inter-px) + py
                                else:
                                    x_inter = (id_y1 - py + fm*px - id_m*id_x1)/(fm-id_m) # <- consider id_m and fm == None!
                                    y_inter = fm*(x_inter-px) + py
                                
                                #check if it intersection is between the st/end of ID line
                                if((x_inter <= max(id_x2,id_x1) and x_inter >= min(id_x2,id_x1)) and
                                    (y_inter <= max(id_y2,id_y1) and y_inter >= min(id_y2,id_y1))):
                                    #record distance, if it is current minimum
                                    d = ((px-x_inter)**2 + (py-y_inter)**2)**.5
                                    if(not found_min_d or d < min_d):
                                        min_d = d
                                        found_min_d = True
                                        inter_i = [h,i]
                                        inter_i_m = id_m
                                        edge_x = px
                                        edge_y = py
                                        min_x_inter = x_inter
                                        min_y_inter = y_inter
                        else: break
            
            #########################################################################################################
            
            #get rest of cluster measurements
            
            #write to file
            if(cross_line_count%2 != 0):
                if(found_min_d):
                    min_d = -1*min_d
            if(found_min_d):
                id_line_coords = id_line_coords_list[inter_i[0]]
                c2_out_file2.write(str(c2_cluster_i) + ',' + root + ',' + str(label) + ',' + str(area) + ',' + str(cx * pixels_to_nm) + ',' + str(cy * pixels_to_nm) + ',' +
                               str(edge_x * pixels_to_nm) + ',' + str(edge_y * pixels_to_nm) + ',' +
                               str(min_x_inter * pixels_to_nm) + ',' + str(min_y_inter * pixels_to_nm) + ',' +
                               str(id_line_coords[inter_i[1]][0] * pixels_to_nm) + ',' + str(id_line_coords[inter_i[1]][1] * pixels_to_nm) + ',' +
                               str(id_line_coords[inter_i[1]+1][0] * pixels_to_nm) + ',' + str(id_line_coords[inter_i[1]+1][1] * pixels_to_nm) + ',' +
                               str(min_d * pixels_to_nm) + '\n')
                #draw d line on image
                if(not touch_line and show_image):
                    rr,cc = draw.line(int(edge_y),int(edge_x),int(min_y_inter),int(min_x_inter))
                    if(cross_line_count%2 != 0):
                        display_image[rr, cc] = [0,0,255]
                    else:
                        display_image[rr, cc] = [0,255,0]
            else: #the end of cell line is unknown near this cluster
                c2_out_file2.write(str(c2_cluster_i) + ',' + root + ',' + str(label) + ',' + str(cx * pixels_to_nm) + ',' + str(cy * pixels_to_nm) + '\n')
                logfile.write('Could not find cluster distance to cell line. File = ' + root + ', Label = ' + str(label) + '\n')
                
            if(show_image):
                #save label and x,y of the cluster for labeling the display image
                labels_array.append([str(label), c2, r2])
                  
            #format overlap entries to print:
            ov_labels = [str(ov[0]) for ov in c2_overlap_dict[label][0]]
            ov_areas = [str(ov[1]*pixels_to_nm**2) for ov in c2_overlap_dict[label][0]]
            
            #output blank for distance if there was an error in finding it
            if(not found_min_d):
                min_d = ''
                min_d_pos = ''
                id_overlap = ''
            else:
                if(min_d > 0): min_d_pos = min_d
                else: min_d_pos = 0
                if(c2_overlap_dict[label][1] > 0 or min_d <= 0): id_overlap = '1'
                else: id_overlap = '0'
            if(min_d == '' or (abs(min_d*pixels_to_nm) <= max_distance_to_ID)): 
            #only write to file if it is <= 600nm from the ID line (max_distance_to_ID)
                c2_out_file.write(str(c2_cluster_i) + ',' + root + ',' + str(label) + ',' + str(min_d*pixels_to_nm) + ',' + str(min_d_pos*pixels_to_nm) +
                                  ',' + id_overlap + ',' + str(area) + ',' + str(cx * pixels_to_nm) + ',' +
                                str(cy * pixels_to_nm) + ',' + str(prop['MajorAxisLength']*pixels_to_nm) + ',' + str(prop['MinorAxisLength']*pixels_to_nm) + ',' +
                                str(perim) + ',' + str(prop['Eccentricity']) + ',' + str(circ) + ',' +
                                str(manual_perim) + ',' + str(manual_circ) + ',')
                
                c2_out_file.write(';'.join(ov_labels) + ',' + ';'.join(ov_areas) + ',' + str(c2_overlap_dict[label][1]*pixels_to_nm**2) + '\n')
            
            c2_cluster_i += 1
            
        #c1 - measure clusters, and output data to the CSV file
        for prop in c1_props:
            #cx,cy is coordinates of centroid of cluster
            cy,cx = prop['Centroid']
            label = prop['Label']
            area = prop['Area'] * pixels_to_nm * pixels_to_nm
            perim = prop['Perimeter'] * pixels_to_nm
            image = prop['Image'] #sliced binary region image same size as bounding box
            if(area <= min_area):
                continue
            circ = 4 * math.pi * ((area - .5*perim) / perim**2)
            if(circ > 1.): circ = 1.
            
            #find perimeter manually, trying to match ImageJ perimater (and therefore circulatiry)
            
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
            
            #write to file
            c1_out_file.write(str(c1_cluster_i) + ',' + root + ',' + str(label) + ',' + str(area) + ',' + str(cx * pixels_to_nm) + ',' +
                            str(cy * pixels_to_nm) + ',' + str(prop['MajorAxisLength']*pixels_to_nm) + ',' + str(prop['MinorAxisLength']*pixels_to_nm) + ',' +
                            str(perim) + ',' + str(prop['Eccentricity']) + ',' + str(circ) + ',' +
                            str(manual_perim) + ',' + str(manual_circ) + ',')
            #format overlap entries to print:
            ov_labels = [str(ov[0]) for ov in c1_overlap_dict[label][0]]
            ov_areas = [str(ov[1]*pixels_to_nm**2) for ov in c1_overlap_dict[label][0]]
            c1_out_file.write(';'.join(ov_labels) + ',' + ';'.join(ov_areas) + ',' + str(c1_overlap_dict[label][1]*pixels_to_nm**2) + '\n')
            
            c1_cluster_i += 1
            
        if(show_image):
            show_image_file = root.replace('Green_Mask', 'Green_Mask_Measurements')
            #show_image_file1 = head + '/' + show_image_file + '.png'
            #io.imsave(show_image_file1, display_image)
            
            plt.imshow(display_image)
            for l in labels_array:
                plt.text(l[1],l[2],l[0],fontsize='1')
            show_image_file1 = head + '/' + show_image_file + '.pdf'
            SaveFigureAsImage(show_image_file1, plt.gcf(), orig_size=(int(display_image.shape[0]), int(display_image.shape[1])))
            plt.clf()
            plt.close()
            
def find_ImageJ_perimeter(image):
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
    return (manual_perim, image_perim)

def filter_mask_for_area(mask, props):
    for prop in props:
        area = prop['Area'] * pixels_to_nm * pixels_to_nm
        if(area <= min_area):
            #remove this cluster from mask since it's smaller than min
            coords = prop['Coordinates']
            for c in coords:
                mask[c[0],c[1]] = False
    return mask

def closest_cluster(perim_coords, bbox, image_mask, image_mask2=[]): 
    #bbox is bounding box around the cluster of interest
    #increase bbox all around, then check if we have found any additional clusters
    #if so, expand to look in a circle around the cluster, and find which surronding cluster is closest
    
    if(len(image_mask2) == 0):
        #we are comparing one image to itself (e.g. red clusters to red clusters)
        label_min = 1
        image_mask2 = image_mask
    else: label_min = 0 #2 images, e.g. red clusters to green clusters
    
    image_width = len(image_mask[0])
    image_height = len(image_mask)
    center_bbox = (bbox[1]+((bbox[3]-bbox[1])/2),bbox[0]+((bbox[2]-bbox[0])/2))
    box_incr = max((bbox[3] - bbox[1]),(bbox[2] - bbox[0])) #increment by max(h,w) of box, to make sure we include the closest cluster (since clusters are irregular shapes)
    x_minus = bbox[1]-box_incr
    x_plus = bbox[3]+box_incr
    y_minus = bbox[0]-box_incr
    y_plus = bbox[2]+box_incr
    ret_val = expand_slice(x_minus, y_minus, x_plus, y_plus, image_width, image_height, image_mask2)
    while(ret_val):
        (labeled_clusters, num_labels) = ret_val
        if(num_labels > label_min):
            #now that we have found surronding cluster(s),
            #increment slice one more time and then
            #expand to look in a circle around cluster, not a rectangle
            
            #increment slice 
            x_minus -= box_incr
            y_minus -= box_incr
            x_plus += box_incr
            y_plus += box_incr
            
            if(x_minus < 0):
                x_minus = 0 
            if(y_minus < 0):
                y_minus = 0
            if(x_plus >= image_width):
                x_plus = image_width-1
            if(y_plus >= image_height):
                y_plus = image_height-1
            
            #cur_h_width = math.ceil(len(labeled_clusters[0])/2.)
            #cur_h_height = math.ceil(len(labeled_clusters)/2.)
            #cur_h_width = math.ceil((x_plus-x_minus)/2.)
            #cur_h_height = math.ceil((y_plus-y_minus)/2.)
            
            #expand to circle radius
            cur_h_width = max((x_plus-center_bbox[0]),(center_bbox[0]-x_minus))
            cur_h_height = max((y_plus-center_bbox[1]),(center_bbox[1]-y_minus))
            diag = math.ceil(((cur_h_width)**2 + (cur_h_height)**2)**.5)
            x_minus -= (diag-cur_h_width)
            y_minus -= (diag-cur_h_height)
            x_plus += (diag-cur_h_width)
            y_plus += (diag-cur_h_height)
            
            ret_val = expand_slice(x_minus, y_minus, x_plus, y_plus, image_width, image_height, image_mask2)
            if(ret_val): (labeled_clusters, num_labels) = ret_val
            else:
                print "Alert! In unfinished else!"
                pass #need to reset x,y plus/minus to entire image (but not over), then rerun expand_slice
            
            #found other surrounding clusters, find which is nearest
            if(x_minus < 0): x_minus = 0 
            if(y_minus < 0): y_minus = 0
            if(x_plus >= image_width): x_plus = image_width-1
            if(y_plus >= image_height): y_plus = image_height-1
            
            #ignore cluster we are measuring from (one image case)
            if(label_min == 1): cluster1_label = labeled_clusters[perim_coords[0][1] - bbox[0]+(bbox[0]-y_minus)][perim_coords[0][0] - bbox[1]+(bbox[1]-x_minus)]
            else: cluster1_label = -1
            
            props = measure.regionprops(labeled_clusters, ['Coordinates', 'Image'])
            min_distance = -1
            closest_cluster_edge_coord = -1
            cluster_perim_coord = -1
            for prop in props:
                label = prop['Label']
                if(label != cluster1_label):
                    #find closest distance for this 'cluster'
                    cur_coords = prop['Coordinates']
                    for cur_c in cur_coords:
                        image_cur_c = (cur_c[1] + bbox[1] - (bbox[1]-x_minus),cur_c[0] + bbox[0] - (bbox[0]-y_minus))
                        for perim_c in perim_coords:
                            distance = find_distance(image_cur_c[0],image_cur_c[1],perim_c[0],perim_c[1])
                            if(distance < min_distance or min_distance < 0):
                                min_distance = distance
                                closest_cluster_edge_coord = image_cur_c
                                cluster_perim_coord = perim_c
            return (min_distance, closest_cluster_edge_coord, cluster_perim_coord)
        else:
            #guard against infinite loop, although it shouldn't happen
            if(x_minus < 0 and y_minus < 0 and x_plus >= image_width and y_plus >= image_height): return False
            
        x_minus -= box_incr
        y_minus -= box_incr
        x_plus += box_incr
        y_plus += box_incr
        ret_val = expand_slice(x_minus, y_minus, x_plus, y_plus, image_width, image_height, image_mask2)
    
def expand_slice(x_minus, y_minus, x_plus, y_plus, image_width, image_height, image_mask):
    s = [[1,1,1],[1,1,1],[1,1,1]]
    
    if(x_minus < 0):
        x_minus = 0 
    if(y_minus < 0):
        y_minus = 0
    if(x_plus >= image_width):
        x_plus = image_width-1
    if(y_plus >= image_height):
        y_plus = image_height-1
    bbox_slice = image_mask[y_minus:y_plus,x_minus:x_plus]
    labeled_clusters, num_labels = ndimage.label(bbox_slice, structure=s)
    return (labeled_clusters, num_labels)
        
def find_distance(x1, y1, x2, y2):
    return ((x1-x2)**2 + (y1-y2)**2)**.5

def cluster_distances2(base_dir):
    #reads MASK of green clusters in cell/reads MASK of red clusters in cell
    #finds distance from each green cluster to the closest red cluster, and vice versa
    #also finds distance from each green cluster to the closest green cluster (and for red)
    #distance is measured in 2 ways - from centroid and from closest edge
    #also calculates measurements of the clusters (and colocalization)
    
    #we don't have an ID line, as in cluster_distances(), which is why we are proceeding as described above

    ###log file
    stamp=datetime.datetime.now().strftime("%Y-%m-%d-%H.%M.%S")
    logfile = open(base_dir + '/log-' + stamp + '.txt', 'w')
    logfile.write('Processing ' + base_dir + '\n')
    ###
    
    file_list = glob.glob(base_dir + '/' + '*.tif')
    if(len(file_list) == 0):
        logfile.write('No tif files found in directory.\n')
        sys.exit(1)
        
    #open CSV file for writing measurements - c2 (GREEN) - includes distances
    c2_out_file = open(base_dir + '/Green_Mask_Measurements.csv', 'w') 
    c2_out_file.write(",File,Label,EdgeDistanceClosestGreen,EdgeLabelClosestGreen,EdgeAreaClosestGreen,EdgeDistanceClosestRed,EdgeLabelClosestRed,EdgeAreaClosestRed,CentroidDistanceClosestGreen,CentroidLabelClosestGreen,CentroidAreaClosestGreen,CentroidDistanceClosestRed,CentroidLabelClosestRed,CentroidAreaClosestRed,Area,CentroidX,CentroidY,MajorAxisLength,MinorAxisLength,Perimeter,Eccentricity,Circularity,OverlapSegmentList,OverlapAreaList,TotalOverlapArea,AreaOverlapClusters\n")
    
    #open CSV file for writing measurements - c1 (RED)
    c1_out_file = open(base_dir + '/Red_Mask_Measurements.csv', 'w') 
    c1_out_file.write(",File,Label,EdgeDistanceClosestRed,EdgeLabelClosestRed,EdgeAreaClosestRed,EdgeDistanceClosestGreen,EdgeLabelClosestGreen,EdgeAreaClosestGreen,CentroidDistanceClosestRed,CentroidLabelClosestRed,CentroidAreaClosestRed,CentroidDistanceClosestGreen,CentroidLabelClosestGreen,CentroidAreaClosestGreen,Area,CentroidX,CentroidY,MajorAxisLength,MinorAxisLength,Perimeter,Eccentricity,Circularity,OverlapSegmentList,OverlapAreaList,TotalOverlapArea,AreaOverlapClusters\n")
        
    c1_cluster_i = 1
    c2_cluster_i = 1
    for image_file in file_list:
        
        #for esperanza
        if('Green_Mask' not in image_file): continue
        (head, tail) = os.path.split(image_file)
        red_image_file = tail.replace('Green_Mask', 'Red_Mask')
        red_image_file = glob.glob(head + '/' + red_image_file)
        if(len(red_image_file) == 1):
            red_image_file = red_image_file[0]
        else:
            logfile.write('Error in finding red image file corresponding to: "' + tail + '".\n')
            continue
        
        #for simone
        #if('(green)' not in image_file): continue
        #(head, tail) = os.path.split(image_file)
        #red_image_file = base_dir + '/' + tail.replace('(green)', '(red)')
        #red_image_file = glob.glob(red_image_file)
        #if(len(red_image_file) == 1):
        #    red_image_file = red_image_file[0]
        #else:
        #    logfile.write('Error in finding red image file in glob: "' + red_image_file + '".\n')
        #    continue
        
        #for the simulations
        #if('c2_mask' not in image_file): continue
        #(head, tail) = os.path.split(image_file)
        #red_image_file = tail.replace('c2_mask', 'c1_mask')
        #red_image_file = glob.glob(head + '/' + red_image_file)
        #if(len(red_image_file) == 1):
        #    red_image_file = red_image_file[0]
        #else:
        #    logfile.write('Error in finding red image file corresponding to: "' + tail + '".\n')
        #    continue
        
        logfile.write('Processing ' + image_file + '...\n')
        
        #read in green cluster mask and red cluster mask
        c2_image = io.imread(image_file, True, plugin) 
        c1_image = io.imread(red_image_file, True, plugin) 
        
        if(reversed_mask): 
            c2_mask = c2_image < 255
        else:
            c2_mask = c2_image > 0
            
        #if(plugin == 'freeimage'): 
        #    c2_mask = c2_image > 0
        #else:
        #    c2_mask = c2_image < 255
        s = [[1,1,1],[1,1,1],[1,1,1]]
        
        #filter out clusters with area less than min area
        c2_labeled_clusters, c2_num_labels = ndimage.label(c2_mask, structure=s)
        c2_props = measure.regionprops(c2_labeled_clusters, ['Area', 'Coordinates'])
        c2_mask = filter_mask_for_area(c2_mask, c2_props)
        
        c2_labeled_clusters, c2_num_labels = ndimage.label(c2_mask, structure=s)
        c2_props = measure.regionprops(c2_labeled_clusters, [
                    'Area',
                    'Centroid',
                    'Eccentricity',
                    'MajorAxisLength',
                    'MinorAxisLength',
                    'Perimeter',
                    'BoundingBox',
                    'Coordinates',
                    'Image'])
        
        #create area dict
        c2_area_dict = {}
        for prop in c2_props:
            c2_area_dict[prop['Label']] = prop['Area']
        
        if(reversed_mask): 
            c1_mask = c1_image < 255
        else:
            c1_mask = c1_image > 0
            
        #if(plugin == 'freeimage'): 
        #    c1_mask = c1_image > 0
        #else:
        #    c1_mask = c1_image < 255
        s = [[1,1,1],[1,1,1],[1,1,1]]
        
        #filter out clusters with area less than min area
        c1_labeled_clusters, c1_num_labels = ndimage.label(c1_mask, structure=s)
        c1_props = measure.regionprops(c1_labeled_clusters, ['Area', 'Coordinates'])
        c1_mask = filter_mask_for_area(c1_mask, c1_props)
        
        c1_labeled_clusters, c1_num_labels = ndimage.label(c1_mask, structure=s)
        c1_props = measure.regionprops(c1_labeled_clusters, [
                    'Area',
                    'Centroid',
                    'Eccentricity',
                    'MajorAxisLength',
                    'MinorAxisLength',
                    'Perimeter',
                    'BoundingBox',
                    'Coordinates',
                    'Image'])
        
        #create area dict
        c1_area_dict = {}
        for prop in c1_props:
            c1_area_dict[prop['Label']] = prop['Area']
        
        #calculate overlap
        #calculate overlap area for each c1 cluster
        c1_overlap_dict = {}
        ip.calculate_overlap(c1_props, c2_labeled_clusters, c1_overlap_dict)
                
        #calculate overlap area for each c2 cluster
        c2_overlap_dict = {}
        ip.calculate_overlap(c2_props, c1_labeled_clusters, c2_overlap_dict)
        
        if(show_image):
            c1_display_image = color.gray2rgb(c1_image)
            c1_display_image_cent = color.gray2rgb(c1_image)
            c1_labels_array = []
            c2_display_image = color.gray2rgb(c2_image)
            c2_display_image_cent = color.gray2rgb(c2_image)
            c2_labels_array = []
        
        (head, tail) = os.path.split(red_image_file)
        (root, ext) = os.path.splitext(tail)
        for prop in c1_props:
            #cx,cy is coordinates of centroid of cluster
            cy,cx = prop['Centroid']
            label = prop['Label']
            area = prop['Area'] * pixels_to_nm * pixels_to_nm
            image = prop['Image'] #sliced binary region image same size as bounding box
            bbox = prop['BoundingBox']
            
            #calculate ImageJ perimeter/circularity:
            (manual_perim, image_perim) = find_ImageJ_perimeter(image)
            manual_circ = 4 * math.pi * ((area - .5*manual_perim) / manual_perim**2)
            if(manual_circ > 1.): manual_circ = 1.
            
            #get perimeter pixels
            perim_coords = []
            for i,val1 in enumerate(image_perim):
                for j,val2 in enumerate(val1):
                    if(val2):
                        #perim_coords.append([j,i])
                        perim_coords.append([j + bbox[1],i + bbox[0]]) #put back in whole image
                        if(show_image):
                            c1_display_image[i + bbox[0]][j + bbox[1]] = [255,0,0]
                            c2_display_image[i + bbox[0]][j + bbox[1]] = [255,0,0]
                            c1_display_image_cent[i + bbox[0]][j + bbox[1]] = [255,0,0]
                            c2_display_image_cent[i + bbox[0]][j + bbox[1]] = [255,0,0]
                            
            #find closest cluster c1 to c1
            if(c1_num_labels > 1): #only look if there is > 1 cluster
                
                #centroid to centroid
                min_cent_distance = -1
                for prop_ in c1_props:
                    cy_,cx_ = prop_['Centroid']
                    label_ = prop_['Label']
                    if(label_ != label):
                        cent_distance = find_distance(cx,cy,cx_,cy_)
                        if(cent_distance < min_cent_distance or min_cent_distance == -1):
                            min_cent_distance = cent_distance
                            closest_cluster_label_cent = label_
                            closest_cluster_pt_cent = (cx_,cy_)
                            closest_cluster_area_cent = c1_area_dict[closest_cluster_label_cent]
                            
                #edge to edge
                (min_distance, closest_cluster_pt, closest_pt) = closest_cluster(perim_coords, bbox, c1_mask)
                closest_cluster_label = c1_labeled_clusters[closest_cluster_pt[1]][closest_cluster_pt[0]]
                closest_cluster_area = c1_area_dict[closest_cluster_label]
                            
                min_distance *= pixels_to_nm
                min_cent_distance *= pixels_to_nm
                closest_cluster_area *= (pixels_to_nm**2)
                closest_cluster_area_cent *= (pixels_to_nm**2)
            else:
                min_distance = ''
                closest_cluster_label = ''
                min_cent_distance = ''
                closest_cluster_label_cent = ''
                closest_cluster_area = ''
                closest_cluster_area_cent = ''
            
            #find closest cluster c1 to c2
            if(c2_num_labels > 0): #must have atleast one c2 found
                if(c1_overlap_dict[label][1] == 0):
                    #edge to edge
                    (min_distance2, closest_cluster_pt2, closest_pt2) = closest_cluster(perim_coords, bbox, c1_mask, c2_mask)
                    closest_cluster_label2 = c2_labeled_clusters[closest_cluster_pt2[1]][closest_cluster_pt2[0]]
                    closest_cluster_area2 = c2_area_dict[closest_cluster_label2]
                    
                    #centroid to centroid
                    min_cent_distance2 = -1
                    for prop_ in c2_props:
                        cy_,cx_ = prop_['Centroid']
                        label_ = prop_['Label']
                        cent_distance = find_distance(cx,cy,cx_,cy_)
                        if(cent_distance < min_cent_distance2 or min_cent_distance2 == -1):
                            min_cent_distance2 = cent_distance
                            closest_cluster_label_cent2 = label_
                            closest_cluster_pt_cent2 = (cx_,cy_)
                            closest_cluster_area_cent2 = c2_area_dict[closest_cluster_label_cent2]
                else:
                    min_distance2 = 0
                    closest_cluster_label2 = ''
                    min_cent_distance2 = 0
                    closest_cluster_label_cent2 = ''
                    closest_cluster_area2 = 0
                    closest_cluster_area_cent2 = 0
                    
                min_distance2 *= pixels_to_nm
                min_cent_distance2 *= pixels_to_nm
                closest_cluster_area2 *= (pixels_to_nm**2)
                closest_cluster_area_cent2 *= (pixels_to_nm**2)
            else:
                min_distance2 = ''
                closest_cluster_label2 = ''
                min_cent_distance2 = ''
                closest_cluster_label_cent2 = ''
                closest_cluster_area2 = ''
                closest_cluster_area_cent2 = ''
            
            #write to file
            c1_out_file.write(str(c1_cluster_i) + ',' + root + ',' + str(label) + ',' +
                            str(min_distance) + ',' + str(closest_cluster_label) + ',' + str(closest_cluster_area) + ',' + 
                            str(min_distance2) + ',' + str(closest_cluster_label2) + ',' + str(closest_cluster_area2) + ',' + 
                            str(min_cent_distance) + ',' + str(closest_cluster_label_cent) + ',' + str(closest_cluster_area_cent) + ',' + 
                            str(min_cent_distance2) + ',' + str(closest_cluster_label_cent2) + ',' + str(closest_cluster_area_cent2) + ',' + 
                            str(area) + ',' + str(cx * pixels_to_nm) + ',' + str(cy * pixels_to_nm) + ',' + str(prop['MajorAxisLength']*pixels_to_nm) +
                            ',' + str(prop['MinorAxisLength']*pixels_to_nm) + ',' + str(manual_perim) + ',' + str(prop['Eccentricity']) + ',' + str(manual_circ) + ',')
            #format overlap entries to print:
            ov_labels = [str(ov[0]) for ov in c1_overlap_dict[label][0]]
            ov_areas = [str(ov[1]*pixels_to_nm**2) for ov in c1_overlap_dict[label][0]]
            
            #areas of overlapping clusters (sum if more than one)
            area_ov_clusters = 0
            for ov in c1_overlap_dict[label][0]:
                area_ov_clusters += c2_area_dict[ov[0]]
            area_ov_clusters *= (pixels_to_nm**2)
            
            c1_out_file.write(';'.join(ov_labels) + ',' + ';'.join(ov_areas) + ',' + str(c1_overlap_dict[label][1]*pixels_to_nm**2) + ',' + str(area_ov_clusters) + '\n')
            
            if(show_image):
                #save label and x,y of the cluster for labeling the display image
                c1_labels_array.append([str(label), bbox[3], bbox[2]])
                
                if(c1_num_labels > 1):
                    #edge distance
                    rr,cc = draw.line(int(closest_cluster_pt[1]),int(closest_cluster_pt[0]),int(closest_pt[1]),int(closest_pt[0]))
                    c1_display_image[rr, cc] = [0,0,255]
                
                    #centroid distance
                    rr,cc = draw.line(int(closest_cluster_pt_cent[1]),int(closest_cluster_pt_cent[0]),int(cy),int(cx))
                    c1_display_image_cent[rr, cc] = [255,0,255]
                
                if(c2_num_labels > 0 and c1_overlap_dict[label][1] == 0):
                    rr,cc = draw.line(int(closest_cluster_pt2[1]),int(closest_cluster_pt2[0]),int(closest_pt2[1]),int(closest_pt2[0]))
                    c2_display_image[rr, cc] = [0,0,255]
                    
                    rr,cc = draw.line(int(closest_cluster_pt_cent2[1]),int(closest_cluster_pt_cent2[0]),int(cy),int(cx))
                    c2_display_image_cent[rr, cc] = [255,0,255]
            
            c1_cluster_i += 1
            
        if(show_image):
            images = (c1_display_image, c2_display_image, c1_display_image_cent, c2_display_image_cent)
            image_names = ['edge','','centroid','']
            for image_i in (0,2):
                plt.imshow(images[image_i])
                for l in c1_labels_array:
                    plt.text(l[1],l[2],l[0],fontsize='1')
                SaveFigureAsImage(head + '/' + root + '-' + image_names[image_i] + '.pdf', plt.gcf(), orig_size=(int(images[image_i].shape[0]), int(images[image_i].shape[1])))
                plt.clf()
                plt.close()
                
                plt.imshow(images[image_i+1])
                for l in c1_labels_array:
                    plt.text(l[1],l[2],l[0],fontsize='1')
                    
                for prop in c2_props:
                    label = prop['Label']
                    bbox = prop['BoundingBox']
                    plt.text(bbox[3],bbox[2],str(label),fontsize='1')
                
                SaveFigureAsImage(head + '/' + root + '-' + image_names[image_i] + '-Green.pdf', plt.gcf(), orig_size=(int(images[image_i+1].shape[0]), int(images[image_i+1].shape[1])))
                plt.clf()
                plt.close()
        
        if(show_image):
            c1_display_image = color.gray2rgb(c1_image)
            c1_display_image_cent = color.gray2rgb(c1_image)
            c1_labels_array = []
            c2_display_image = color.gray2rgb(c2_image)
            c2_display_image_cent = color.gray2rgb(c2_image)
            c2_labels_array = []
            
        (head, tail) = os.path.split(image_file)
        (root, ext) = os.path.splitext(tail)
        for prop in c2_props:
            #cx,cy is coordinates of centroid of cluster
            cy,cx = prop['Centroid']
            label = prop['Label']
            area = prop['Area'] * pixels_to_nm * pixels_to_nm
            image = prop['Image'] #sliced binary region image same size as bounding box
            bbox = prop['BoundingBox']
            
            #calculate ImageJ perimeter/circularity:
            (manual_perim, image_perim) = find_ImageJ_perimeter(image)
            manual_circ = 4 * math.pi * ((area - .5*manual_perim) / manual_perim**2)
            if(manual_circ > 1.): manual_circ = 1.
            
            #get perimeter pixels
            perim_coords = []
            for i,val1 in enumerate(image_perim):
                for j,val2 in enumerate(val1):
                    if(val2):
                        #perim_coords.append([j,i])
                        perim_coords.append([j + bbox[1],i + bbox[0]]) #put back in whole image
                        if(show_image):
                            c2_display_image[i + bbox[0]][j + bbox[1]] = [0,255,0]
                            c1_display_image[i + bbox[0]][j + bbox[1]] = [0,255,0]
                            c2_display_image_cent[i + bbox[0]][j + bbox[1]] = [0,255,0]
                            c1_display_image_cent[i + bbox[0]][j + bbox[1]] = [0,255,0]
            
            #find closest cluster c2 to c2
            if(c2_num_labels > 1):
                #edge to edge               
                (min_distance, closest_cluster_pt, closest_pt) = closest_cluster(perim_coords, bbox, c2_mask)
                closest_cluster_label = c2_labeled_clusters[closest_cluster_pt[1]][closest_cluster_pt[0]]
                closest_cluster_area = c2_area_dict[closest_cluster_label]
                
                #centroid to centroid
                min_cent_distance = -1
                for prop_ in c2_props:
                    cy_,cx_ = prop_['Centroid']
                    label_ = prop_['Label']
                    if(label_ != label):
                        cent_distance = find_distance(cx,cy,cx_,cy_)
                        if(cent_distance < min_cent_distance or min_cent_distance == -1):
                            min_cent_distance = cent_distance
                            closest_cluster_label_cent = label_
                            closest_cluster_pt_cent = (cx_,cy_)
                            closest_cluster_area_cent = c2_area_dict[closest_cluster_label_cent]
                            
                min_distance *= pixels_to_nm
                min_cent_distance *= pixels_to_nm
                closest_cluster_area *= (pixels_to_nm**2)
                closest_cluster_area_cent *= (pixels_to_nm**2)
            else:
                min_distance = ''
                closest_cluster_label = ''
                min_cent_distance = ''
                closest_cluster_label_cent = ''
                closest_cluster_area = ''
                closest_cluster_area_cent = ''
            
            #find closest cluster c2 to c1
            if(c1_num_labels > 0): #must have atleast one c1 found
                if(c2_overlap_dict[label][1] == 0):
                    #edge to edge
                    (min_distance2, closest_cluster_pt2, closest_pt2) = closest_cluster(perim_coords, bbox, c2_mask, c1_mask)
                    closest_cluster_label2 = c1_labeled_clusters[closest_cluster_pt2[1]][closest_cluster_pt2[0]]
                    closest_cluster_area2 = c1_area_dict[closest_cluster_label2]
                    
                    #centroid to centroid
                    min_cent_distance2 = -1
                    for prop_ in c1_props:
                        cy_,cx_ = prop_['Centroid']
                        label_ = prop_['Label']
                        cent_distance = find_distance(cx,cy,cx_,cy_)
                        if(cent_distance < min_cent_distance2 or min_cent_distance2 == -1):
                            min_cent_distance2 = cent_distance
                            closest_cluster_label_cent2 = label_
                            closest_cluster_pt_cent2 = (cx_,cy_)
                            closest_cluster_area_cent2 = c1_area_dict[closest_cluster_label_cent2]
                else:
                    min_distance2 = 0
                    closest_cluster_label2 = ''
                    min_cent_distance2 = 0
                    closest_cluster_label_cent2 = ''
                    closest_cluster_area2 = 0
                    closest_cluster_area_cent2 = 0
            
                min_distance2 *= pixels_to_nm
                min_cent_distance2 *= pixels_to_nm
                closest_cluster_area2 *= (pixels_to_nm**2)
                closest_cluster_area_cent2 *= (pixels_to_nm**2)
            else:
                min_distance2 = ''
                closest_cluster_label2 = ''
                min_cent_distance2 = ''
                closest_cluster_label_cent2 = ''
                closest_cluster_area2 = ''
                closest_cluster_area_cent2 = ''
            
            #write to file
            c2_out_file.write(str(c2_cluster_i) + ',' + root + ',' + str(label) + ',' +
                            str(min_distance) + ',' + str(closest_cluster_label) + ',' + str(closest_cluster_area) + ',' + 
                            str(min_distance2) + ',' + str(closest_cluster_label2) + ',' + str(closest_cluster_area2) + ',' + 
                            str(min_cent_distance) + ',' + str(closest_cluster_label_cent) + ',' + str(closest_cluster_area_cent) + ',' + 
                            str(min_cent_distance2) + ',' + str(closest_cluster_label_cent2) + ',' + str(closest_cluster_area_cent2) + ',' + 
                            str(area) + ',' + str(cx * pixels_to_nm) + ',' + str(cy * pixels_to_nm) + ',' + str(prop['MajorAxisLength']*pixels_to_nm) + ',' +
                            str(prop['MinorAxisLength']*pixels_to_nm) + ',' + str(manual_perim) + ',' + str(prop['Eccentricity']) + ',' + str(manual_circ) + ',')
            #format overlap entries to print:
            ov_labels = [str(ov[0]) for ov in c2_overlap_dict[label][0]]
            ov_areas = [str(ov[1]*pixels_to_nm**2) for ov in c2_overlap_dict[label][0]]
            
            #areas of overlapping clusters (sum if more than one)
            area_ov_clusters = 0
            for ov in c2_overlap_dict[label][0]:
                area_ov_clusters += c1_area_dict[ov[0]]
            area_ov_clusters *= (pixels_to_nm**2)
            
            c2_out_file.write(';'.join(ov_labels) + ',' + ';'.join(ov_areas) + ',' + str(c2_overlap_dict[label][1]*pixels_to_nm**2) + ',' + str(area_ov_clusters) + '\n')
            
            if(show_image):
                #save label and x,y of the cluster for labeling the display image
                c2_labels_array.append([str(label), bbox[3], bbox[2]])
                
                #c2 to c2
                if(c2_num_labels > 1):
                    #edge distance
                    rr,cc = draw.line(int(closest_cluster_pt[1]),int(closest_cluster_pt[0]),int(closest_pt[1]),int(closest_pt[0]))
                    c2_display_image[rr, cc] = [255,0,0]
                    
                    #centroid distance
                    rr,cc = draw.line(int(closest_cluster_pt_cent[1]),int(closest_cluster_pt_cent[0]),int(cy),int(cx))
                    c2_display_image_cent[rr, cc] = [255,0,255]
                
                #c2 to c1
                if(c1_num_labels > 0 and c2_overlap_dict[label][1] == 0):
                    #edge distance
                    rr,cc = draw.line(int(closest_cluster_pt2[1]),int(closest_cluster_pt2[0]),int(closest_pt2[1]),int(closest_pt2[0]))
                    c1_display_image[rr, cc] = [0,0,255]
                    
                    #centroid distance
                    rr,cc = draw.line(int(closest_cluster_pt_cent2[1]),int(closest_cluster_pt_cent2[0]),int(cy),int(cx))
                    c1_display_image_cent[rr, cc] = [255,0,255]
            
            c2_cluster_i += 1
            
        if(show_image):
            images = (c1_display_image, c2_display_image, c1_display_image_cent, c2_display_image_cent)
            image_names = ['edge','','centroid','']
            for image_i in (0,2):
                plt.imshow(images[image_i+1])
                for l in c2_labels_array:
                    plt.text(l[1],l[2],l[0],fontsize='1')
                SaveFigureAsImage(head + '/' + root + '-' + image_names[image_i] + '.pdf', plt.gcf(), orig_size=(int(images[image_i+1].shape[0]), int(images[image_i+1].shape[1])))
                plt.clf()
                plt.close()
                
                plt.imshow(images[image_i])
                for l in c2_labels_array:
                    plt.text(l[1],l[2],l[0],fontsize='1')
                    
                for prop in c1_props:
                    label = prop['Label']
                    bbox = prop['BoundingBox']
                    plt.text(bbox[3],bbox[2],str(label),fontsize='1')
                
                SaveFigureAsImage(head + '/' + root + '-' + image_names[image_i] + '-Red.pdf', plt.gcf(), orig_size=(int(images[image_i].shape[0]), int(images[image_i].shape[1])))
                plt.clf()
                plt.close()
            
######################################################################################################################
#if calling program on the command line, the following can be un-commented

#cluster distances to ROI line - type = 1
#cluster distances - closest cluster - type = 2

if len(sys.argv) == 1:
    base_dir = 'E:/Sarah_2018/Rothenberg_archive/_CellsForAnalysis-anonymous/testpython' #Cre_1'
    run_type = 1
else:
    base_dir = sys.argv[1]
    run_type = int(sys.argv[2])
    
if(run_type == 1):
    cluster_distances(base_dir)
elif(run_type == 2):
    cluster_distances2(base_dir)
else:
    print "Invalid run type."


#cluster_distances(base_dir)

#base_dir = 'C:\\Projects\\Rothenberg\\Simone\\10_24_2014\\WT Masks' #R33Q Masks' #WT Masks' #'C:\\Projects\\Rothenberg\\Delmar\\Data\\Esperanza - 9_30_2014\\tissue_WT_Nav1.5-NCad\\20140410' #20140410' #20140226' #20131210' #20131023'
#base_dir = 'C:\\Projects\\Rothenberg\\Delmar\\Data\\Esperanza - 12_2_2014\\HL_1_cells_peptide\20141021_HL1_peptide_Nav1.5\\High stain Scramble'

#base_dir = 'C:/Projects/Rothenberg/Delmar/Data/Esperanza - 3_11_2014/D378_Nav1.5-Ncad/pkp2Het_1'
#base_dir = 'C:/Projects/Rothenberg/Delmar/Data/Esperanza - 3_11_2014/D378_Nav1.5-Ncad/wt_1'
#base_dir = 'C:/Projects/Rothenberg/Delmar/Data/Esperanza - 3_11_2014/D378_Nav1.5-Ncad/wt_2'
#base_dir = 'C:/Projects/Rothenberg/Delmar/Data/Esperanza - 5_19_2014\D378_AnkG/D378_1' #Cre_1'
#base_dir = 'C:\\Projects\\Rothenberg\\Delmar\\Data\\Esperanza - 12_5_2014\\TAC_study\\Sham 21_9-23-2014' #TAC 8_9-23-2014' #TAC 8_9-22-2014' #Sham 21_9-23-2014'

#base_dir = 'C:\\Projects\\Rothenberg\\Delmar\\Data\\Esperanza - 9_30_2014\\tissue_WT_Nav1.5-NCad\\20140226'
#base_dir = 'C:\\Projects\\Rothenberg\\Delmar\\Data\\Esperanza - 12_12_2014\\correlative_tissue_WT_Nav1.5-NCad'
#base_dir = 'C:\\Projects\\Rothenberg\\Delmar\\Data\\Esperanza - 12_12_2014\\Simulations\\1'
#cluster_distances2(base_dir)
#cluster_distances(base_dir)
