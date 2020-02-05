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

import pandas as pd
import glob
import os
import numpy as np
import sys

def fix_cols_order(df, color):
    cols = df.columns.tolist()
    cols.insert(0, cols.pop(cols.index('CentroidVolumeClosestRed')))
    cols.insert(0, cols.pop(cols.index('CentroidDistanceClosestRed')))
    cols.insert(0, cols.pop(cols.index('CentroidLabelClosestRed')))
    if(color == 'RED'):
        cols.insert(0, cols.pop(cols.index('CentroidVolumeClosestOverlapRed')))
        cols.insert(0, cols.pop(cols.index('CentroidDistanceClosestOverlapRed')))
        cols.insert(0, cols.pop(cols.index('CentroidLabelClosestOverlapRed')))
        
    cols.insert(0, cols.pop(cols.index('CentroidVolumeClosestGreen')))
    cols.insert(0, cols.pop(cols.index('CentroidDistanceClosestGreen')))
    cols.insert(0, cols.pop(cols.index('CentroidLabelClosestGreen')))
    if(color == 'GREEN'):
        cols.insert(0, cols.pop(cols.index('CentroidVolumeClosestOverlapGreen')))
        cols.insert(0, cols.pop(cols.index('CentroidDistanceClosestOverlapGreen')))
        cols.insert(0, cols.pop(cols.index('CentroidLabelClosestOverlapGreen')))
    
    cols.insert(0, cols.pop(cols.index('EdgeVolumeClosestRed')))
    cols.insert(0, cols.pop(cols.index('EdgeDistanceClosestRed')))
    cols.insert(0, cols.pop(cols.index('EdgeLabelClosestRed')))
    if(color == 'RED'):
        cols.insert(0, cols.pop(cols.index('EdgeVolumeClosestOverlapRed')))
        cols.insert(0, cols.pop(cols.index('EdgeDistanceClosestOverlapRed')))
        cols.insert(0, cols.pop(cols.index('EdgeLabelClosestOverlapRed')))
    
    cols.insert(0, cols.pop(cols.index('EdgeVolumeClosestGreen')))
    cols.insert(0, cols.pop(cols.index('EdgeDistanceClosestGreen')))
    cols.insert(0, cols.pop(cols.index('EdgeLabelClosestGreen')))
    if(color == 'GREEN'):
        cols.insert(0, cols.pop(cols.index('EdgeVolumeClosestOverlapGreen')))
        cols.insert(0, cols.pop(cols.index('EdgeDistanceClosestOverlapGreen')))
        cols.insert(0, cols.pop(cols.index('EdgeLabelClosestOverlapGreen')))
    
    cols.insert(0, cols.pop(cols.index('VolumeOverlapClusters')))
    cols.insert(0, cols.pop(cols.index('TotalOverlapVolume')))
    cols.insert(0, cols.pop(cols.index('OverlapVolumeList')))
    cols.insert(0, cols.pop(cols.index('OverlapSegmentList')))
    cols.insert(0, cols.pop(cols.index('Volume (nm^3)')))
    cols.insert(0, cols.pop(cols.index('Label')))
    cols.insert(0, cols.pop(cols.index('File')))
    cols.insert(0, cols.pop(cols.index('Dir')))
    return df[cols] #apply re-ordering of cols

def fill_volume_cols(df_1, color_1, df_2, color_2):
    #fill green volume columns
    df_1['VolumeOverlapClusters'] = ''
    df_1['CentroidVolumeClosest' + color_1] = ''
    df_1['CentroidVolumeClosest' + color_2] = ''
    df_1['EdgeVolumeClosest' + color_1] = ''
    df_1['EdgeVolumeClosest' + color_2] = ''
    
    df_1['CentroidVolumeClosestOverlap' + color_1] = ''
    df_1['EdgeVolumeClosestOverlap' + color_1] = ''
    
    for i in df_1.index:
        if((not np.isnan(df_1['TotalOverlapVolume'][i])) and df_1['TotalOverlapVolume'][i] > 0):
            #overlapping red
            ov_labels = df_1['OverlapSegmentList'][i]
            ov_labels = str(ov_labels) #it's float type if there's only one label, so convert
            ov_labels = ov_labels.split(';')
            #add volumes of all labels if it overlaps more than one
            vol = 0
            for label in ov_labels:
                label = float(label)
                label = int(label)
                vol += df_2['Volume (nm^3)'][label-1] #index is one less than label
            df_1['VolumeOverlapClusters'][i] = vol
            
            #get volume of closest c1 that is overlapping with a c2
            label = df_1['CentroidLabelClosestOverlap' + color_1][i]
            if(not np.isnan(label)):
                df_1['CentroidVolumeClosestOverlap' + color_1][i] = df_1['Volume (nm^3)'][label-1]
            label = df_1['EdgeLabelClosestOverlap' + color_1][i]
            if(not np.isnan(label)):
                df_1['EdgeVolumeClosestOverlap' + color_1][i] = df_1['Volume (nm^3)'][label-1]
        else:
            #closest c2 - centroid and edge
            label = df_1['CentroidLabelClosest' + color_2][i]
            if(not np.isnan(label)):
                df_1['CentroidVolumeClosest' + color_2][i] = df_2['Volume (nm^3)'][label-1]
            label = df_1['EdgeLabelClosest' + color_2][i]
            if(not np.isnan(label)):
                df_1['EdgeVolumeClosest' + color_2][i] = df_2['Volume (nm^3)'][label-1]
        #closest c1
        label = df_1['CentroidLabelClosest' + color_1][i]
        if(not np.isnan(label)):
            df_1['CentroidVolumeClosest' + color_1][i] = df_1['Volume (nm^3)'][label-1]
        label = df_1['EdgeLabelClosest' + color_1][i]
        if(not np.isnan(label)):
            df_1['EdgeVolumeClosest' + color_1][i] = df_1['Volume (nm^3)'][label-1]
            
def reprocess_3D(start_dir):
    df_all_red = pd.DataFrame()
    df_all_green = pd.DataFrame()
    
    #read in red and green table, for all files in all directories in the start_dir
    #find volumes of colocalized and/or closest cluster (red and green) add these as columns to the data frame
    #keep one big data frame with all the information, for all files
    
    for root, dirs, files in os.walk(start_dir):
        if(os.path.basename(root).startswith("_")): continue
        file_list = glob.glob(root + '/' + '*_More.csv')
        for file_name in file_list:
            (head, tail) = os.path.split(file_name)
            (root, ext) = os.path.splitext(tail)
            if('GREEN' in root):
                #green file
                cur_green_df = pd.read_csv(file_name, index_col=0)
                
                #red file
                red_file_name = tail.replace('GREEN', 'RED', 1)
                cur_red_df = pd.read_csv(head + '/' + red_file_name, index_col=0)
                
                fill_volume_cols(cur_green_df, 'Green', cur_red_df, 'Red')
                fill_volume_cols(cur_red_df, 'Red', cur_green_df, 'Green')
                
                cur_green_df['Dir'] = os.path.basename(head)
                cur_red_df['Dir'] = os.path.basename(head)
                cur_green_df['File'] = tail
                cur_red_df['File'] = red_file_name
                
                df_all_green = pd.concat([df_all_green,cur_green_df])  
                df_all_red = pd.concat([df_all_red,cur_red_df])
                
    #save df's with all the data from all files in the folder
    
    #fix columns
    #green
    df_all_green.rename(columns={' ': 'Label'}, inplace=True)
    df_all_green = fix_cols_order(df_all_green, 'GREEN')
    
    #red
    df_all_red.rename(columns={' ': 'Label'}, inplace=True)
    df_all_red = fix_cols_order(df_all_red, 'RED')
    
    df_all_green.to_csv(start_dir + '/Green_Mask_Measurements.csv', index=False)
    df_all_red.to_csv(start_dir + '/Red_Mask_Measurements.csv', index=False)
    
    
if len(sys.argv) == 1:
    start_dir = "C:\Projects\Rothenberg\Delmar\Data\3D Data\Analysis_3D_adult_cardiomyocytes-NO_MIN_VOL\test"
    #'C:\\Projects\\Rothenberg\\Delmar\\Data\\3D Data\\Analysis_3D_adult_cardiomyocytes'
else:
    start_dir = sys.argv[1]
    
#start_dir = 'C:\\Projects\\Rothenberg\\Delmar\\Data\\3D Data\\Analysis_3D_adult_cardiomyocytes' #20140804_3D_iPSCs_DelmarSetup' #20140810_3D_iPSCs_DelmarSetup' #Analysis_3D_adult_cardiomyocytes'
reprocess_3D(start_dir)

#start_dir = ('C:\\Projects\\Rothenberg\\Delmar\\Data\\3D Data\\Analysis_3D_adult_cardiomyocytes\\min_max_runs\\min-150000_max-1000000',)

#start_dir = ('C:\\Projects\\Rothenberg\\Delmar\\Data\\3D Data\\Analysis_3D_adult_cardiomyocytes\\min_max_runs\\min-400000_nomax',
#             'C:\\Projects\\Rothenberg\\Delmar\\Data\\3D Data\\Analysis_3D_adult_cardiomyocytes\\min_max_runs\\min-1000000_nomax',
#             'C:\\Projects\\Rothenberg\\Delmar\\Data\\3D Data\\Analysis_3D_adult_cardiomyocytes\\min_max_runs\\nomin_max-400000',
#             'C:\\Projects\\Rothenberg\\Delmar\\Data\\3D Data\\Analysis_3D_adult_cardiomyocytes\\min_max_runs\\nomin_max-1000000')
#for cur_dir in start_dir:    
#    reprocess_3D(cur_dir)

