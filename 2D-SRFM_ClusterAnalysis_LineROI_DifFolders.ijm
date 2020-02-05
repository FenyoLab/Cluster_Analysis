///    Copyright (C) 2020  Esperanza Agullo-Pascual
///
///    This program is free software: you can redistribute it and/or modify
///    it under the terms of the GNU General Public License as published by
///    the Free Software Foundation, either version 3 of the License, or
///    (at your option) any later version.
///
///    This program is distributed in the hope that it will be useful,
///    but WITHOUT ANY WARRANTY; without even the implied warranty of
///    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
///    GNU General Public License for more details.
///
///    You should have received a copy of the GNU General Public License
///    along with this program.  If not, see <http://www.gnu.org/licenses/>.
///

///Analizes particles in 2D-SRFM (2color) images given an ROI file or a zip file with different ROIs
///Input: folder containing folder/s with tif images (2D-SRFM green/red) and roi/zip file with ROIs for each image (same name with .roi or .zip extension)NEEDS a LINE ROI, not a polygon
///*** change the minimum particle size for analysis if needed -> LINE 13
///Output:  - for each image -> particles table (red and green), mask of particles detected, zip file with particle's ROI. 
///			- summary of analysis for each image: number clusters, total area, % area

macro "2D-SRFM_ClusterAnalysis" {

	///Set measurements for analysis
	run("Set Measurements...", "area centroid center perimeter bounding fit shape feret's display redirect=None decimal=3");

	///Set minimum particle size for analysis
	min = 1000;
	
	ProcessingDir = getDirectory("Choose Directory with folders to process");

	//1)Get folders list and loop through the folders
	folders = getFileList(ProcessingDir);
	 
    for (f=0; f<folders.length; f++) {
    	if (endsWith(folders[f], "/")){
		
		FolderToProcess = folders[f];
		folder_name = substring(FolderToProcess, 0, lengthOf(FolderToProcess)-1);
		
		print("Folder to process is "+folder_name);

		////
		////Analysis for each folder////
		////
		
		//2)Get list of images to analyze in the folder
		
		dir = ProcessingDir+FolderToProcess;
	
		list = getFileList(dir);

		///GREEN CHANNEL
	    for (i=0; i<list.length; i++) {
	    	
	    	if (endsWith(list[i], ".roi") || endsWith(list[i], ".zip")){
	    		
			    //3)Get RGB image name
				
				image = list[i];
				print(image);
				len = lengthOf(image);
				image_name = substring(image, 0, len-4); //substring 
				print(image_name);
				
				//4)Open RGB image
			
				open(dir+"//"+image_name + ".tif");
		
				//5)Image processing / particle analysis
		
				run("Split Channels");
				selectWindow(image_name + ".tif (blue)");
				close();
				selectWindow(image_name + ".tif (red)");
				close();
				selectWindow(image_name + ".tif (green)");
				setAutoThreshold("Default dark");
				
				//6)Open .zip or .roi file
				if (File.exists(dir+image_name+".zip")){
				roiManager("Open", dir+image_name+".zip");
	    		/*
	    		array1 = newArray("0");; 
				for (k=0;k<roiManager("count");k++){ 
	        		array1 = Array.concat(array1,k); 
	        		Array.print(array1);
				roiManager("select", array1);
				run("Line to Area");
	    		*/
	    		roiManager("Combine");
					} 
				
				else {
					roiManager("Open", dir+image_name+".roi");
					roiManager("select", 0);
					run("Line to Area");	
				}
				
				run("Analyze Particles...", "size="+min+"-Infinity circularity=0.00-1.00 show=Masks display exclude summarize add");
				selectWindow("Mask of " + image_name + ".tif (green)");

				//7)Save Mask and Results
				saveAs("Tiff", dir+image_name+"_Green_Mask.tif");
				close();
				selectWindow("Results");
		    	saveAs("Results", dir + image_name + "_Green_Clusters.xls");
				run("Close");
				
				selectWindow("ROI Manager");
		    	run("Close");
		    	selectWindow(image_name + ".tif (green)");
		    	close();
		    	}}
		    	
		    	//8)Save Summary with all data from different images (GREEN only)
				selectWindow("Summary");
		    	saveAs("Results", dir + folder_name + "_Green_Clusters_Summary.xls");
		    	run("Close");

		////RED CHANNEL	
		for (i=0; i<list.length; i++) {
	    	
	    	if (endsWith(list[i], ".roi") || endsWith(list[i], ".zip")){
	    		
			    //3)Get RGB image name
				
				image = list[i];
				print(image);
				len = lengthOf(image);
				image_name = substring(image, 0, len-4); //substring 
				print(image_name);
				
				//4)Open RGB image
			
				open(dir+"//"+image_name + ".tif");
		
				//5)Image processing / particle analysis
		
				run("Split Channels");
				selectWindow(image_name + ".tif (blue)");
				close();
				selectWindow(image_name + ".tif (green)");
				close();
				selectWindow(image_name + ".tif (red)");
				setAutoThreshold("Default dark");
				
				//6)Open .zip or .roi file
				if (File.exists(dir+image_name+".zip")){
				roiManager("Open", dir+image_name+".zip");
	    		/*
	    		array1 = newArray("0");; 
				for (k=0;k<roiManager("count");k++){ 
	        		array1 = Array.concat(array1,k); 
	        		Array.print(array1);
				roiManager("select", array1);
				//run("Line to Area");
	    		*/
	    		roiManager("Combine");
					}
				
				else {
					roiManager("Open", dir+image_name+".roi");
					roiManager("select", 0);
					run("Line to Area");	
				}
				
				run("Analyze Particles...", "size="+min+"-Infinity circularity=0.00-1.00 show=Masks display exclude summarize add");
				selectWindow("Mask of " + image_name + ".tif (red)");

				//7)Save Mask and Results
				saveAs("Tiff", dir+image_name+"_Red_Mask.tif");
				close();
				selectWindow("Results");
		    	saveAs("Results", dir + image_name + "_Red_Clusters.xls");
				run("Close");
				
				selectWindow("ROI Manager");
		    	run("Close");
		    	selectWindow(image_name + ".tif (red)");
		    	close();
		    	}}
		    	
		    	//8)Save Summary with all data from different images (GREEN only)
				selectWindow("Summary");
		    	saveAs("Results", dir + folder_name + "_Red_Clusters_Summary.xls");
		    	run("Close");

    	}}
    	print("MACRO FINISHED");
}