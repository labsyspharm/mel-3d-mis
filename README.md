# Highly Multiplexed 3D Profiling of Cell States and Immune Niches in Human Tumours

**Authors:** Clarence Yapp, Ajit J. Nirmal, Felix Zhou, Alex Y.H. Wong, Juliann B. Tefft, Yi Daniel Lu, Zhiguo Shang, Zoltan Maliga, Paula Montero Llopis, George F Murphy, Christine G Lian, Gaudenz Danuser, Sandro Santagata, and Peter K. Sorger

*This respository ([labsyspharm/mel-3d-mis](https://github.com/labsyspharm/mel-3d-mis)) hosts original code associated with the above publication for 3D image registration, intensity quantification, cell shape & orientation analysis, collagen-to-cell distance quantification, and cell interaction analysis.  Further details about accessing the primary data associated with the paper will be published here after acceptance.*

## System Requirements
### Hardware requirements
+ Cell type interaction analysis, intensity quantification, and cell shape & orientation analysis require a standard computer with 8GB RAM.
+ 3D image registration and collagen-to-cell distance quantification require significant amounts of RAM (tested on 600-700GB of RAM on a high performance cluster) and depends on the size of the image ROI. Updated versions will be less RAM intensive and may also feature parallelization over tiles.

### Software requirements
MATLAB version 2021/2022 (or newer) on Windows 10 with the following toolboxes:
+ Curve Fitting Toolbox (tested on version 3.6)
+ Statistics and Machine Learning Toolbox (tested on version 12.2)
+ Image Processing Toolbox (tested on version 11.4)
+ Parallel Computing Toolbox (OPTIONAL tested on version 7.5)
Install time is typically 1 hour or less but will vary depending on network speed and computer specifications.

# Demo on running cell interaction analysis on two cells.
## Overview
This demo illustrates how to assess cell type interactions between neighboring cells by measuring membrane intensities along horizontal line intensity profiles. A CD4 and CD8 T cell are used as an example here and is illustrated in Figure 6p-r of corresponding manuscript. <br/> <img src="https://github.com/labsyspharm/mel-3d-mis/blob/main/Demo/overview.PNG" width="400">

## Instructions to run on included data
Download code and source image (Figure 6.tif) found in demo folder. The .tif file is a truncated version of Figure 6p from manuscript but with only 4 z planes and 2 channels (CD4 and CD8 respectively). Load the the CD4 and CD8 channels by running the first 2 lines. Choose either tight or loose interaction, and run the corresponding code block. Run remainder of script.
## Expected output and run time
When complete, the demo will output a graph with two intensity profiles and their fitted polynomial curves. The x-axis is in microns along intensity profile(s). The y-axis is the normalized intensity where 1 is the maximum intensity along the line profile. Red 'x' denote the peak of each polynomial curve and estimates the location of the membrane at sub-pixel accuracy. The separation between both peaks (denoted by red 'x') along x-axis is a measure of how 'tight' the membrane interaction is. For example, the distance between cell 2 and 3 is ~30nm and represents a tight interaction (type I). <br/>
<img src="https://github.com/labsyspharm/mel-3d-mis/blob/main/Demo/tight.PNG" width="400">

The distance between cell 4 and 5 is ~150nm. This is a type II interaction. <br/>
<img src="https://github.com/labsyspharm/mel-3d-mis/blob/main/Demo/loose.PNG" width="400"> <br/>
Expected run time will be on the order of seconds especially if data is being loaded from local disk. If data is on network, run time will be longer.


## Instructions to run on YOUR data
Image can be a multichannel 3D file or separated channels. Images should be generated with considerations outlined in Supplementary section of corresponding manuscript (ie. high signal-to-noise ratio, validated antibodies, etc).
1. update file paths to locations of channel A and B (line 1 and 2). 
2. images may need to be rotated in order to draw horizontal line intensity profiles. Change variable *rot_angle* to an angle that rotates the image so that your line profile will horizontally through the cell membrane of interest.
3. enter XY coordinates of the start of line profile (*xstart, ystart*) and the end (*xend, yend*). Select optimal Z plane index, *planeZ*.
4. update *pixelSize* to reflect the lateral pixel sizes (microns/pixel) in your image data.
5. (optional), change colors of curves following RGB convention on line 38.










