BoneMeasurement
===============

Computes, combines and plots long bone end point measurements from
2D DXA scan images with a single Python script.

Measurements made by an expert on a model image are transferred
to assay images through region of interest extraction and
non-linear registration.
Computed measurements which consist of locations and confidence figures
(based on local image cross correlation) are output as a CSV file.
You can expect there to be points with large errors and points not found;
these depend on variations of pose, variations in intensity, missing limbs,
occlusion, etc.... One way to reduce errors and increase the number of
points found  is to make separate runs using different models and combine
these later using the confidence values to determine which measurements
are kept.

The measurements may also be plotted for visual assessment.

The pipeline used to process the scan is:

1. Pre-process images
   * convert to Woolz format from DICOM
   * remove background borders 
   * histogram match assay images to model image
   * create a smoothed image for finding ROI centres
   * convert to NIfTI format (for ANTs)
   * find ROI centres
     + scan from head to toes looking a 1D mage intensity profiles
     + save ROI centres to file
   * extract regions of interest about ROI centres
     + create mask
       - cut ROI image
       - copy and smooth image
       - apply spatial window function
       - threshold
       - perform connected component labeling
       - discard small and other regions not connected to largest
       - apply resulting mask to ROI image
     + apply Sobel filter to enhance edges
     + smooth
     + shift coordinates ROI to origin
2. Register assay ROIs to model ROIs using ANTs
3. Read precomputed measurements.
4. Map measurement locations from model to assays using ROI offset and ANTs
   transforms appending to measurements
5. Combine measurements keeping assay duplicates with highest confidence
   which are above the confidence threshold
6. Save measurements
7. Plot measurements

Many of these pipeline stages are optional and are controlled by the command
line arguments.  Run using -h or --help for usage.

This software relies on:

- _Python 2.7_ (including _numpy_, _pydicom_ and _matplotlib_)
- _ANTs_ for image registration - https://stnava.github.io/ANTs/
- _Imagemagick convert_ - http://www.imagemagick.org/
- _Woolz_ image processing - https://github.com/billhill00/Woolz
- _PyWoolz_ Python binding to Woolz - https://github.com/billhill00/PyWoolz

Example plot (excluding head):

![eample plot][exp]

[exp]: example_plot.png
