# geos-tough

This repo contains scripts used for the coupled GEOS-TOUGH simulation.

## geos-converter

These scripts are used to generate a Python module that can read raw GEOS outputs. It requires silo and hdf5. Unfortunately, these scripts are not compatible with the latest version of Python (3.9)

## meshmaker

These scripts are used to generate an executable of the latest mesh maker (developed by Dr. George Moridis at LBNL) for TOUGH.

## tough-preprocessor

These scripts are used to generate input files for TOUGH simulations. 

### Manipulating GEOS outputs

```
script_geos2tough.py
```
Top-level GEOS-converting script. It outputs each GEOS field as an N by 4 array. N is the number of data points. The 4 columns are [x, y, z, s], where s could be any output variable (aperture, pressure, etc.). 

```
geos_convert.py
```
Contains functions that handle the conversion process. 

```
script_geos_compare.py
```
Plot GEOS aperture from different GEOS simulations with different post-processing steps (propped, unpropped, modified, etc.).
