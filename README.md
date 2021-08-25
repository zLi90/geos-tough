# geos-tough

This repo contains scripts used for the coupled GEOS-TOUGH simulation.

## geos-converter

These scripts are used to generate a Python module that can read raw GEOS outputs. It requires silo and hdf5. Unfortunately, these scripts are not compatible with the latest version of Python (3.9)

## meshmaker

These scripts are used to generate an executable of the latest mesh maker (developed by Dr. George Moridis at LBNL) for TOUGH.

## tough-preprocessor

These scripts are used to generate input files for TOUGH simulations. Some scripts call auxiliary functions that can be found in `mm_input.py`

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
Plot GEOS aperture from different GEOS simulations with different post-processing steps (propped, unpropped, modified, etc.). An example figure looks like this:

<img width="696" alt="Screen Shot 2021-08-25 at 7 44 14 PM" src="https://user-images.githubusercontent.com/18181974/130784951-b06a5e2c-8bdb-49c5-9fe5-9611f329a5c0.png">


```
script_sensntwk.py
```
Generate input files (INPUT, MESH, INCON) for TOUGH simulations. This entire automated process contains the following steps:
* Create an input file for mesh maker, then call the mesh maker to generate the MESH file. Options are available to create MESH with MINC.
* Insert GEOS outputs into the INCON file to represent non-uniform properties of the hydraulic fractures.
* Apply non-uniform permeability for the MINC fractures in INCON. This is used to represent the secondary fractures.
* Automatically write the INPUT file for TOUGH.




