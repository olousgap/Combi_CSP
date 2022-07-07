# Combi_CSP
CombiCSP is an open source software for dynamic modelling of concentrating solar energy power plants. CombiCSP utilizes solar resource, system engineering inputs as well as financial tools to provide dynamic simulations and annual yields of concentrating solar power plants. It readily provides modelling of plants based on solar power tower and parabolic trough collectors and it can be extended to novel solar energy modeling approaches and analyses as needed.
The relevant work pertaining to this software is : https://doi.org/10.1016/j.apenergy.2022.119450

# Quick-Start

the `Combi_CSP.ipynb` generates a graphs using functions from the code. 

Additionally, example cases are scripted in the following files:
- CSPCret.py: heliostat and CSP power and energy outputs in a location in Crete, Greece
- CSP50Compare.py: Combined heliostat and CSP power and energy outputs in a location in Crete, Greece


# installation/setup

## requirements

The following mainstream packages are required for this library.
- matplotlib
- scipy
- pandas
- ipykernel

Additional libraries are:

- pvlib_python 
- iapws (The InternationalAssociation for the Properties of Water and Steam) 
- numpy-financial 

## Conda installation

The following describes a minimum environment using conda. 

(Optional) Preferably create a new environment for the packages

`conda create -n combicsp python=3"

installation requires:
- matplotlib
- scipy
- pandas
- ipykernel

`conda install matplotlib scipy pandas ipykernel`

- pvlib_python:


`conda install -c conda-forge pvlib-python`

- iapws: The InternationalAssociation for the Properties of Water and Steam

`conda install -c conda-forge iapws`


- numpy-financial 

`conda install -c conda-forge numpy-financial`
