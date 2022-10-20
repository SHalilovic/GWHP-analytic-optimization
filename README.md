# GWHP-analytic-optimization
Supporting code for the publication: Optimizing the spatial arrangement of groundwater heat pumps and their well locations

Installation: 
Download all the files from the repository and place them into one folder.

The files are the following:
- 'optimization_steady_case.py': main Python script for optimizing the GWHP locations and their well layouts in the steady case.
- 'optimization_dynamic_case.py': main Python script for optimizing the GWHP locations and their well layouts in the dynamic case.
- 'support_functions.py': Python script with the implementation of supporting functions. 
- 'GWHP_data' : folder that contains all the required input data 
    
Usage:
1) Install the required Python packages: "numpy", "csv" and "scipy".
2) Install Python-MIP package (https://www.python-mip.com/) and Gurobi solver (https://www.gurobi.com/documentation/quickstart.html) for the best performance. 
3) Run the script 'optimization_steady_case.py' or 'optimization_dynamic_case.py'.
4) The result files will be saved in the folder 'GWHP_data' as csv files with the ending "_optimal_steady" or "_optimal_dynamic". 
 
