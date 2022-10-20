# -*- coding: utf-8 -*-
"""
Optimization with analytical formulas for groundwater thermal plumes
- steady case

@author: Smajil Halilovic
"""

import numpy as np
import csv
from mip import Model, xsum, maximize, BINARY
from support_functions import *

# Reference point (center) for the coordinate system
ref_pt_x = 679700.0
ref_pt_y = 5335550.0
# Direction of groundwater flow - angle
direction = 90-31.9  # deg
direction = direction*np.pi/180 # rad
# Parameters required in LAHM - dictionary:
params_LAHM = {
    "n" : 0.3,  # [-]
    "b" : 8.51,  # [m]
    "v_a" : 4 * 10 ** (-5),  # [m/s]
    "alpha_L" : 5,  # [m]
    "alpha_T" : 0.5,  # [m]
    "Cw" : 4.185 * 10 ** 6,  # [J/Km**3]
    "Cm" : 2.888 * 10 ** 6,  # [J/Km**3]
    "R" : 2.3, # = Cm / (n * Cw)  # [-]
    "DT_inj" : 5  # [K]
}

# Define paths for input files:
file_pump_steady = 'GWHP_data/pump_rate_per_plot.csv'
file_well_positions = 'GWHP_data/well_positions_dist_5m.csv'

## Read data from input files
# Read the pumping rates of all regions (parcels)
regions_q = {}
with open(file_pump_steady, 'r') as file:
    reader = csv.reader(file)
    next(reader) # skip the first row
    for row in reader:
        id = row[0]
        q_demand = row[5]
        regions_q[id] = 2*float(q_demand)* 10 ** (-3) # [m3/s] # multiplied with 2 for winter case!

# Read positions of wells
well_ext_ids = {}
well_inj_ids = {}
with open(file_well_positions, 'r') as file:
    reader = csv.reader(file)
    next(reader)  # skip the first row
    x_inj = []
    y_inj = []
    x_ext = []
    y_ext = []
    q_rates = []
    for row in reader:
        id = row[0]
        x_coord = float(row[2])-ref_pt_x
        y_coord = float(row[3])-ref_pt_y
        # rotate the coordinate system
        x_coord, y_coord = rotation(x_coord, y_coord, direction)
        if row[1]=="extraction":
            # extraction wells
            # lists of coordinates
            x_ext.append(x_coord)
            y_ext.append(y_coord)
            # dictionary that relates regions and wells
            if id in well_ext_ids:
                well_ext_ids[id].append(x_ext.index(x_coord))
            else:
                well_ext_ids[id] = []
                well_ext_ids[id].append(x_ext.index(x_coord))
        else:
            # injection wells
            # lists of coordinates
            x_inj.append(x_coord)
            y_inj.append(y_coord)
            # list of pumping rates
            q_rates.append(regions_q[id])
            # dictionary that relates regions and wells
            if id in well_inj_ids:
                well_inj_ids[id].append(x_inj.index(x_coord))
            else:
                well_inj_ids[id] = []
                well_inj_ids[id].append(x_inj.index(x_coord))

# Indices for potential wells
Inj = range(len(x_inj))
print("Nr. of potential injection wells:")
print(len(x_inj))
Ext = range(len(x_ext))
print("Nr. of potential extraction wells:")
print(len(x_ext))

# Current time [s]
t_days = 120 # days (4 months)
t_current = t_days*24*60*60 # seconds

# Generate optimization model
m = Model(solver_name="GRB") # Gurobi solver is used!

# Add optimization variables
di = [m.add_var(var_type=BINARY) for i in Inj]
de = [m.add_var(var_type=BINARY) for e in Ext]

## Objective function
# All combinations for E_extracted
pe = [E_extracted(q_rates[j], t_current, params_LAHM) for j in Inj]
# maximize energy extracted from groundwater, i.e. its thermal potential:
m.objective = maximize(xsum(pe[k] * di[k] for k in Inj))

## Optimization constraints
# The number of wells per region is set to max. 1
for key in well_ext_ids:
    m += xsum(de[i] for i in well_ext_ids[key]) <= 1
# Number of extraction and injection wells is equal per region
for key in well_inj_ids:
    m += xsum(di[i] for i in well_inj_ids[key]) == xsum(de[i] for i in well_ext_ids[key])

# Constraint about minimum distance between installed extraction and injection wells on each plot
for key in well_ext_ids:
    for i in well_ext_ids[key]:
        for j in well_inj_ids[key]:
            dist_ij = distance_pts(x_ext[i], y_ext[i], x_inj[j], y_inj[j])
            if dist_ij < 10:
                m += de[i] + di[j] <= 1

# Constraint about negative thermal interference: DT_ext <= 1 [K]
for i in Ext:
    dT_ext = [delta_temp(x_ext[i]-x_inj[j], y_ext[i]-y_inj[j], q_rates[j], t_current, params_LAHM) for j in Inj]
    m += xsum(di[j] * dT_ext[j] for j in Inj) <= (-99 * de[i] + 100)

# Solve the optimization problem
m.optimize()
print('Optimal solution cost {} found'.format(m.objective_value))

# Read results from the optimal solution:
opt_vars = [v.x for v in m.vars]
opt_i = opt_vars[:len(x_inj)]
opt_di = [int(round(k)) for k in opt_i]
opt_e = opt_vars[len(x_inj):len(x_inj)+len(x_ext)]
opt_de = [int(round(k)) for k in opt_e]
opt_xe = [opt_de[i]*x_ext[i] for i in Ext if opt_de[i]>0]
opt_ye = [opt_de[i]*y_ext[i] for i in Ext if opt_de[i]>0]
opt_xi = [opt_di[i]*x_inj[i] for i in Inj if opt_di[i]>0]
opt_yi = [opt_di[i]*y_inj[i] for i in Inj if opt_di[i]>0]

# Print results:
print("Number of optimal extraction wells:", len(opt_xe))
print("Optimal extraction wells - (x, y):")
for i in range(len(opt_xe)):
    print('(',opt_xe[i], ',', opt_ye[i], ')')
print("Number of optimal injection wells:", len(opt_xi))
print("Optimal injection wells - (x, y):")
for i in range(len(opt_xi)):
    print('(',opt_xi[i], ',', opt_yi[i], ')')

# Save results:
with open(file_well_positions, 'r') as file:
    reader = csv.reader(file)
    next(reader) # skip the first row
    opt_wells = []
    for row in reader:
        id = row[0]
        x_coord = float(row[2]) - ref_pt_x
        y_coord = float(row[3]) - ref_pt_y
        # rotate the coordinate system
        x_coord, y_coord = rotation(x_coord, y_coord, direction)
        if row[1]=="extraction":
            # lists of coordinates
            opt_wells.append(opt_de[x_ext.index(x_coord)])
        else:
            # injection wells
            # lists of coordinates
            opt_wells.append(opt_di[x_inj.index(x_coord)])
# resulting list of installed (optimal) wells:
installed_wells = ['installed']
installed_wells.extend(opt_wells)
# add list as column to csv:
original_well_positions_file = file_well_positions
optimal_well_positions_file = file_well_positions[:-4] + "_optimal_steady.csv"
add_column_in_csv(original_well_positions_file, optimal_well_positions_file, lambda row, line_num: row.append(installed_wells[line_num-1]))
