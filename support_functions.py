# -*- coding: utf-8 -*-
"""
Definitions of supporting functions

@author: Smajil Halilovic
"""

import numpy as np
import csv
from scipy import special

def add_column_in_csv(input_file, output_file, transform_row):
    """ Append a column in existing csv using csv.reader / csv.writer classes """
    # Open the input_file in read mode and output_file in write mode
    with open(input_file, 'r') as read_obj, \
            open(output_file, 'w', newline='') as write_obj:
        # Create a csv.reader object from the input file object
        csv_reader = csv.reader(read_obj)
        # Create a csv.writer object from the output file object
        csv_writer = csv.writer(write_obj)
        # Read each row of the input csv file as list
        for row in csv_reader:
            # Pass the list / row in the transform function to add column text for this row
            transform_row(row, csv_reader.line_num)
            # Write the updated row / list to the output file
            csv_writer.writerow(row)

def rotation(x, y, theta):
    """ Rotate the coordinate system """
    x_r = x*np.cos(theta) + y*np.sin(theta)
    y_r = -x*np.sin(theta) + y*np.cos(theta)
    return x_r, y_r

def distance_pts(x1, y1, x2, y2):
    """ Compute Euclidean distance between two points - 2D space """
    d = np.sqrt((x1-x2)**2+(y1-y2)**2)
    return d

def arrange_list(list_org, start_pos):
    """ Rearrange elements in list - needed to e.g. start simulation in August instead of January """
    # list_org: original list
    # start_pos: starting position for the new list: if 1 -> the same as old one
    # list_new: new rearranged list
    list_new = []
    ind = start_pos-1
    for i in range(len(list_org)):
        if ind < len(list_org):
            list_new.append(list_org[ind])
            ind += 1
        else:
            ind = 0
            list_new.append(list_org[ind])
            ind += 1
    return list_new

def delta_temp(x, y, q_inj, t, params_LAHM):
    """ Analytical formula for thermal plumes - LAHM """
    r = np.sqrt(x ** 2 + (y ** 2) * (params_LAHM["alpha_L"] / params_LAHM["alpha_T"]))
    ampl_term = 1 / (4 * params_LAHM["n"] * params_LAHM["b"] * params_LAHM["v_a"] * np.sqrt(np.pi * params_LAHM["alpha_T"]))
    exp_term = np.exp((x - r) / (2 * params_LAHM["alpha_L"]))
    erfc_term = special.erfc((r - params_LAHM["v_a"] * t / params_LAHM["R"]) / (2 * np.sqrt(params_LAHM["v_a"] * params_LAHM["alpha_L"] * t / params_LAHM["R"])))
    delta_T = q_inj * params_LAHM["DT_inj"] * ampl_term * exp_term * (1 / np.sqrt(r)) * erfc_term

    return delta_T

def E_extracted(q_inj, time_period, params_LAHM):
    """ Formula for the extracted energy from the aquifer - objective function """
    return q_inj*params_LAHM["Cw"]*params_LAHM["DT_inj"]*time_period
