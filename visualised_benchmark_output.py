#!/usr/bin/python3
#===============================================================
# v1.0 - Initial version, Adrian Jackson
#===============================================================
#
#----------------------------------------------------------------------
# Copyright 2021 EPCC, The University of Edinburgh
#
# svd-benchmark is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# eigen-benchmark is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with svd-benchmark.  If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------
#

# This program parses the output for runs of the svd benchmark to produce graphs of the runtime
# compared to the matrix size and process decomposition for the experiments undertaken.
#
# This is a very fragile program as it is closely tied to both the output of the benchmark
# and the output of batch scripts to run the program. In particular it expects the eigen
# benchmark to output lines like this:
#
# Matrix size is 16 x 16
# Process grid is 2 x 2 (0 are idle)
# Blacs gather time is: 0.000182152
# Blacs scatter time is: 0.000117779
# ...
# Scalapack time is: 0.00342679
# Eigen time is: 0.000210047
#
# Extraneous lines in the output file processed by this program, or lines with differing formats,
# will likely cause this program to crash or work incorrectly.
# A more robust approach to this would be to use some self consistent output language for the
# benchmark data, such as XML or YAML.


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys


def readfile(filename):
    array = []
    temparray = None;
    gather_time = None;
    scatter_time = None
    f = open(filename, "r")
    for line in f:
        line = line.rstrip()
        if("Matrix size" in line):
            if not (temparray is None):
                total_time = gather_time + scatter_time
                temparray.append(total_time)
                array.append(temparray)
            temparray = []
            text = line.split(" ")
            temparray.append([text[3],text[5]])
        elif("Process grid" in line):
            text = line.split(" ")
            temparray.append([text[3],text[5]])
        elif("Blacs gather time" in line):
            text = line.split(" ")
            gather_time = float(text[4])
            temparray.append(text[4])
        elif("Blacs scatter time" in line):
            text = line.split(" ")
            scatter_time = float(text[4])
            temparray.append(text[4])
        elif("Scalapack time" in line):
            text = line.split(" ")
            temparray.append(text[3])
        elif("Eigen time" in line):
            text = line.split(" ")
            temparray.append(text[3])
            
    if not (temparray is None):
        total_time = gather_time + scatter_time
        temparray.append(total_time)
        array.append(temparray)

    return array
            

if(len(sys.argv) != 2):
    print("You need to provide the name of the file with the data in it to be processed")
    exit(0)

filename = sys.argv[1]

print("Processing " + filename)

data_array = readfile(filename)

output_data = []
matrix_size = None
data = None

for place, element in enumerate(data_array):
    temp_size = int(element[0][0])*int(element[0][1])

    if(matrix_size is None):
        matrix_size = temp_size
        data = []
    if(temp_size != matrix_size):
        output_data.append([matrix_size,data])
        matrix_size = temp_size
        data = []
    data.append([element[1:]])
    if(place == len(data_array)-1):
       output_data.append([matrix_size,data])

eigen_graphs = []
scalapack_graphs = []
blacs_graphs = []
x_points = []
for place in output_data:
    for element in place[1:]:
        for next_element in element:
            graph_name = next_element[0][0][0] + 'x' + next_element[0][0][1]
            found = False
            found_index = 0
            for index,graph_element in enumerate(scalapack_graphs):
                if(graph_element[0] == graph_name):
                    found = True
                    found_index = index
                    break
            if(not found):
                scalapack_graphs.append([graph_name,float(next_element[0][3])])
                if(len(next_element[0]) == 6):
                    eigen_graphs.append([graph_name,float(next_element[0][4])])
                blacs_graphs.append([graph_name,float(next_element[0][-1])])
            else:
                scalapack_graphs[found_index].append(float(next_element[0][3]))
                if(len(next_element[0]) == 6):
                    eigen_graphs[found_index].append(float(next_element[0][4]))
                blacs_graphs[found_index].append(float(next_element[0][-1]))
        x_points.append(place[0])

plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Set2.colors)
fig = plt.figure()
ax = plt.axes()
        
for elements in scalapack_graphs:
    label = 'Scalapack (' + elements[0] + ')'
    data_part = elements[1:]
    y_points = []
    error_bars = [[],[]]
    for first, second, third in zip(*[iter(data_part)]*3):
        arith_mean = (first +  second + third)/3
        y_points.append(arith_mean)
        error_bars[0].append(arith_mean-min(first, second, third))
        error_bars[1].append(max(first, second, third)-arith_mean)
    plt.errorbar(x_points, y_points, yerr=error_bars, label=label)


label = 'Eigen'
elements = eigen_graphs[0]
data_part = elements[1:]
y_points = []
error_bars = [[],[]]
for first, second, third in zip(*[iter(data_part)]*3):
    arith_mean = (first +  second + third)/3
    y_points.append(arith_mean)
    error_bars[0].append(arith_mean-min(first, second, third))
    error_bars[1].append(max(first, second, third)-arith_mean)
plt.errorbar(x_points, y_points, yerr=error_bars, label=label)

 
plt.yscale('log')
plt.xscale('log')
plt.ylabel('Runtime (seconds)')
plt.xlabel('Matrix size')
plt.xticks(rotation=90)
#ax.ticklabel_format(style='plain')
plt.legend(loc="upper left")

plt.tight_layout()

plt.savefig(filename+"_output.png",format="png")


    #    data_part = element[1:]
#    row_labels.append(element[0])
#    if(first_run):
#        x_points = np.array([item[0] for item in data_part], dtype=int)
#        first_run = False
#    runtimes = np.array([item[2] for item in data_part], dtype=float)
#    table.append(runtimes)


        
print(x_points)
print(scalapack_graphs)
print(eigen_graphs)
print(blacs_graphs)
#fig = plt.figure()
#ax = plt.axes()
    
#plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.Set2.colors)

#x_points = None
#first_run = True
#table = []
#row_labels = []
#rows = len(data_array)

#colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


#plt.clf()
#plt.cla()
#plt.close()
#ax = plt.axes()
#ax.axis('off')
#ax.axis('tight')
#table_pic = plt.table(table, cellColours=None, cellLoc='right', colWidths=None, rowLabels=row_labels, rowColours=colors, rowLoc='left', colLabels=x_points, colColours=['grey']*len(table[0]), colLoc='center', loc='center', edges='closed')

#table_pic.auto_set_font_size(False)
#table_pic.set_fontsize(10)
#table_pic.scale(2,2)

#plt.tight_layout()

#plt.savefig(filename+"_table_output.png",bbox_inches='tight',format="png")
