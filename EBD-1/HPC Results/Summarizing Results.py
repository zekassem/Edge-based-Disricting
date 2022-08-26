import time
import numpy as np
import csv
## End of Data
# Input parameters to the districting problem
# Input parameters to the districting problem
no_dist=[2,4,6,8,10] # the number of districts
tol=[0.08,0.1,0.2,0.3,1]# tolerance parameter into the balancing constraints
probs=[87,108,113,122]
with open('Results/Summary_Exact.csv', 'w') as newFile:
    newFileWriter = csv.writer(newFile, lineterminator='\n')
    newFileWriter.writerow(
        ['Dataset', 'No. of Nodes', 'No. of Edges', 'No. of Districts', 'Tolerance', 'Time', 'Objective Function'])
    i=1
    for prob in probs:
        graph_file = 'data_' + str(prob) + '.csv'
        no_nodes = np.loadtxt(graph_file, skiprows=1, max_rows=1, usecols=0,
                              dtype=int)  # skipping lines to get the number of nodes in the planar graph
        no_nodes1 = int(no_nodes)
        coordinates = np.loadtxt(graph_file, skiprows=4, max_rows=no_nodes1, usecols=(1, 2),
                                 delimiter=',')  # skipping lines to get the x & y coordinates of the nodes
        s = no_nodes + 4 + 2  # Calculating the lines needed to be skipped to get the number of edges
        no_edges = np.loadtxt(graph_file, skiprows=s, max_rows=1, usecols=0,
                              dtype=int)  # skipping lines to get the
        no_edges1 = int(no_edges)
        for num in no_dist:
            for t in tol:
                results_file = 'Exact_no_dist_' + str(num) + '_tol_' + str(t) + '_prob_' + str(prob) + '.csv'
                time = np.loadtxt(results_file, skiprows=1, max_rows=1, usecols=0, dtype=float, delimiter=',')
                obj= np.loadtxt(results_file, skiprows=1, max_rows=1, usecols=1, dtype=float, delimiter=',')
                newFileWriter.writerow([i, no_nodes1, no_edges1, num, t, time, obj])
                i+=1










