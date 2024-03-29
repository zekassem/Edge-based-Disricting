import numpy as np
import csv
## End of Data
# Input parameters to the districting problem
no_dist=[2,4,6,8,10] # the number of districts
tol=[0.02,0.04,0.06]# tolerance parameter into the balancing constraints
probs=[87,108,113,122,130,149,151,158,171,224,350,400,500]

with open('Results/Summary_1st_Logic_2nd_wo_z_2.csv', 'w') as newFile:
    newFileWriter = csv.writer(newFile, lineterminator='\n')
    newFileWriter.writerow(
        ['Dataset', 'No. of Nodes', 'No. of Edges', 'No. of Districts', 'Tolerance', 'Time wo Cuts', 'Time w Cuts',
         'Objective Function wo Cuts', 'Objective w Cuts', 'No_Iterations_wo_logic', 'No_Cuts_wo_logic',
         'No_Iterations_w_logic', 'No_Cuts_w_logic'])
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
                results_file = '1st_Logic_2nd_wo_z_variables_no_dist_' + str(num) + '_tol_' + str(t) + '_prob_' + str(prob) + '.csv'
                try:
                    time_wo_cuts = np.loadtxt(results_file, skiprows=1, max_rows=1, usecols=0, dtype=float,
                                              delimiter=',')
                    obj_wo_cuts = np.loadtxt(results_file, skiprows=1, max_rows=1, usecols=1, dtype=float,
                                             delimiter=',')
                    No_Iterations_wo_logic = np.loadtxt(results_file, skiprows=3, max_rows=1, usecols=1, dtype=float,
                                                        delimiter=',')
                    No_cuts_wo_logic = np.loadtxt(results_file, skiprows=3, max_rows=1, usecols=2, dtype=float,
                                                  delimiter=',')
                    time_w_cuts = np.loadtxt(results_file, skiprows=5, max_rows=1, usecols=0, dtype=float,
                                             delimiter=',')
                    obj_w_cuts = np.loadtxt(results_file, skiprows=5, max_rows=1, usecols=1, dtype=float, delimiter=',')
                    No_Iterations_w_logic = np.loadtxt(results_file, skiprows=7, max_rows=1, usecols=1, dtype=float,
                                                       delimiter=',')
                    No_cuts_w_logic = np.loadtxt(results_file, skiprows=7, max_rows=1, usecols=2, dtype=float,
                                                 delimiter=',')
                    newFileWriter.writerow(
                        [i, no_nodes1, no_edges1, num, t, time_wo_cuts, time_w_cuts, obj_wo_cuts, obj_w_cuts,
                         No_Iterations_wo_logic, No_cuts_wo_logic, No_Iterations_w_logic, No_cuts_w_logic])
                    i+=1
                except:
                    time_wo_cuts = np.loadtxt(results_file, skiprows=1, max_rows=1, usecols=0, dtype=str,
                                              delimiter=',')
                    if time_wo_cuts == 'Infeasible':
                        newFileWriter.writerow(
                            [i, no_nodes1, no_edges1, num, t, 'Infeasible', 'Infeasible', 'Infeasible', 'Infeasible',
                             'Infeasible', 'Infeasible', 'Infeasible', 'Infeasible'])
                        i += 1
                    pass
                else:
                    pass










