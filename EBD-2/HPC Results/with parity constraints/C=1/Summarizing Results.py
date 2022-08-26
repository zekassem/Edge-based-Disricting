import numpy as np
import csv
## End of Data
# Input parameters to the districting problem
no_dist=[2,4,6,8,10] # the number of districts
tol=[0.1,0.2,0.3]# tolerance parameter into the balancing constraints
tol_2=[0.1,0.2,0.3,0.5]
probs=[122,130]

with open('Results/Summary_1st_Logic.csv', 'w') as newFile:
    newFileWriter = csv.writer(newFile, lineterminator='\n')
    newFileWriter.writerow(
        ['Dataset', 'No. of Nodes', 'No. of Edges', 'No. of Districts', 't_1', 't_2', 'Time wo Cuts', 'Time w Cuts',
         'Improvement in Time (1st_Logic)', 'Objective Function wo Cuts', 'Objective w Cuts', 'Gap'])
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
                for t_2 in tol_2:
                    results_file = '1st_Logic_no_dist_'+str(num)+'_tol_1_'+str(t)+'_tol_2_'+str(t_2)+'_prob_'+str(prob)+'.csv'
                    time_wo_cuts = np.loadtxt(results_file, skiprows=1, max_rows=1, usecols=0, dtype=float, delimiter=',')
                    obj_wo_cuts = np.loadtxt(results_file, skiprows=1, max_rows=1, usecols=1, dtype=float, delimiter=',')
                    time_w_cuts = np.loadtxt(results_file, skiprows=5, max_rows=1, usecols=0, dtype=float, delimiter=',')
                    obj_w_cuts = np.loadtxt(results_file, skiprows=5, max_rows=1, usecols=1, dtype=float, delimiter=',')
                    gap = ((obj_w_cuts - obj_wo_cuts) / (obj_wo_cuts))  # calculating gap percentage
                    speedup = time_wo_cuts / time_w_cuts  # calculating xspeedup
                    newFileWriter.writerow([i, no_nodes1, no_edges1, num, t,t_2, time_wo_cuts, time_w_cuts, speedup, obj_wo_cuts, obj_w_cuts, gap])
                    i+=1










