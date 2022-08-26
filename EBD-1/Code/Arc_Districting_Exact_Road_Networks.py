import cplex
import time
import numpy as np
import collections
import csv
import copy

# Loading Data
class Data():
    def __init__(self,file):
        self.file = file
        self.no_nodes = np.loadtxt(file, skiprows=1, max_rows=1, usecols=0,
                                   dtype=int)  # skipping lines to get the number of nodes in the planar graph
        self.no_nodes1 = int(self.no_nodes)
        self.coordinates = np.loadtxt(file, skiprows=4, max_rows=self.no_nodes1, usecols=(1, 2),
                                      delimiter=',')  # skipping lines to get the x & y coordinates of the nodes
        self.s = self.no_nodes + 4 + 2  # Calculating the lines needed to be skipped to get the number of edges
        self.no_edges = np.loadtxt(file, skiprows=self.s, max_rows=1, usecols=0, dtype=int)  # skipping lines to get the
        self.no_edges1 = int(self.no_edges)
        self.z = self.s + 4  # Calculating the lines needed to be skipped to get adjacency list
        self.adj_list = np.loadtxt(file, skiprows=self.z, max_rows=self.no_edges1, usecols=(0, 1), dtype=int,
                                   delimiter=',')
        self.node = list(range(1, (self.no_nodes) + 1))
        self.from_node_i = self.adj_list[:, 0].tolist()
        self.to_node_j = self.adj_list[:, 1].tolist()
        self.edges_v = list(range(1, self.no_edges + 1))  # Edges Vector

## End of Data
# Input parameters to the districting problem
no_dist=[10] # the number of districts
tol=[0.2]# tolerance parameter into the balancing constraints


def constraints_wo_cuts(c,edges_v,node,edges, no_nodes,no_edges,no_dist,t,x_v,w_v,s_v,f_v_n): #  This function adds all of the constraints for the original districting problem
    # sum i e V (xie)=1 for e e E: each edge is assigned to one territory
    sum_x = []
    coeff = []
    for e in edges_v:
        sum_x = x_v[(e - 1):no_edges * no_nodes:no_edges]
        coeff = [1 for y in range(no_nodes)]
        c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["E"], rhs=[1])
        sum_x = []
        coeff = []

    # sum w_ii=p
    sum_x = []
    coeff = []
    sum_x = w_v[0:no_nodes * no_nodes:no_nodes + 1]
    coeff = [1 for g in range(no_nodes)]
    c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["E"],
                             rhs=[no_dist])
    sum_x = []
    coeff = []

    # Connectivity Constraint # 1: sik<=wik for all i,k e V
    sum_x = []
    coeff = []
    s_variables = []
    for i in node:
        for k in node:
            sum_x.append(w_v[(i - 1) * no_nodes + (k - 1)])
            coeff.append(-1)
            s_variables = s_v[(i - 1) * no_nodes + (k - 1)]
            sum_x.append(s_variables)
            coeff.append(1)
            c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["L"],
                                     rhs=[0])
            sum_x = []
            coeff = []
            s_variables = []

    sum_x = []
    coeff = []
    # Balancing Constraints
    # sum e e E xie <= (|E|/p) * (1+tau) wii for i e V
    rhs_1 = (no_edges / no_dist) * (1 + t)

    for i in node:
        sum_x.append(w_v[(i - 1) * (no_nodes + 1)])
        coeff.append(-rhs_1)
        sum_x.extend(x_v[(i - 1) * no_edges:i * no_edges])
        coeff2 = [1 for q in range(no_edges)]
        coeff.extend(coeff2)
        c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["L"], rhs=[0])
        coeff = []
        sum_x = []

    # sum e e E xie >= (|E|/p) * (1-tau) wii for i e V
    sum_x = []
    coeff = []
    coeff2 = []
    rhs_2 = (no_edges / no_dist) * (1 - t)
    for i in node:
        sum_x.append(w_v[(i - 1) * (no_nodes + 1)])
        sum_x.extend(x_v[(i - 1) * no_edges:i * no_edges])
        coeff.append(-rhs_2)
        coeff2 = [1 for q in range(no_edges)]
        coeff.extend(coeff2)
        c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["G"], rhs=[0])

        coeff = []
        sum_x = []
    # sum k e V sik=wii for every district there is one sink
    sum_x = []
    coeff = []
    for i in node:
        sum_x.append(w_v[(i - 1) * (no_nodes + 1)])
        coeff.append(-1)
        sum_x.extend(s_v[(i - 1) * no_nodes:i * no_nodes])
        coeff2 = [1 for q in range(no_nodes)]
        coeff.extend(coeff2)
        c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["E"],
                                 rhs=[0])
        sum_x = []
        coeff = []

    sum_x = []
    coeff = []
    # sum (j,k) e delta(j) x_i,(j,k)<=|delta(j)| wij for i,j e V
    for i in node:
        for j in node:
            sum_x.append(w_v[(i - 1) * no_nodes + (j - 1)])
            incident_edges=edges[j]
            coeff.append(-len(incident_edges))
            indicies = [(q - 1) + (i - 1) * no_edges for q in incident_edges]
            x_variables = [x_v[i] for i in indicies]
            sum_x.extend(x_variables)
            coeff2 = [1 for r in range(len(incident_edges))]
            coeff.extend(coeff2)
            c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["L"],
                                     rhs=[0])
            sum_x=[]
            coeff=[]

    # sum (j,k) e delta(j) x_i,(j,k)>=wij for i,j e V
    for i in node:
        for j in node:
            sum_x.append(w_v[(i - 1) * no_nodes + (j - 1)])
            incident_edges=edges[j]
            coeff.append(-1)
            indicies = [(q - 1) + (i - 1) * no_edges for q in incident_edges]
            x_variables = [x_v[i] for i in indicies]
            sum_x.extend(x_variables)
            coeff2 = [1 for r in range(len(incident_edges))]
            coeff.extend(coeff2)
            c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["G"],
                                     rhs=[0])
            sum_x=[]
            coeff=[]


    # Calculating M (|E|/p) * (1+tau)+1 (maximum no. of nodes when we have a tree. No. of edges in a tree=No. of Nodes -1)
    M = rhs_1 + 1
    # Connectivity Constraint # 2: sum j|(j,k) fjki<=(M-1)wik for i,k e V
    sum_x = []
    coeff = []
    f_variables = []
    for i in node:
        for k in node:
            sum_x.append(w_v[(i - 1) * no_nodes + (k - 1)])
            coeff.append(-M + 1)
            arcs_to = arcsto[k]
            indicies = [(q - 1) * no_nodes + (i - 1) for q in arcs_to]
            f_variables=[f_v_n[l] for j,l in enumerate(indicies)]
            sum_x.extend(f_variables)
            coeff1 = [1 for r in range(len(arcs_to))]
            coeff.extend(coeff1)
            c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["L"],
                                     rhs=[0])
            sum_x = []
            coeff = []
            f_variables = []

    # Connectivity Constraint # 3: sum j|(k,j) fkji -sum j|(j,k) fjki>= wik- Msik
    f_variables_1 = []
    f_variables_2 = []
    sum_x = []
    coeff = []
    for i in node:
        for k in node:
            sum_x.append(w_v[(i - 1) * no_nodes + (k - 1)])
            coeff.append(-1)
            s_variables = s_v[(i - 1) * no_nodes + (k - 1)]
            sum_x.append(s_variables)
            coeff.append(M)
            arcs_to = arcsto[k]
            arcs_from = arcsfrom[k]
            indicies_1 = [(q - 1) * no_nodes + (i - 1) for q in arcs_from]
            indicies_2 = [(q - 1) * no_nodes + (i - 1) for q in arcs_to]
            f_variables_1 = [f_v_n[l] for j, l in enumerate(indicies_1)]
            sum_x.extend(f_variables_1)
            f_variables_2 = [f_v_n[l] for j, l in enumerate(indicies_2)]
            sum_x.extend(f_variables_2)
            coeff1 = [1 for r in range(len(arcs_to))]
            coeff2 = [-1 for r in range(len(arcs_from))]
            coeff.extend(coeff1)
            coeff.extend(coeff2)
            c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["G"],
                                     rhs=[0])

            sum_x = []
            coeff = []
            f_variables_1 = []
            f_variables_2 = []

    # f_kji <= (M-1) xi,(j,k) for all i,j,k
    sum_x = []
    coeff = []
    for i, j in enumerate(x_v):
        sum_x.append(f_v_n[(edge_e[i] - 1) * no_nodes + node_i[i] - 1])
        sum_x.append(x_v[i])
        coeff.append(1)
        coeff.append(-M + 1)
        c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["L"],
                                 rhs=[0])
        sum_x = []
        coeff = []
    for i, j in enumerate(x_v):
        sum_x.append(f_v_n[(no_nodes * no_edges) + (edge_e[i] - 1) * no_nodes + node_i[i] - 1])
        sum_x.append(x_v[i])
        coeff.append(1)
        coeff.append(-M + 1)
        c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["L"],
                                 rhs=[0])
        sum_x = []
        coeff = []


class Model: # This is the model where we add the variables and the constraints
    def __init__(self,d_nodetoedge,x_v,w_v,s_v,f_v_n,edges_v,no_nodes,no_edges,num,t,edges,node):
        c = cplex.Cplex()
        # Setting the objective function to be Minmization
        c.objective.set_sense(c.objective.sense.minimize)
        # Declare decision variables (first argument is decision variables names, second argument is type of decision variables,
        # third argument is objective function coefficients)
        c.variables.add(names=x_v, types=["B"] * len(x_v), obj=d_nodetoedge)
        c.variables.add(names=w_v, types=["B"] * len(w_v))
        c.variables.add(names=s_v, types=["B"] * len(s_v))
        c.variables.add(names=f_v_n, types=["C"] * len(f_v_n))
        constraints_wo_cuts(c,edges_v,node,edges, no_nodes,no_edges,num,t,x_v,w_v,s_v,f_v_n)
        self.c=c
    def solve_districting(self):
        c=self.c
        t0 = time.time()
        c.solve()
        l_1 = round(time.time() - t0, 2)
        sol_1 = [x_v[i] for i, j in enumerate(c.solution.get_values(x_v)) if j > 0]
        sol_2 = [w_v[i] for i, j in enumerate(c.solution.get_values(w_v)) if j > 0]
        sol_3=[s_v[i] for i, j in enumerate(c.solution.get_values(s_v)) if j > 0]
        sol_4 = [f_v_n[i] for i, j in enumerate(c.solution.get_values(f_v_n)) if j > 0]
        obj = c.solution.get_objective_value()
        print("gap tolerance = ", c.parameters.mip.tolerances.mipgap.get())
        print(sol_1)
        print(sol_2)
        print(sol_3)
        print(sol_4)
        print(obj)
        return obj, l_1, c

probs=[78]
for prob in probs:
    file='data_'+str(prob)+'.csv'
    data1=Data(file)
    node=data1.node
    edges_v=data1.edges_v
    no_nodes=data1.no_nodes
    no_edges=data1.no_edges
    from_node_i=data1.from_node_i
    to_node_j=data1.to_node_j

    # creating a vector of index for every pair of node (i) and edge (e)
    node_i = []
    edge_e = []
    for i in node:
        for j in edges_v:
            node_i.append(i)
            edge_e.append(j)

    # Creating a vector for every pair of node (i) and node (j)
    node_i_n = []
    node_j = []
    for i in node:
        for j in node:
            node_i_n.append(i)
            node_j.append(j)
    # getting the flow
    from_node_a = copy.copy(from_node_i)
    to_node_a = copy.copy(to_node_j)
    for i in to_node_j:
        from_node_a.append(i)

    for i in from_node_i:
        to_node_a.append(i)

    # Creating a list of variable names (x_i,(j,k)): binary variable whether the edge (j,k) is assigned to district i
    x_v = []
    for i in range(len(node_i)):
        x = 'x' + str(node_i[i]) + '_' + str(edge_e[i])
        x_v.append(x)

    # Creating a list of variable names (wij): binary variable whether the node i is the center or not
    w_v = []
    for i in range(len(node_i_n)):
        w = 'w' + str(node_i_n[i]) + '_' + str(node_j[i])
        w_v.append(w)

    # Creating a list of variable names (sik): binary variable whether the node k is the sink for the center node i
    s_v = []
    for i in range(len(node_i_n)):
        s = 's' + str(node_i_n[i]) + '_' + str(node_j[i])
        s_v.append(s)

    # Creating a list of variable names (fkj): amount of flow from node k to node j
    f_v = []
    for i in range(len(from_node_a)):
        f = 'f' + str(from_node_a[i] + 1) + '_' + str(to_node_a[i] + 1)
        f_v.append(f)

    # Creating a list of variable names (fkji): amount of flow from node k to node j for district with center node i
    f_v_n = []

    for i in range(len(f_v)):
        for k in node:
            f_v_n.append((f_v[i] + '_' + str(node[k - 1])))

    # importing the distance from every node to every edge
    d_nodetoedge=np.loadtxt(open("nodetoedge_djalg_"+str(prob)+".csv", "rb"), delimiter=",", skiprows=1,max_rows=len(node_i),usecols=2).tolist()

    # Finding the arcs from each node:
    arcsfrom = collections.defaultdict(list)
    for i in range(len(from_node_a)):
        arcsfrom[from_node_a[i] + 1].append(i + 1)

    # Finding the arcs to each node:
    arcsto = collections.defaultdict(list)
    for i in range(len(to_node_a)):
        arcsto[to_node_a[i] + 1].append(i + 1)
    # Finding Neighboring nodes for each node (using only one loop)
    adj=collections.defaultdict(list)
    for i in range(len(from_node_i)):
        adj[from_node_i[i]+1].append(to_node_j[i]+1)
        adj[to_node_j[i]+1].append(from_node_i[i]+1)

    # Finding incident edges for every node
    edges=collections.defaultdict(list)
    for i in range(len(from_node_i)):
        edges[from_node_i[i]+1].append(i+1)
        edges[to_node_j[i]+1].append(i+1)

    # Creating a list of the degree for every node
    degree=[]
    for i in node:
        degree.append(len(edges[i]))

    for num in no_dist:
        for t in tol:
            print('Parameters')
            print(f'no of districts{num}')
            print(f'Tolerance{t}')
            distr_v = list(range(1, num + 1))  # Districts vector
            with open('Results/Exact_no_dist_'+str(num)+'_tol_'+str(t)+'_prob_'+str(prob)+'.csv','w') as newFile:
                newFileWriter = csv.writer(newFile, lineterminator='\n')
                newFileWriter.writerow(['Computation Time_WO_Logic_Cuts', 'Objective Function'])
                model=Model(d_nodetoedge,x_v,w_v,s_v,f_v_n,edges_v,no_nodes,no_edges,num,t,edges,node)
                obj,l_1,c=model.solve_districting()
                lower_bound_wo_cuts = c.solution.MIP.get_best_objective()
                newFileWriter.writerow([l_1 , round(obj, 2)])
                newFileWriter.writerow(['Lower_Bound'])
                newFileWriter.writerow([round(lower_bound_wo_cuts,2)])



