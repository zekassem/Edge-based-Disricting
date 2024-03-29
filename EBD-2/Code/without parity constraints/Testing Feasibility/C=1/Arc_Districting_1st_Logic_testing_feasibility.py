import cplex
import time
import numpy as np
import collections
import csv

# Loading Data
class Data():
    def __init__(self,file):
        self.file=file
        self.no_nodes=np.loadtxt(file,skiprows=1,max_rows=1,usecols=0,dtype=int) # skipping lines to get the number of nodes in the planar graph
        self.no_nodes1=int(self.no_nodes)
        self.coordinates= np.loadtxt(file,skiprows=4,max_rows=self.no_nodes1,usecols=(1,2), delimiter=',') # skipping lines to get the x & y coordinates of the nodes
        self.s=self.no_nodes+4+2 # Calculating the lines needed to be skipped to get the number of edges
        self.no_edges=np.loadtxt(file,skiprows=self.s,max_rows=1,usecols=0,dtype=int) # skipping lines to get the
        self.no_edges1=int(self.no_edges)
        self.z=self.s+4 # Calculating the lines needed to be skipped to get adjacency list
        self.adj_list=np.loadtxt(file,skiprows=self.z,max_rows=self.no_edges1,usecols=(0,1),dtype=int, delimiter=',')
        self.node=list(range(1,(self.no_nodes)+1))
        self.from_node_i=self.adj_list[:,0].tolist()
        self.to_node_j=self.adj_list[:,1].tolist()
        self.edges_v=list(range(1,self.no_edges+1)) #Edges Vector

## End of Data
# Input parameters to the districting problem
no_dist=[2,4,6,8,10] # the number of districts
tol=[0.02,0.04,0.06]# tolerance parameter into the balancing constraints

# Breadth First Function
def BFS(edges_v,index, result, l):
    # defining the arrays that will be in the algorithm
    pred = []  # predecessor vector
    color = []  # color vector, 0 means white, 1 means gray ,2 means black
    d = []  # distance vector
    Q = []  # set of gray vectors
    s = l  # the source edge index in the edge vector
    for e in edges_v:
        color.append(int(0))  # having all the colors of edges to be white
        d.append(int(0))
        pred.append(int(0))

    color[s - 1] = 1
    Q.append(s)
    current_dis = 0

    while len(Q) != 0:  # while the cardinality of set of gray edges is not equal to zero
        u = Q.pop(0)  # Dequeue the first edge
        edges_nei=edges_neighboring[u] # Neighboring edges

        edges_neigh_n = list(set(edges_nei).intersection(set(result[index]))) # We're only considering the edges that are selected in the territory that's why we're doing the intersection
        for i in edges_neigh_n:  # This is the main loop
            if color[i - 1] == 0:
                color[i - 1] = 1
                d[i - 1] = current_dis + 1
                pred[i - 1] = u
                Q.append(i)  # appending the set of gray nodes

        color[u - 1] = 2

    b = color
    return b

def constraints_wo_cuts(c,edges_v, no_nodes,no_edges,no_dist,t,x_v,w_v): #  This function adds all of the constraints for the original districting problem
    # input edges_v, no_nodes,no_edges,no_dist,t
    # sum i e V (xie)=1 for e e E: each edge is assigned to one territory
    for e in edges_v:
        sum_x = x_v[(e - 1):no_edges * no_nodes:no_edges]
        coeff = [1 for i in range(no_nodes)]
        c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["E"], rhs=[1])
        sum_x = []
        coeff = []

    # sum i e V w_i=p
    sum_x = w_v
    coeff = [1 for i in range(no_nodes)]
    c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["E"],
                             rhs=[no_dist])
    sum_x = []
    coeff = []
    # Balancing Constraints
    # sum e e E xie <= (|E|/p) * (1+tau) wi for i e V
    rhs_1 = (no_edges / no_dist) * (1 + t)
    for i in node:
        sum_x.append(w_v[i - 1])
        coeff.append(-rhs_1)
        sum_x.extend(x_v[(i - 1) * no_edges:i * no_edges])
        coeff2 = [1 for q in range(no_edges)]
        coeff.extend(coeff2)
        c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["L"], rhs=[0])
        coeff = []
        sum_x = []

    # sum e e E xie >= (|E|/p) * (1-tau) wi for i e V
    sum_x = []
    coeff = []
    rhs_2 = (no_edges / no_dist) * (1 - t)
    for i in node:
        sum_x.append(w_v[i - 1])
        sum_x.extend(x_v[(i - 1) * no_edges:i * no_edges])
        coeff.append(-rhs_2)
        coeff2 = [1 for q in range(no_edges)]
        coeff.extend(coeff2)
        c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["G"], rhs=[0])

        coeff = []
        sum_x = []



def logic_cuts(c,w_v,edges,x_v,no_edges): # function that adds the logic cuts
    # Logic Cut sum (i,k) e delta(i) xi,(i,k) >= 2 * wi for i e V, logic cuts are added as user cuts
    sum_x = []
    coeff = []
    for i in node:
        sum_x.append(w_v[i - 1])
        coeff.append(-2)
        incident_edges = edges[i]
        indicies = [(q - 1) + (i - 1) * no_edges for q in incident_edges]
        x_variables = [x_v[i] for i in indicies]
        sum_x.extend(x_variables)
        coeff2 = [1 for r in range(len(incident_edges))]
        coeff.extend(coeff2)
        c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["G"],
                                 rhs=[0])
        sum_x = []
        coeff = []
        x_variables = []

class Model: # This is the model where we add the variables and the constraints
    def __init__(self,d_nodetoedge,x_v,w_v,edges_v,no_nodes,no_edges,num,t):
        c = cplex.Cplex()
        # Setting the objective function to be Minmization
        c.objective.set_sense(c.objective.sense.minimize)
        # Declare decision variables (first argument is decision variables names, second argument is type of decision variables,
        # third argument is objective function coefficients)
        c.variables.add(names=x_v, types=["B"] * len(x_v), obj=d_nodetoedge)
        c.variables.add(names=w_v, types=["B"] * len(w_v))
        constraints_wo_cuts(c, edges_v, no_nodes, no_edges, num, t,x_v,w_v)
        self.c=c
    def branch_cut(self,edges_v): # This is the method that has the branch and cut
        c=self.c
        c.parameters.clocktype.set(2)
        c.parameters.timelimit.set(43200)
        t0 = time.time()
        c.solve()
        result=[]
        center_node1 = [i for i, val in enumerate(c.solution.get_values(w_v)) if val > 0.01]
        center_node = [i + 1 for i in center_node1]
        for o in center_node:
            result.append([edge_e[i] for i,val in enumerate(c.solution.get_values(x_v)) if val>0.01 and node_i[i]==o])
        a = 1
        num_iterations=0 # counter for the number of iterations
        num_cuts=0 # counter for the number of
        while a > 0:
            num_iterations+=1
            C = []
            index = 0
            for o in center_node:  # This loop is to detect whether there are disconnected districts or not
                explored_edges_s1 = []  # edges that are found using BFS
                R = []  # Set of Edges that need to be explored using BFS
                R.extend(result[index])
                l = R[0]  # the source edge from which we will start exploring
                b = BFS(edges_v, index, result, l)
                explored_edges_s1.extend([i + 1 for i, val in enumerate(b) if val == 2])
                explored_edges_s = list(set(explored_edges_s1))
                unexplored_edges = list(
                    set(R).difference(set(explored_edges_s)))  # list of unexplored edges within a district
                if len(unexplored_edges) > 0:
                    C.append(0)
                else:
                    C.append(1)
                explored_edges_s1 = []  # list of explored edges for every district to keep track of all connected components
                while len(
                        unexplored_edges) > 0:  # This while loop to find all of the different connected components and add all of the needed cuts
                    explored_edges_s1.extend([i + 1 for i, val in enumerate(b) if val == 2])
                    explored_edges_s = list(set(explored_edges_s1))
                    unexplored_edges = list(
                        set(R).difference(set(explored_edges_s)))  # list of unexplored edges within a district
                    # Find the connected component
                    # Find the neighboring edges to the connected component (excluding edges in the connected component)
                    # Add the needed cuts
                    connect_edges = [j for i, j in enumerate(R) if b[j - 1] == 2]  # Connected Component (Set Sk)
                    connect_edges_neigh_nested = [edges_neighboring[i] for i in connect_edges]
                    connect_edges_neighboring1 = [j for sublist in connect_edges_neigh_nested for j in sublist]
                    connect_edges_neighboring = set(connect_edges_neighboring1)  # Removing duplicated edges
                    connect_edges_neighboring_n = list(connect_edges_neighboring.difference(
                        set(
                            connect_edges)))  # Neighboring edges to connected component excluding edges in the connected component
                    s = (1 - len(connect_edges))  # 1-|S|
                    sum_x = []
                    coeff = []
                    indicies_1 = [(o - 1) * no_edges + (q - 1) for q in connect_edges]
                    x_variables = [x_v[i] for i in indicies_1]
                    sum_x.extend(x_variables)
                    coeff1 = [-1 for r in range(len(indicies_1))]
                    coeff.extend(coeff1)
                    indicies_2 = [(o - 1) * no_edges + (q - 1) for q in connect_edges_neighboring_n]
                    x_variables_1 = [x_v[i] for i in indicies_2]
                    coeff2 = [1 for r in range(len(indicies_2))]
                    coeff.extend(coeff2)
                    sum_x.extend(x_variables_1)
                    if len(sum_x)>0:
                        num_cuts+=1
                    c.linear_constraints.add(lin_expr=[cplex.SparsePair(sum_x, coeff)], senses=["G"], rhs=[s])
                    if len(unexplored_edges) > 0:  # finding the next connected component
                        l = unexplored_edges[0]
                        b = BFS(edges_v, index, result, l)
                index = index + 1
            if sum(C) < len(center_node):
                a = 1
                c.solve()
                center_node1 = [i for i, val in enumerate(c.solution.get_values(w_v)) if val > 0.01]
                center_node = [i + 1 for i in center_node1]
                result = []
                for o in center_node:
                    result.append(
                        [edge_e[i] for i, val in enumerate(c.solution.get_values(x_v)) if val > 0.01 and node_i[i] == o])
            else:
                a = 0
        l_1 = round(time.time() - t0, 2)
        sol_1 = [x_v[i] for i, j in enumerate(c.solution.get_values(x_v)) if j > 0.01]
        sol_2 = [w_v[i] for i, j in enumerate(c.solution.get_values(w_v)) if j > 0.01]
        obj = c.solution.get_objective_value()
        print("gap tolerance = ", c.parameters.mip.tolerances.mipgap.get())
        print(sol_1)
        print(sol_2)
        print(obj)

        return obj,l_1,c,num_iterations,num_cuts

class Model_W_Cuts(Model): # creating a separate model with cuts (inheritance from Model)
    def __init__(self,edges):
        super().__init__(d_nodetoedge,x_v,w_v,edges_v,no_nodes,no_edges,num,t)
        c=self.c
        logic_cuts(c,w_v,edges,x_v,no_edges)
    def branch_cut_logic(self):
        super().branch_cut(edges_v) # will exactly the same

probs=[87,108,113,122,130,149,151,158]
for prob in probs:
    file='data_'+str(prob)+'.csv'
    data1=Data(file)
    node=data1.node
    edges_v=data1.edges_v
    no_nodes=data1.no_nodes
    no_edges=data1.no_edges
    # creating a vector of index for every pair of node (i) and edge (e)
    node_i = []
    edge_e = []
    for i in node:
        for j in edges_v:
            node_i.append(i)
            edge_e.append(j)


    # Creating Variables

    # Creating a list of variable names (x_i,(j,k)): binary variable whether the edge (j,k) is assigned to district i
    x_v = []
    for i in range(len(node_i)):
        x = 'x' + str(node_i[i]) + '_' + str(edge_e[i])
        x_v.append(x)


    # Creating a list of variable names (wi): binary variable whether the node i is the center or not
    w_v = []
    for i in range(len(node)):
        w = 'w' + str(node[i])
        w_v.append(w)

    # creating a list of variable names (zi): binary variable whether there are adjacent center in the neighborhood of i or not
    z_v = []
    for i in range(len(node)):
        z = 'z' + str(node[i])
        z_v.append(z)
    # importing the distance from every node to every edge
    d_nodetoedge=np.loadtxt(open("nodetoedge_djalg_"+str(prob)+".csv", "rb"), delimiter=",", skiprows=1,max_rows=len(node_i),usecols=2).tolist()

    from_node_i=data1.from_node_i
    to_node_j=data1.to_node_j
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


    # Creating a dictionary for neighboring edges for each edge
    edges_neighboring= collections.defaultdict(list)
    for e in edges_v:
        nodes_corr=[]
        nodes_corr.append(from_node_i[e - 1]+1)
        nodes_corr.append(to_node_j[e - 1]+1)
        edges_nei=[]
        for i in nodes_corr:
            edges_nei.extend(edges[i])
        edges_nei1=list(set(edges_nei))
        l=[]
        l.append(e)
        edges_nei=[x for x in edges_nei1 if x not in l]
        edges_neighboring[e].extend(edges_nei)
    for num in no_dist:
        for t in tol:
            print('Parameters')
            print(f'no of districts{num}')
            print(f'Balancing_Tolerance{t}')
            distr_v = list(range(1, num + 1))  # Districts vector
            try:
                with open('Results/1st_Logic_no_dist_'+str(num)+'_tol_'+str(t)+'_prob_'+str(prob)+'.csv','w') as newFile:
                    newFileWriter = csv.writer(newFile, lineterminator='\n')
                    newFileWriter.writerow(['Computation Time_WO_Logic_Cuts', 'Objective Function'])
                    model=Model(d_nodetoedge,x_v,w_v,edges_v,no_nodes,no_edges,num,t)
                    obj,l_1,c,num_iterations_wo_cuts,num_cuts_wo_logic=model.branch_cut(edges_v)
                    lower_bound_wo_cuts=c.solution.MIP.get_best_objective()
                    newFileWriter.writerow([l_1 , round(obj, 2)])
                    newFileWriter.writerow(['Lower_Bound', 'Number of Iterations','Number of Cuts'])
                    newFileWriter.writerow([round(lower_bound_wo_cuts, 2), num_iterations_wo_cuts,num_cuts_wo_logic])
                    model_w_cuts=Model_W_Cuts(edges)
                    obj,l_1,c,num_iterations_w_cuts,num_cuts_w_logic=model_w_cuts.branch_cut(edges_v)
                    lower_bound_w_cuts = c.solution.MIP.get_best_objective()
                    newFileWriter.writerow(['Computation Time_W_Logic_Cuts', 'Objective Function'])
                    newFileWriter.writerow([l_1, round(obj, 2)])
                    newFileWriter.writerow(['Lower_Bound','Number of Iterations','Number of Cuts'])
                    newFileWriter.writerow([round(lower_bound_w_cuts,2),num_iterations_w_cuts,num_cuts_w_logic])
            except:
                with open('Results/1st_Logic_no_dist_'+str(num)+'_tol_'+str(t)+'_prob_'+str(prob)+'.csv','w') as newFile:
                    newFileWriter = csv.writer(newFile, lineterminator='\n')
                    newFileWriter.writerow(['Computation Time_WO_Logic_Cuts', 'Objective Function'])
                    newFileWriter.writerow(['Infeasible'])
                pass


