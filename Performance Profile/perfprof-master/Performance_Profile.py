from perfprof import *
import numpy as np
import time
t0 = time.time()
alg_file = 'Algorithm_Time.csv'
Time_wo_Cut = np.loadtxt(alg_file, skiprows=1, max_rows=281, usecols=0,delimiter=',',dtype=float).tolist()
Time_1st_Cut = np.loadtxt(alg_file, skiprows=1, max_rows=281, usecols=1,delimiter=',',dtype=float).tolist()
Time_1st_Cut_w_z=np.loadtxt(alg_file, skiprows=1, max_rows=281, usecols=2,delimiter=',',dtype=float).tolist()
Time_1st_Cut_wo_z=np.loadtxt(alg_file, skiprows=1, max_rows=281, usecols=3,delimiter=',',dtype=float).tolist()


linespecs = ['r--', 'b-','g:','m-.']
labels = ['0', '1','1,2a','1,2b']
data = np.vstack((Time_wo_Cut,Time_1st_Cut, Time_1st_Cut_w_z,Time_1st_Cut_wo_z)).T
perfprof(data, linespecs=linespecs, legendnames=labels, usetex=True)
plt.savefig(r'fig/Newest_performance.png')
l_1 = round(time.time() - t0, 2)
print(l_1)
