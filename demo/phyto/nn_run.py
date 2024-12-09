import h5py
from scipy.sparse import csc_matrix
import time
import os
import numpy as np
from math import ceil
cxfel_root=os.environ['CXFEL_ROOT']
startup_file = cxfel_root+'misc_tools/startup.py'
exec(open(startup_file).read())

rawfile = '/home/uwm/huang229/Data/huang229/phyto/dataPhyto_upto1ps_int_sortdelay_unifdelay_LPF_DRL_SCL_BST_nS68084_nBrg62530.mat'
t0 = time.time()
f = open('phyto_nn_time.log','w')
jc = h5py.File(rawfile)['T_lpf_drl_scl_bst']['jc']
ir = h5py.File(rawfile)['T_lpf_drl_scl_bst']['ir']
data = h5py.File(rawfile)['T_lpf_drl_scl_bst']['data']
T = csc_matrix((data,ir,jc),shape=(68084,62530)).toarray()
jc = h5py.File(rawfile)['M_drl']['jc']
ir = h5py.File(rawfile)['M_drl']['ir']
data = h5py.File(rawfile)['M_drl']['data']
M = csc_matrix((data,ir,jc),shape=(68084,62530)).toarray()

N,D = T.shape
nN = 30000
c = 16384
sigma_factor = 10
nEigs = 5

num_worker = 32
n = ceil(N/num_worker)
alpha= 1.0
t1 = time.time()

from misc_tools import read_h5,write_h5
data_file = 'data_file_for_sna.h5'
variable_name = 'NxD_matrix_for_sna'
sqDist_file = 'sqDist.h5'
sq_code = cxfel_root+"misc_tools/prepare_squared_distance_file_.py"
write_h5(data_file,T,variable_name)
write_h5(data_file,M,'mask')
print('Preparation took {0:.2f} seconds'.format(t1-t0),file=f)

import subprocess 
subprocess.run(["mpiexec","-N",str(num_worker),"python",sq_code,data_file,variable_name,str(N),str(D),'dSqM',str(c),"True","False",str(n),str(nN),sqDist_file,"True","True","True","False"])

yVal = read_h5(sqDist_file,'yVal')
yCol = read_h5(sqDist_file,'yCol')
yRow = read_h5(sqDist_file,'yRow')
yVal[np.where(yCol==yRow)]=0
zero_diag_sq = 'newsq.h5'
write_h5(zero_diag_sq,yVal,'yVal')
write_h5(zero_diag_sq,yCol,'yCol')
write_h5(zero_diag_sq,yRow,'yRow')
t2 = time.time()
print("sqDist Done in {0:.2f} seconds".format(t2-t1),file=f)

ferguson_code = cxfel_root+"/ferguson/run_ferguson_.py"
subprocess.run(["python",ferguson_code,zero_diag_sq])
t3 = time.time()
print("Ferguson Analysis Done in {0:.2f} seconds".format(t3-t2),file=f)

diffmap_code = cxfel_root+"/diffmap/run_diffmap_.py"
subprocess.run(["python",diffmap_code,zero_diag_sq,str(sigma_factor),str(nEigs),str(alpha)])
eigVal = read_h5("eigVec_eigVal.h5",'eigVal')
eigVec = read_h5('eigVec_eigVal.h5','eigVec')
mu = (eigVec[:,0])*(eigVec[:,0])
psi = eigVec[:,1:].T/eigVec[:,0]
mu_psi_file = 'mu_psi.h5'
write_h5(mu_psi_file,mu,'mu')
write_h5(mu_psi_file,psi,'psi')
t4 = time.time()
print("Diffusion Map Done in {0:.2f} seconds".format(t4-t3),file=f)

dotp_code = cxfel_root+"/sna/run_sna_.py"
dotp_file = 'xtx.h5'
subprocess.run(["mpiexec","-N",str(num_worker),"python",dotp_code,data_file,variable_name,str(N),str(D),'dotM',str(c),"True","False",str(n),str(nN),dotp_file,"True","True","True"])
t5 = time.time()
print("Dot Product Done in {0:.2f} seconds".format(t5-t4),file=f)

cleanup_code = cxfel_root+"/sna/post_sna_cleanup_.py"
subprocess.run(["python",cleanup_code,'data_chunk',data_file,str(n),str(c)])
if c>1:
    subprocess.run(["python",cleanup_code,'pipe','dummy',str(n),str(c)])
t6 = time.time()
print("Post SnA Cleanup Done in {0:.2f} seconds".format(t6-t5),file=f)

ell = nEigs
num_copy = 2
nlsa_code = cxfel_root+"/nlsa/nlsa_v.py"
subprocess.run(["python",nlsa_code,data_file,variable_name,mu_psi_file,str(ell),str(N),str(D),str(n),str(c),str(num_copy)])
subprocess.run(["python",cleanup_code,'square','dummy',str(n),str(c)])
t7=time.time()
print("NLSA Done in {0:.2f} seconds".format(t7-t6),file=f)
f.close()