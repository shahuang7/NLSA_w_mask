# Perform initialization &
#   add cxfel code to system path.
import os

cxfel_root = os.environ['CXFEL_ROOT']
startup_file = cxfel_root+'/misc_tools/startup.py'
exec(open(startup_file).read())

# Read demo data, & set algorithmic parameters.
from misc_tools import read_h5

xyz = read_h5(cxfel_root+'/data/DM_Filled_Circle.mat','Y',h5=False,transpose=True)

yRow_yCol_yVal_file = 'sqDist.h5'
nN = 100
sigma_factor = 2.0
nEigs = 20
alpha = 1.0

# Calculate pairwise squared Euclidean distances using scipy.spatial.distance.cdist.
from scipy.spatial.distance import cdist

sqDist = cdist(xyz,xyz,'sqeuclidean')

# Truncate, symmetrize, and write out squared distances.
from sparsify import smallest_items_in_each_row_of_table,symmetrize
from misc_tools import write_h5

yRow,yCol,yVal = smallest_items_in_each_row_of_table(sqDist,nN)
yRow_symm,yCol_symm,yVal_symm = symmetrize(yRow,yCol,yVal)
write_h5(yRow_yCol_yVal_file,yRow_symm,'yRow')
write_h5(yRow_yCol_yVal_file,yCol_symm,'yCol')
write_h5(yRow_yCol_yVal_file,yVal_symm,'yVal')

# Calculate characteristic length scale of the data.
from ferguson import analyze

sigma_opt,_ = analyze(yRow_yCol_yVal_file)

# Calculate Diffusion Map embedding.
from diffmap import analyze,plot2D

sigma = sigma_factor*sigma_opt
h5_eigVec_eigVal = analyze(yRow_yCol_yVal_file,sigma,nEigs,alpha)

# Visualize the manifold.
eigVec = read_h5(h5_eigVec_eigVal,'eigVec')
for j in range(1,9):
  # colored based on \psi_j
  write_h5('colorcode.h5',eigVec[:,j]/eigVec[:,0],'colorcode')
  figure_name = plot2D(h5_eigVec_eigVal,[1,2],s=20)
  os.rename(figure_name,'diffmap_2D_psi_{}_colored.jpg'.format(j))
