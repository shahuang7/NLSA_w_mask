# Perform initialization &
#   add cxfel code to system path.
cxfel_root = getenv('CXFEL_ROOT')
addpath([cxfel_root '/wrappers'])

# Read demo data, & set algorithmic parameters.
data_file='DM_Spiral.mat';variable_name='Y';N=2000;D=3;
h5='False';transpose='True';
system(['ln -s ' cxfel_root '/data/' data_file ' ./' data_file])

yRow_yCol_yVal_file = 'sqDist.h5'
nN = 100
sigma_factor = 2.0
nEigs = 20
alpha = 1.0

# Calculate symmetrized, truncated, pairwise squared Euclidean distances using
# the Shift-and-Add library.
c=1;n=500;cleanup='True';no_block='True';run_mpi='False';num_worker=1;
run_prepare_squared_distance_file_py(data_file,variable_name,N,D,'dSq',c,...
  h5,transpose,n,nN,yRow_yCol_yVal_file,cleanup,...
  no_block,run_mpi,num_worker)
system(['rm ' data_file])

# Calculate characteristic length scale of the data.
run_ferguson_py(yRow_yCol_yVal_file)
load('sigma_opt.h5','sigma_opt')

# Calculate Diffusion Map embedding.
sigma = sigma_factor*sigma_opt
run_diffmap_py(yRow_yCol_yVal_file,sigma_factor,nEigs,alpha)

# Visualize the manifold.
load('eigVec_eigVal.h5','eigVec')
eigVec = eigVec';
x = eigVec(:,2)./eigVec(:,1); my_xlabel = '\psi_1';
y = eigVec(:,3)./eigVec(:,1); my_ylabel = '\psi_2';
z = eigVec(:,4)./eigVec(:,1); my_zlabel = '\psi_3';
s = 50;
color = x;
h = figure('position',[300,200,1500,500]);
hsp = subplot(1,3,1); scatter(x,y,s,color), axis equal,
  xlabel(my_xlabel,'fontsize',15), ylabel(my_ylabel,'fontsize',15)
hsp = subplot(1,3,2); scatter(x,z,s,color), axis equal,
  xlabel(my_xlabel,'fontsize',15), ylabel(my_zlabel,'fontsize',15)
hsp = subplot(1,3,3); scatter(y,z,s,color), axis equal,
  xlabel(my_ylabel,'fontsize',15), ylabel(my_zlabel,'fontsize',15)
colormap hsv
print('diffmap_2D_psi_1_colored.jpg')
