# Perform initialization &
#   add cxfel code to system path.
cxfel_root = getenv('CXFEL_ROOT')
addpath([cxfel_root '/wrappers'])
addpath([cxfel_root '/misc_tools'])
addpath([cxfel_root '/nlsa'])

##### System Specification #####
load([cxfel_root '/data/Eduardo_cryoEM_N100001_D10.mat'],'y')
id = @(t) [floor(t*1000)+1];
testFunc = @(t) y(id(t),1:10);
D = 10;

##### Simulation Specification #####
N = 11600;
delta_t = 50e-3; % fs

arg_list = argv();
command = sprintf('sigma_jitter = %s',arg_list{1}); eval(command)
command = sprintf('num_worker = %s',arg_list{2}); eval(command)
command = sprintf('no_block = ''%s''',arg_list{3}); eval(command)
command = sprintf('run_mpi = ''%s''',arg_list{4}); eval(command)

if (num_worker>1), run_mpi = 'True'; end

##### samples #####
# k and jitter are in unit of timesteps:
k = ([1:N]'-1);
jitter = randn(N,1)*sigma_jitter/delta_t;
 
sigma_jitter_max = 10;
jitter_max = 6*sigma_jitter_max/delta_t;
jitter_range = 2*jitter_max;
a = 100/(N-1+jitter_range); b = a*jitter_max;
k = a*k+b;                  jitter = a*jitter;

T_jitter_free = testFunc(k);
T_jittered = testFunc(k+jitter);

T_jittered_saved = T_jittered;

save('T_jittered_saved.mat','T_jittered_saved','-v7')

# Read demo data, & set algorithmic parameters.
data_file='T_jittered_saved.mat';
variable_name='T_jittered_saved',h5='False',transpose='False'


yRow_yCol_yVal_file = 'sqDist.h5'
nN = 1000
sigma_factor = 2.0
nEigs = 30
alpha = 1.0
c = 5800

# Calculate symmetrized, truncated, pairwise squared Euclidean distances using
# the Shift-and-Add library.
n=500,cleanup='True'

run_prepare_squared_distance_file_py(data_file,variable_name,N,D,'dSq',c,h5,transpose,n,nN,yRow_yCol_yVal_file,cleanup,no_block,run_mpi,num_worker)

# Calculate characteristic length scale of the data.
run_ferguson_py(yRow_yCol_yVal_file)
load('sigma_opt.h5','sigma_opt')


# Calculate Diffusion Map embedding.
run_diffmap_py(yRow_yCol_yVal_file,sigma_factor,nEigs,alpha)

load('eigVec_eigVal.h5','eigVec')
eigVec = eigVec';
mu = eigVec(:,1).^2;
psi = bsxfun(@rdivide,eigVec(:,2:end),eigVec(:,1));

# Calculate pairwise dot products using
# the Shift-and-Add library.

run_sna_py(data_file,variable_name,N,D,'dot',c,h5,transpose,n,nN,yRow_yCol_yVal_file,cleanup,no_block,run_mpi,num_worker)
delete('.sna_run_info')
[filepath,name,ext] = fileparts(data_file);
run_post_sna_cleanup_py('data_chunk',[name ext],n,c)
if (c>1), run_post_sna_cleanup_py('pipe','dummy',n,c), end

# nlsa.
ell = nEigs;
num_copy = 2;

load(data_file,variable_name)
eval(sprintf('X1 = %s'';',variable_name))

[U_NLSA,S_NLSA,V_NLSA] = extract_topos_chronos(ell,X1,mu,psi,D,N,n,c,num_copy);

run_post_sna_cleanup_py('square','dummy',n,c)

calc_name_template = 'NLSA_cryoEM_spikeAI_sigma_jitter_%.2f';
myMAT  = sprintf([calc_name_template '.mat'],sigma_jitter);
myJPEG = sprintf([calc_name_template '.jpg'],sigma_jitter);
save(myMAT,'c','D','myJPEG','N','nEigs','num_copy','T_jitter_free','T_jittered','U_NLSA','S_NLSA','V_NLSA')

if exist('.NLSA_only','file'), return, end

system(['octave NLSA_recon_comp_fit_multi_pixel.m ' myMAT]);
