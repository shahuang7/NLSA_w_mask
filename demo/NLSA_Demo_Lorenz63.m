# Perform initialization &
#   add cxfel code to system path.
cxfel_root = getenv('CXFEL_ROOT')
addpath([cxfel_root '/wrappers'])
addpath([cxfel_root '/misc_tools'])
addpath([cxfel_root '/nlsa'])

##### System Specification #####
load([cxfel_root '/data/Lorenz63_x.mat'],'x_fine')
id = @(t) [floor(t*1000)+1];
testFunc = @(t) [x_fine(id(t))];
D = 1;

##### Simulation Specification #####
N = 11600;
delta_t = 50e-3; % fs

arg_list = argv();
command = sprintf('sigma_jitter = %s',arg_list{1});
eval(command)

##### samples #####
# k and jitter are in unit of timesteps:
k = ([1:N]'-1);
jitter = randn(N,1)*sigma_jitter/delta_t;
 
# k and jitter are in the normalized unit used in the differential equations:
# interested in 20 < t < 24
k = k/N*4+20;
jitter = jitter/N*4;

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
n=500,cleanup='True',no_block='True',run_mpi='False',num_worker=1
run_prepare_squared_distance_file_py(data_file,variable_name,N,D,'dSq',c,...
h5,transpose,n,nN,yRow_yCol_yVal_file,cleanup,no_block,run_mpi,num_worker)



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

cleanup='False',no_block='True',run_mpi='False',num_worker=1
run_sna_py(data_file,variable_name,N,D,'dot',c,...
h5,transpose,n,nN,yRow_yCol_yVal_file,cleanup,no_block,run_mpi,num_worker)

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

U_NLSA = real(U_NLSA);
S_NLSA = real(diag(S_NLSA));
V_NLSA = real(V_NLSA);

##### region of interest (both for display and for fitting) #####
x_time = [0:N-c+1-num_copy];
ROI = (c+num_copy)/2+x_time;
T_jittered = T_jittered(ROI,:);
T_jitter_free = T_jitter_free(ROI,:);

# reconstruction & fit to known

num_mode = nEigs;
X = reconstruct(D,[1:num_mode],x_time,[1:num_copy],U_NLSA,eye(num_mode),V_NLSA);
Z = permute(X,[2,3,1]);
W = ones(1,length(x_time));
[best_fit_NLSA,coeff]=Component_Fit(Z,T_jitter_free',diag(W));

err = best_fit_NLSA-T_jitter_free';
err = err(:);
err = sqrt(mean(err.*err));

##### graphics-related #####

hFigure = figure(1);
set(hFigure,'color','w')
set(hFigure,'resize','off')
Pix_SS = get(0,'screensize');
screenWidth = Pix_SS(3);
screenHeight = Pix_SS(4);
pos = [10 900 screenWidth/2 0.8*screenHeight];
try
  warning('off','Octave:abbreviated-property-match')
catch
end
set(hFigure,'pos',pos)

figure(hFigure)
 
hsp = subplot(2,2,1);
x = x_time;             my_xlabel = 'Time (\Deltat)';
y = T_jitter_free(:,1); my_ylabel = '';
y_min = min(y); y_max = max(y);
my_title = 'Jitter-free Signal';
plotRF(hsp,x,y,my_xlabel,my_ylabel,my_title,'b-')
set(hsp,'ylim',[y_min,y_max]+[-1,1]*0.1*(y_max-y_min))

hsp = subplot(2,2,2);
x = x_time;             my_xlabel = 'Time (\Deltat)';
y = T_jittered(:,1);    my_ylabel = '';
y_min = min(y); y_max = max(y);
my_title = 'Input Signal to NLSA';
plotRF(hsp,x,y,my_xlabel,my_ylabel,my_title,'b.')
hLine = get(hsp,'children');
set(hLine,'markerSize',4)
set(hsp,'ylim',[y_min,y_max]+[-1,1]*0.1*(y_max-y_min))

hsp = subplot(2,2,3);
num_S = min([nEigs c*D]);
x = [1:num_S];          my_xlabel = 'Singular Value #';
y = S_NLSA(1:num_S);    my_ylabel = 'Magnitude (A.U.)';
y = y/y(1);
y_min = 0.0; y_max = 1.0;
my_title = 'Singular Value Spectrum';
plotRF(hsp,x,y,my_xlabel,my_ylabel,my_title,'b-o')
set(hsp,'xlim',[0.8 num_S],'ylim',[y_min,y_max],'xTick',[5:5:num_S])

hsp = subplot(2,2,4);
x = x_time;             my_xlabel = 'Time (\Deltat)';
y = T_jitter_free(:,1); my_ylabel = '';
y_min = min(y); y_max = max(y);
my_title = 'Reconstructed from NLSA Output';
plotRF(hsp,x,y,my_xlabel,my_ylabel,my_title,'b-')
y = best_fit_NLSA(1,:);
addplotRF(hsp,x,y,'r-')
set(hsp,'ylim',[y_min,y_max]+[-1,1]*0.1*(y_max-y_min))
legend({'Jitter-free Signal','Reconstruction'},'location','northeast');

myJPEG = sprintf('NLSA_Demo_Lorenz63_sigma_jitter_%.2f.jpg',sigma_jitter);
print(myJPEG)
