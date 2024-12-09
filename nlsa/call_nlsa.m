cxfel_root = getenv('CXFEL_ROOT')
addpath([cxfel_root '/wrappers'])
addpath([cxfel_root '/misc_tools'])
addpath([cxfel_root '/nlsa'])
arg_list = argv();
command = sprintf('ell = %s',arg_list{1}); eval(command)
command = sprintf('D = %s',arg_list{2}); eval(command)
command = sprintf('N = %s',arg_list{3}); eval(command)
command = sprintf('n = %s',arg_list{4}); eval(command)
command = sprintf('c = %s',arg_list{5}); eval(command)
command = sprintf('num_copy = %s',arg_list{6}); eval(command)
%command = sprintf('norm = %s',arg_list{7}); eval(command)

data_file = 'data_file_for_sna.h5';
variable_name = 'NxD_matrix_for_sna';

load(data_file,variable_name);
eval(sprintf('X1 = %s'';',variable_name))
X1 = X1';

load('mu_psi.h5','mu');
mu = mu';
load('mu_psi.h5','psi');

[U_NLSA,S_NLSA,V_NLSA] = extract_topos_chronos(ell,X1,mu,psi,D,N,n,c,num_copy);
run_post_sna_cleanup_py('square','dummy',n,c)
save('uv.mat','-v7','U_NLSA','S_NLSA','V_NLSA')