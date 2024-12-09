function run_diffmap_py(yRow_yCol_yVal_file,sigma_factor,nEigs,alpha)
  
  library_root = getenv('CXFEL_ROOT');
  library_path = [library_root '/diffmap'];
  
  command = sprintf('python3 %s/run_diffmap_.py %s %.4f %d %.4f',library_path,yRow_yCol_yVal_file,...
    sigma_factor,nEigs,alpha);
  system_with_anaconda3_path(command);

