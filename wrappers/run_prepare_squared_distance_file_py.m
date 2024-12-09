function run_prepare_squared_distance_file_py(data_file,variable_name,N,D,measure_type,c,h5,transpose,n,nN,sqDist_file,cleanup,no_block,run_mpi,num_worker)
  
  library_root = getenv('CXFEL_ROOT');
  library_path = [library_root '/misc_tools'];
  
  command = sprintf(['python3 %s/prepare_squared_distance_file_.py ',...
    '%s %s %d %d %s %d %s %s %d %d %s %s %s %s'],...
    library_path,...
    data_file,variable_name,N,D,measure_type,c,h5,transpose,n,nN,sqDist_file,...
    cleanup,no_block,run_mpi);
  if strcmp(run_mpi,'True')
    mpi_prefix = sprintf('mpiexec -N %d ',num_worker);
    command = [mpi_prefix command];
  end
  
  system_with_anaconda3_path(command);
