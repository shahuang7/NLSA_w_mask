function run_ferguson_py(yRow_yCol_yVal_file)
  
  library_root = getenv('CXFEL_ROOT');
  library_path = [library_root '/ferguson'];
  
  command = sprintf('python3 %s/run_ferguson_.py %s',library_path,yRow_yCol_yVal_file);
  system_with_anaconda3_path(command);

