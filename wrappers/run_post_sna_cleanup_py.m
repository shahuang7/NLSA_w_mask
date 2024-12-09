function run_post_sna_cleanup_py(block_type,data_file,n,c)
  
  library_root = getenv('CXFEL_ROOT');
  library_path = [library_root '/sna'];
  command = sprintf('python3 %s/post_sna_cleanup_.py %s %s %d %d',...
              library_path,block_type,data_file,n,c);
  system_with_anaconda3_path(command);
