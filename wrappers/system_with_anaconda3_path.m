function system_with_anaconda3_path(command)
  
  command = ['export PATH="/opt/anaconda3/bin:$PATH" && ' command];
  system(command);
