function measure=read_square_block(n,row,col,c)
  
  template = '.square_block_n%d_row%d_col%d_c%d.h5';
  fileName = sprintf(template,n,row-1,col-1,c);
  fid = fopen(fileName,'r');
  if (fid<0)
    measure = zeros(n);
  else
    load(fileName,'measure')
    measure = measure';
  end
  fclose(fid);
