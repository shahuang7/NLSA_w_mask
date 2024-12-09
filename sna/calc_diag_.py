################################################################################
def calc_diag(n,diag):
# 
# copyright (c) Russell Fung 2023
################################################################################
  
  from sna import diag_file_template,read_block,read_run_info,write_diag
  
  import os
  
  diag_file = diag_file_template().format(n,diag,1)
  if os.path.exists(diag_file): return
  
  import numpy as np
  
  _,_,N,_,_,_,_,_ = read_run_info()
  c = 1
  num_super = N-c+1
  num_block,leftover = np.divmod(num_super,n)
  num_block += (leftover>0)
  num_block -= diag
  measure = np.zeros((num_block*n,n))
  mask = np.zeros((num_block*n,n))
  
  for block in range(num_block):
    row = block
    col_triu = diag+block
    col_tril = col_triu+1
    
    triu_square_block, triu_mask = read_block('square',n,row,col_triu,1)
    tril_square_block, tril_mask = read_block('square',n,row,col_tril,1)
    
    sandbox = np.hstack((triu_square_block,tril_square_block))
    sandbox_mask = np.hstack((triu_mask,tril_mask))

    for jj in range(n):
      measure[block*n+jj,:] = sandbox[jj,jj:jj+n]
      mask[block*n+jj,:] = sandbox_mask[jj,jj:jj+n]
  measure = measure[:num_super,:]
  mask = mask[:num_super,:]
  
  write_diag(n,diag,1,measure,mask)

