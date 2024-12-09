################################################################################
def plot(log_sigma,log_sum_kernel,x_linear,y_linear,sigma_opt,dimensionality):
# 
# copyright (c) Laura Williams & Russell Fung 2018
# updated April 13, 2022
################################################################################
  
  import socket
  if ("compute-" in socket.gethostname()):
    return
  
  import matplotlib.pyplot as plt
  import numpy as np
  
  x = log_sigma
  y = log_sum_kernel
  y_min = np.min(y)
  y_max = np.max(y)
  
  x_sigma_opt = np.max(x_linear)
  y_sigma_opt = y_min
  x_dimensionality = x_sigma_opt
  y_dimensionality = 0.5*(y_min+y_max)
  
  plt.figure()
  plt.plot(x,y,'bo-',linewidth=2.0,fillstyle='none')
  plt.plot(x_linear,y_linear,'r-',linewidth=2.0)
  plt.xlabel('$\ln \sigma$',fontsize=15)
  plt.ylabel('$\ln \Sigma A_{ij}$',fontsize=15)
  plt.text(x_sigma_opt,y_sigma_opt,'$\sigma_{opt}=%6.4f$'%sigma_opt,fontsize=15)
  plt.text(x_dimensionality,y_dimensionality,'dimensionality=%6.2f'%dimensionality,fontsize=15)
  plt.show(block=False)
  figure_name = 'ferguson.jpg'
  plt.savefig(figure_name)
  plt.close()
  
  return figure_name

