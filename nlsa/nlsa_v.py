def reconstruct(D,k_interest,t_interest,copy_interest,U,V):
    import numpy as np
    from math import floor
    num_k = len(k_interest)
    num_t = len(t_interest)
    num_cc = len(copy_interest)
    X = np.zeros((num_k,D,num_t))
    for ik in range(num_k):
        Xk = np.zeros((D,num_t))
        print('Reconstructing Mode {}'.format(ik))
        for cc in copy_interest:
            Xk = Xk + np.outer(U[cc*D:(cc+1)*D,ik],V[cc:cc+num_t,ik])
        X[ik,:,:] = Xk/num_cc
    return X

def nlsa(data_file,variable_name,mu_psi_file,ell,N,D,n,c,num_copy):
    import time
    from misc_tools import read_h5,write_h5
    from sna import read_block
    from math import ceil
    
    t0 = time.time()
    X1 = read_h5(data_file,variable_name).T
    mu = read_h5(mu_psi_file,'mu')
    psi = read_h5(mu_psi_file,'psi')
    nS = N-c+1
    psi_1 = np.concatenate((np.ones((1,nS)),psi))
    psi = psi_1
    del psi_1
    mu_psi = mu*psi
    mu_psi = mu_psi.T
    

    ATA = np.zeros((ell+1,ell+1))
    numRow = ceil(nS/n)
    numCol = numRow
    numFiles = numRow*numCol
    fileNum = 0
    for row in range(numRow):
        j0 = row*n
        j1 = j0+n
        j1 = min(j1,nS)
        XcTXc_mu_psi = np.zeros((n,ell+1))
        for col in range(numCol):
            i0 = col*n
            i1 = i0+n
            i1 = min(i1,nS)

            if (row>col):
                XcTXc,_ = read_block('square',n,col,row,c)
                XcTXc = XcTXc.T
            else:
                XcTXc,_ = read_block('square',n,row,col,c)
            
            if (col==numCol-1):
                n1 = nS%n
                XcTXc = XcTXc[:,:n1]
            
            XcTXc_mu_psi = XcTXc_mu_psi+np.matmul(XcTXc,mu_psi[i0:i1,:])
            fileNum = fileNum+1
            msg = 'processing files {0:.1f}% ...'.format(fileNum/numFiles*100)
            if (fileNum%np.floor(0.1*numFiles)==0):
                print(msg)
            
        if (row==numRow-1):
            n1 = nS%n
            XcTXc_mu_psi = XcTXc_mu_psi[:n1,:]

        ATA = ATA+np.matmul(mu_psi[j0:j1,:].T,XcTXc_mu_psi)

    [S_sq,EV]=eigh(ATA)
    order = np.argsort(S_sq)[::-1]
    S_sq = S_sq[order]
    V = EV[:,order]

    S = np.sqrt(S_sq)
    invS = np.diag(1./S)
    S = np.diag(S)

    U = np.zeros((D*num_copy,ell+1))
    for cc in range(num_copy):
        U[cc*D:(cc+1)*D,:] = X1[:,c-int((c-num_copy)/2)-cc:N-int((c-num_copy)/2)-cc+1]@mu_psi@V@invS

    V = np.matmul(psi.T,V) # V has a shape of nS*(ell+1)

    write_h5('usv.h5',V,'V')
    write_h5('usv.h5',S,'S')
    write_h5('usv.h5',U,'U')
    t1 = time.time()
    print('Elapsed time {0:.2f} second'.format(t1-t0))
    
    X = reconstruct(D,np.arange(ell),np.arange(N-c+1-num_copy+1),np.arange(num_copy),U,V)
    write_h5('usv_nlsa.h5',X,'recon')

# Initialization
import os
import numpy as np
from scipy.linalg import eigh
import sys 

cxfel_root = os.environ['CXFEL_ROOT']
startup_file = cxfel_root+'misc_tools/startup.py'
exec(open(startup_file).read())

if __name__ == "__main__":
    data_file = sys.argv[1]
    variable_name = sys.argv[2]
    mu_psi_file = sys.argv[3]
    ell = sys.argv[4]
    N = sys.argv[5]
    D = sys.argv[6]
    n = sys.argv[7]
    c = sys.argv[8]
    num_copy = sys.argv[9]

    ell = int(ell)
    N = int(N)
    D = int(D)
    n = int(n)
    c = int(c)
    num_copy = int(num_copy)

    nlsa(data_file,variable_name,mu_psi_file,ell,N,D,n,c,num_copy)