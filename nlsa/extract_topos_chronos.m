function [U,S,V]=extract_topos_chronos(ell,X1,mu,psi,D,N,n,c,num_copy)
  
  psi = double(psi);
  mu = double(mu);
  nS = size(psi,1);
  psi = [ones(nS,1) psi(:,1:ell)];
  mu_psi = bsxfun(@times,mu,psi);
  
  ATA = find_ATA(ell,mu_psi,N,n,c);
  
  [EV,S_sq] = eig(ATA);
  [S_sq,order] = sort(diag(S_sq),'descend');
  V = EV(:,order);
  
  S = sqrt(S_sq);
  invS = diag(1./S);
  S = diag(S);
  
  U = zeros(D*num_copy,ell+1);
  for cc=1:num_copy
    U((cc-1)*D+[1:D],:) = X1(:,[c:end]-(c-num_copy)/2-cc+1)*mu_psi*V*invS;
  end
  
  V = psi*V;
% end function extract_topos_chronos
