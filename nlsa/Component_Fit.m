function [best_fit,coeff]=Component_Fit(Z,T,W)
  [num_pixel,wnd_width,num_mode] = size(Z);
  LHS = zeros(num_mode);
  RHS = zeros(num_mode,1);
  for dd=1:num_pixel
    Zd = squeeze(Z(dd,:,:));
    LHS = LHS+Zd'*W*Zd;
    RHS = RHS+Zd'*W*T(dd,:)';
  end
  LHS_01 = squeeze(sum(bsxfun(@times,Z,sum(W)),2))';
  if (num_pixel==1), LHS_01 = LHS_01'; end
  LHS_11 = sum(W(:))*eye(num_pixel);
  RHS_1 = T*diag(W);
  LHS = [LHS LHS_01; LHS_01' LHS_11];
  RHS = [RHS;RHS_1];
  coeff = LHS\RHS;
  best_fit = zeros(num_pixel,wnd_width);
  for jj=1:num_mode
    best_fit = best_fit+squeeze(Z(:,:,jj))*coeff(jj);
  end
  best_fit = bsxfun(@plus,best_fit,coeff(num_mode+1:end));
%end function Component_Fit
