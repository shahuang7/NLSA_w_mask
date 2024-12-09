function X=reconstruct(D,kOfInterest,tOfInterest,copyOfInterest,U,S,V)
  
  num_k = length(kOfInterest);
  num_t = length(tOfInterest);
  num_cc = length(copyOfInterest);
  X = zeros(num_k,D,num_t);
  for i3=1:num_t
    t = tOfInterest(i3);
    mesg = sprintf('time point# %d of %d ...',i3,num_t);
    if (mod(i3,floor(0.1*num_t))==0), disp(mesg); end
    for i1=1:num_k
      k = kOfInterest(i1);
      Xk = zeros(D,1);
      for cc=copyOfInterest
        Xk = Xk+U((cc-1)*D+[1:D],k)*S(k,k)*V(t+cc,k);
      end
      X(i1,:,i3) = Xk/num_cc;
    end
  end
% end function reconstruct
