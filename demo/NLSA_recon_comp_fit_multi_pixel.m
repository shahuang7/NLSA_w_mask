# Perform initialization &
#   add cxfel code to system path.
cxfel_root = getenv('CXFEL_ROOT')
addpath([cxfel_root '/misc_tools'])
addpath([cxfel_root '/nlsa'])

arg_list = argv();
command = sprintf('myMAT = "%s"',arg_list{1});
eval(command)

load(myMAT,'c','D','myJPEG','N','nEigs','num_copy','T_jitter_free','T_jittered','U_NLSA','S_NLSA','V_NLSA')

U_NLSA = real(U_NLSA);
S_NLSA = real(diag(S_NLSA));
V_NLSA = real(V_NLSA);

##### region of interest (both for display and for fitting) #####
x_time = [0:N-c+1-num_copy];
ROI = (c+num_copy)/2+x_time;
T_jittered = T_jittered(ROI,:);
T_jitter_free = T_jitter_free(ROI,:);

# reconstruction & fit to known

num_mode = nEigs;
X = reconstruct(D,[1:num_mode],x_time,[1:num_copy],U_NLSA,eye(num_mode),V_NLSA);
Z = permute(X,[2,3,1]);
W = ones(1,length(x_time));
[best_fit_NLSA,coeff]=Component_Fit(Z,T_jitter_free',diag(W));

err = best_fit_NLSA-T_jitter_free';
err = err(:);
err = sqrt(mean(err.*err));

##### graphics-related #####

hFigure = figure(1);
set(hFigure,'color','w')
set(hFigure,'resize','off')
Pix_SS = get(0,'screensize');
screenWidth = Pix_SS(3);
screenHeight = Pix_SS(4);
pos = [10 900 screenWidth/2 0.8*screenHeight];
try
  warning('off','Octave:abbreviated-property-match')
catch
end
set(hFigure,'pos',pos)

figure(hFigure)

numCol = ceil(sqrt(D+2));
numRow = ceil((D+2)/numCol);

x = x_time;             my_xlabel = '';

for d=1:D
  hsp = subplot(numRow,numCol,d);
  y = T_jittered(:,d);    my_ylabel = '';
  y_min = min(y); y_max = max(y);
  my_title = '';
  plotRF(hsp,x,y,my_xlabel,my_ylabel,my_title,'b.')
  hLine = get(hsp,'children');
  set(hLine,'markerSize',4)
  set(hsp,'ylim',[y_min,y_max]+[-1,1]*0.1*(y_max-y_min))
  
  y = T_jitter_free(:,d);
  addplotRF(hsp,x,y,'y-')
  hLine = get(hsp,'children');
  set(hLine,'markerSize',6)
  y = best_fit_NLSA(d,:);
  addplotRF(hsp,x,y,'r-')
  hLine = get(hsp,'children');
  set(hLine,'markerSize',8)
  
  set(hsp,'xtick',[])
end

hsp = subplot(numRow,numCol,numRow*numCol);
num_S = min([nEigs c*D]);
x = [1:num_S];          my_xlabel = 'Singular Value #';
y = S_NLSA(1:num_S);    my_ylabel = 'Magnitude (A.U.)';
y = y/y(1);
y_min = 0.0; y_max = 1.0;
my_title = 'Singular Value Spectrum';
plotRF(hsp,x,y,my_xlabel,my_ylabel,my_title,'b-o')
set(hsp,'xlim',[0.8 num_S],'ylim',[y_min,y_max],'xTick',[5:5:num_S])

print(myJPEG)
