% Calculate the local correlation for each window
% u: the original wave
% stp_size: the size to shift the window, set as half w_size in recon_NF
% w_size: window size
% Wwin: the weight of window

function cu=getWinCorrW(u,stp_size,w_size,Wwin,ramp_correction)
if ~exist('Wwin','var')
    Wwin=ones(w_size,w_size);
end
hw_size=floor(w_size/2);
Wwin=padarray(Wwin,hw_size*[1,1]);
ww_size=size(Wwin,1);
hww_size=floor(ww_size/2);
[on1,on2]=size(u);
u=padarray(u,(hww_size)*[1,1]);
[n1,n2]=size(u);
hn1=floor(n1/2);
hn2=floor(n2/2);
hon1=floor(on1/2);
hon2=floor(on2/2);

grid_sx=[0:stp_size:hon1];
grid_sx=[-grid_sx(end:-1:2),grid_sx]+hn1+1;
grid_sy=[0:stp_size:hon2];
grid_sy=[-grid_sy(end:-1:2),grid_sy]+hn2+1;

Nsx=length(grid_sx);
Nsy=length(grid_sy);

if strcmp(class(u),'gpuArray')
    cu=zeros(Nsx,Nsy,w_size,w_size,'gpuArray');
else
cu=zeros(Nsx,Nsy,w_size,w_size);
end
grid_w=[-hww_size:hww_size];
grid_ws=[-hw_size:hw_size];
grid_ws_c=grid_ws+hww_size+1;

for j1=1:Nsx
    for j2=1:Nsy  
        tu=u(grid_sx(j1)+grid_w,grid_sy(j2)+grid_w);
        tuW=tu.*Wwin;
        % calculate local correlation and correct the phase
        tc=fftCorr_s(tuW,tu).*ramp_correction(hon1+grid_w,hon2+grid_w);
        cu(j1,j2,:,:)=reshape(tc(grid_ws_c,grid_ws_c),1,1,w_size,w_size);
    end
end

