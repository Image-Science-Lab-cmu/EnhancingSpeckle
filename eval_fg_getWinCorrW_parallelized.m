% This function evaluates the error function and its gradient with respect to the current solution. 
% The error is the squared difference between a target local correlation
% array and the one produced by the current guess.
function [f,g]=eval_fg_getWinCorrW_parallelized(u,stp_size,w_size,Wwin,targ,bdr_size,suppressDC)
if ~exist('Wwin','var')
    Wwin=ones(w_size,w_size);
end
hw_size=floor(w_size/2);
Wwin=padarray(Wwin,hw_size*[1,1]);
ww_size=size(Wwin,1);
hww_size=floor(ww_size/2);
u   = padarray(u,(bdr_size)*[1,1]);
on  = size(u); hon = floor(on/2);
u   = padarray(u,(hww_size)*[1,1]);
n   = size(u); hn  = floor(n/2);

grid_sx=[0:stp_size:hon(1)];
grid_sx=[-grid_sx(end:-1:2),grid_sx]+hn(1)+1;
grid_sy=[0:stp_size:hon(2)];
grid_sy=[-grid_sy(end:-1:2),grid_sy]+hn(2)+1;

Nsx=length(grid_sx);
Nsy=length(grid_sy);
if strcmp(class(u),'gpuArray')
    cu=zeros(w_size,w_size,'gpuArray');
    g=zeros(n,'gpuArray');
else
    cu=zeros(w_size,w_size);
    g=zeros(n);
end
grid_w=[-hww_size:hww_size];
grid_ws=[-hw_size:hw_size];
grid_ws_c=grid_ws+hww_size+1;
f=0;

flipfull = @(x) fliplr(flipud(x));

%Nsx
batch_size = 1;
parallized= true;
M = ones(2*hw_size+1, 2*hw_size+1,'gpuArray');

% The DC component in autocorrealtion is not important
% Since beads have blurry sizes, we suppressed a  blurred region of DC
if suppressDC
    %M(hw_size+1,hw_size+1) = 0;
    M(hw_size+1+[-2:2],hw_size+1+[-2:2])= 1-conv2(ones(3),ones(3))/9;
end


if parallized
    tus   =  zeros(2*hww_size+1,2*hww_size+1, batch_size,'gpuArray');
    targs =  zeros(2*hw_size+1, 2*hw_size+1,  batch_size,'gpuArray');
    errf  =  zeros(2*hww_size+1,2*hww_size+1, batch_size,'gpuArray');
    
    [jx,jy] = meshgrid(1:Nsx,1:Nsy);
    
    for batch_id = 1: ceil(length(jx(:))/batch_size)           
        for jid= 1:batch_size
            full_id = jid + (batch_id-1)*batch_size;
            if(full_id<= Nsx* Nsy)
                j1= jx(full_id); j2=jy(full_id);
                %fprintf('%d,%d\n',j1,j2);                
                tus(:,:,jid)   = u(grid_sx(j1)+grid_w,grid_sy(j2)+grid_w);
                targs(:,:,jid) = squeeze(targ(j1,j2,:,:));                
            else                
                tus(:,:,jid:batch_size)   = 0;
                targs(:,:,jid:batch_size) = 0;
                jid = jid-1;
                break;
            end
        end
        tuW     = tus.*Wwin;        
        errf(grid_ws_c,grid_ws_c,:) = M.*(fftCorr_s_ndim(tuW,tus,grid_ws_c) - targs);
        
        f       = f+sum(abs(errf(:)).^2);        
        tgs     = 2*( Wwin.*fftCorr_s_ndim(errf,flipfull(tus))+ flipfull(fftCorr_s_ndim(errf,tuW)));

        for jid= 1:batch_size
            full_id = jid + (batch_id-1)*batch_size;            
            if(full_id<= Nsx* Nsy)            
                j1= jx(full_id); j2=jy(full_id);
                g(grid_sx(j1)+grid_w,grid_sy(j2)+grid_w)=g(grid_sx(j1)+grid_w,grid_sy(j2)+grid_w)+tgs(:,:,jid);                
            else
                break;
            end
        end    
    end
else
    for j1=1:Nsx
        for j2=1:Nsy
            %fprintf('%d,%d\n',j1,j2);
            tu=u(grid_sx(j1)+grid_w,grid_sy(j2)+grid_w);
            tuW=tu.*Wwin;        
            tc = fftCorr_s(tuW,tu);
            errf= tc(grid_ws_c,grid_ws_c) - squeeze(targ(j1,j2,:,:));
            f = f+sum(abs(errf(:)).^2);
            errf=padarray(errf,(hww_size-hw_size)*[1,1]);

            tg=2*Wwin.*fftCorr_s(errf,flip(tu));
            tg=tg+2*flip(fftCorr_s(errf,Wwin.*tu));             
            
            g(grid_sx(j1)+grid_w,grid_sy(j2)+grid_w)=g(grid_sx(j1)+grid_w,grid_sy(j2)+grid_w)+tg;
        end
    end
end
dn = (n-on)/2+bdr_size;

g=g(dn(1)+1:end-dn(1),dn(2)+1:end-dn(2));

g=real(g);
f=real(f);

end
