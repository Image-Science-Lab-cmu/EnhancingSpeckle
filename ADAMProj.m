
% an ADAM optimizer to optimize the function
% f: function to be optimized
% x: initial values
% stp_size: initial step size
% maxItr: max iterations

function [x,e]=ADAMProj(f, x, stp_size, maxItr,varargin)

thr_norm= 1e-6;

beta1=0.9;
beta2=0.999;

m=zeros(size(x));
v=zeros(size(x));
eps=0.1^7;

stp_size0=stp_size;
pe=1e16; g=zeros(size(x));
px=x;
tic

errV=zeros(1,maxItr,'gpuArray');
itrCtrV=zeros(1,maxItr,'gpuArray');
stp_size_t=stp_size;
for itr=1:maxItr
    fprintf("iter = %d,\n",itr);
    m=beta1*m+(1-beta1)*g;
    v=beta2*v+(1-beta2)*abs(g).^2;
    if itr>1
        m_hat=m./(1-beta1.^(itr-1));
        v_hat=v./(1-beta2.^(itr-1));
    else
        m_hat=m; v_hat=v;
    end
     
    tg=stp_size*m_hat./sqrt(v_hat+eps);
    x=max(x-tg,0);
    [e,g]=f(x);
    eps_t=eps; stp_size_t=stp_size;
    
    %if the error has increased, reduce step size
    while ((e>pe*1.000001)|(norm(x(:))^2<thr_norm))
        itrCtrV(itr)=itrCtrV(itr)+1;
        eps_t=eps_t*2;
        stp_size_t=stp_size_t/2;
        fprintf('reduced step size %f\n',stp_size_t);
        m=pm; v=pv; g=pg; x=px; 
        m=beta1*m+(1-beta1)*g;
        v=beta2*v+(1-beta2)*abs(g).^2;
        m_hat=m/(1-beta1.^(itr-1));
        v_hat=v/(1-beta2.^(itr-1));
        tg=stp_size_t*m_hat./sqrt(v_hat+eps_t);
        %constraint solution to be positive.
        %stp_size=stp_size_t;
        x=max(x-tg,0);
        [e,g]=f(x);
    end
    
    errV(itr)=e;
    
    r=norm(x(:)-px(:))/norm(tg(:));
    if r<0.1^4
        r
        return
    end

    px=x; pg=g;
    pm=m; pv=v;
    pe=e;
    
    
    if (mod(itr,5)==0)
        tsCtr=sum(itrCtrV(itr-5+1:itr));
        %if error did not increase in any of the last iterations we can
        %safely increase step size
        if tsCtr==0
            stp_size=stp_size*1.4
        end
        %if error increased in recent iterations, need to decrease step
        %size. 
        if tsCtr>5
            stp_size=stp_size/2
        end
    end
    
    if (mod(itr,50)==1)
        toc
        tic
        itr
        r
        e
        if (itr>5)
            sum(itrCtrV(itr-5+1:itr))
        end
    end

    
    figure(41), imagesc(gather(real(x))); colormap(hot); drawnow;

end












