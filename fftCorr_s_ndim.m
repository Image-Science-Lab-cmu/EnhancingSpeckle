function cI=fftCorr_s_ndim(I,I2,grid_ws_c)
    
    [n1,n2,~]=size(I);

%     if max(n1,n2)<=20
%         cI=conv2(I,rot90(I2,2),'same');
%         return
%     end
        
    pn1=2^nextpow2(n1);
    pn2=2^nextpow2(n2);
    
    d1=ceil((pn1-n1)/2);
    d2=ceil((pn2-n2)/2);
          
    fI=fft2((I),pn1,pn2);
    fI2=fft2((I2),pn1,pn2);        
    cI=fftshift2(ifft2((fI).*conj(fI2)));
    %cI=fftshift(ifft2(fft2(I,pn1,pn2).*conj(fft2(I2,pn1,pn2))));
    
    cI=cI(d1+[1:n1],d2+[1:n2],:);

    if exist('grid_ws_c','var')
       cI = cI(grid_ws_c,grid_ws_c,:);
    end    
