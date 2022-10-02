function cI=fftCorr_s(I,I2)

    
    [n1,n2]=size(I);
    
    if max(n1,n2)<=20
        cI=conv2(I,rot90(I2,2),'same');
        return
    end
    
    
    pn1=2^nextpow2(n1);
    pn2=2^nextpow2(n2);
    
    d1=ceil((pn1-n1)/2);
    d2=ceil((pn2-n2)/2);
    
    
    %I=padarray(I,[ceil((n1-1)/2),ceil((n2-1)/2)]);
    %I2=padarray(I2,[ceil((n1-1)/2),ceil((n2-1)/2)]);

    
    
   %fI=fft2(ifftshift(I),pn1,pn2);
   %fI2=fft2(ifftshift(I2),pn1,pn2);
   fI=fft2((I),pn1,pn2);
   fI2=fft2((I2),pn1,pn2);
   cI=fftshift(ifft2((fI).*conj(fI2)));
   
   cI=cI(d1+[1:n1],d2+[1:n2]);
   %cI = cI(grid_ws_c,grid_ws_c);
    

