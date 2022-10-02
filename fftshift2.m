function x = fftshift2(x)
    x = fftshift(fftshift(x,1),2); %used to shift 2d in multi dimension array
end