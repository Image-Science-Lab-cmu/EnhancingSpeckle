function [x, m, fsize] = padarrays(x, m, shape)
% Pad arrays to make them the same size and allow for boundary effects
xsize = size(x);
msize = size(m);
switch shape
    
    case 'wrap'
        fsize = xsize;
        % ensure x no smaller than m
        if any(msize > xsize)  && ~isempty(x)
            x = repmat(x, ceil(msize ./ size(x)));
            xsize = size(x);
        end
        % pad m with zeros
        if any(msize < xsize)  % test, as user may have optimised already
            m = exindex(m, 1:xsize(1), 1:xsize(2), {0});
        end
        % recentre m so that y(1,1) corresponds to mask centred on x(1,1)
        mc = 1 + floor(msize/2);
        me = mc + xsize - 1;
        m = exindex(m, mc(1):me(1), mc(2):me(2), 'circular');
    
    case 'full'
        fsize = xsize + msize - 1;  % enough room for no overlap
                x = exindex(x, 1:fsize(1), 1:fsize(2), {0});
        m = exindex(m, 1:fsize(1), 1:fsize(2), {0});
    
    case 'valid'
        fsize = xsize - msize + 1;
        % pad m with zeros (don't test first, as likely to be needed)
        m = exindex(m, 1:xsize(1), 1:xsize(2), {0});
        % shift m so that y(1,1) corresponds to mask just inside x
        me = msize + xsize - 1;
        m = exindex(m, msize(1):me(1), msize(2):me(2), 'circular');
    case 'same'
        fsize = xsize;
        mmid = floor(msize/2);
        xsize = xsize + mmid;   % border to avoid edge effects
        x = exindex(x, 1:xsize(1), 1:xsize(2), {0});
        m = exindex(m, 1:xsize(1), 1:xsize(2), {0});
        % recentre m so that y(1,1) corresponds to mask centred on x(1,1)
        mc = 1 + mmid;
        me = mc + xsize - 1;
        m = exindex(m, mc(1):me(1), mc(2):me(2), 'circular');
        
    case 'reflect'
        fsize = xsize;
        xsize = xsize + msize - 1;   % border to avoid edge effects
        xc = 1 - floor((msize-1)/2);
        xe = xc + xsize - 1;
        x = exindex(x, xc(1):xe(1), xc(2):xe(2), 'symmetric');
        m = exindex(m, 1:xsize(1), 1:xsize(2), {0});
        % recentre m so that y(1,1) corresponds to mask centred on x(1,1)
        me = msize + xsize - 1;
        m = exindex(m, msize(1):me(1), msize(2):me(2), 'circular');
    otherwise
        error('conv_fft2:badshapeopt', 'Unrecognised shape option: %s', shape);
end
end
