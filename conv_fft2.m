function y = conv_fft2(x, m, shape)
%CONV_FFT2 Two dimensional convolution via the FFT.
%   Y = CONV_FFT2(X, M) performs the 2-D convolution of matrices X and M.
%   If [mx,nx] = size(X) and [mm,nm] = size(M), then size(Y) =
%   [mx+mm-1,nx+nm-1]. Values near the boundaries of the output array are
%   calculated as if X was surrounded by a border of zero values.
%
%   Y = CONV_FFT2(X, M, SHAPE) where SHAPE is a string returns a
%   subsection of the 2-D convolution with size specified by SHAPE:
%
%       'full'    - (default) returns the full 2-D convolution,
%       'same'    - returns the central part of the convolution
%                   that is the same size as X (using zero padding),
%       'valid'   - returns only those parts of the convolution
%                   that are computed without the zero-padded
%                   edges, size(Y) = [mx-mm+1,nx-nm+1] when
%                   size(X) > size(M),
%       'wrap'    - as for 'same' except that instead of using
%                   zero-padding the input X is taken to wrap round as
%                   on a toroid.
%       'reflect' - as for 'same' except that instead of using
%                   zero-padding the input X is taken to be reflected
%                   at its boundaries.
%
%   For shape options 'full', 'same' and 'valid' this function should
%   return the same result as CONV2, to within a small tolerance. For all
%   shape options, this function should return the same result as CONVOLVE2
%   (available via the MATLAB Central File Exchange), with no TOL argument,
%   to within a small tolerance.
%
%   CONV_FFT2 uses multiplication in the frequency domain to compute the
%   convolution. It may be faster than CONV2 and CONVOLVE2 for masks above
%   a certain size. This should be checked experimentally for any given
%   application and system.
%
%   See also CONV2, CONVOLVE2, FILTER2.
% Copyright David Young, April 2011
error(nargchk(2,3,nargin,'struct'));
if nargin < 3
    shape = 'full';    % shape default as for CONV2
end;
[x, m, fsize] = padarrays(x, m, shape);
% no need to trap case of real x and m - fft2 handles efficiently
y = ifft2(fft2(x) .* fft2(m));   % central operation, basic form
% trim to correct output size
if ~isequal(fsize, size(y))
    y = y(1:fsize(1), 1:fsize(2));
end
end