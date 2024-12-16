function h_filt = Get4dSmoothingKernel(smoothing_stencil, filter_type)

if nargin<2
    filter_type = 'avg';
end

if strcmp(filter_type,'avg')
    x1 = -(smoothing_stencil(1)-1)/2:(smoothing_stencil(1)-1)/2;
    x2 = -(smoothing_stencil(2)-1)/2:(smoothing_stencil(2)-1)/2;
    x3 = -(smoothing_stencil(3)-1)/2:(smoothing_stencil(3)-1)/2;
    x4 = -(smoothing_stencil(4)-1)/2:(smoothing_stencil(4)-1)/2;
    [X1, X2, X3, X4] = ndgrid(x1,x2,x3,x4);
    h_filt = ((X1./max(X1(:))).^2+...
        (X2./max(X2(:))).^2+...
        (X3./max(X3(:))).^2)<=1;
    %            (X4./max(X4(:))).^2)<=1;
    % avg_filt = ones(smoothing_stencil);
    h_filt = h_filt./sum(h_filt(:));
elseif strcmp(filter_type,'gaussian')
    % Centered Gaussian lowpass
    % filter of size SIZE with standard deviations defined as
    % SIZE/(4*sqrt(2*log(2))) so that FWHM equals half filter size
    % (http://en.wikipedia.org/wiki/FWHM). Such a FWHM-dependent standard
    % deviation yields a congruous Gaussian shape (what should be expected
    % for a Gaussian filter!).
    
    sig = smoothing_stencil/(4*sqrt(2*log(2)));
    siz   = (smoothing_stencil-1)/2;
    
    [x,y,z,t] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3),-siz(4):siz(4));
    h_filt = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2 + t.*t/2/sig(4)^2));
    h_filt = h_filt/sum(h_filt(:));
else
    error('Invalid filter type')
end
return
