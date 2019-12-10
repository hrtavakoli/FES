
function saliency  = calculateImageSaliency(img, nSample, radius, sigma0, sigma1, ph1)
%
% function saliency  = calculateImageSaliency(image, nSample, radius, sigma0, sigma1, ph1)
%
% calculates the saliency of an image at a given scale
%
% @param
%   img : an image
%   nSample : number of samples in the center
%   radius : radius of the circular samples are going to be take from
%   sigma0 : standard deviation of kernels in surround
%   sigma1 : standard deviation of kernel in center
%   ph1 : P(1|x)
% 
% please refer to the following paper for details
% Rezazadegan Tavakoli H, Rahtu E & Heikkilä J, 
% "Fast and efficient saliency detection using sparse sampling and kernel density estimation."
% Proc. Scandinavian Conference on Image Analysis (SCIA 2011), 2011, Ystad, Sweden.
%
% The code has been tested on Matlab 2010a (32-bit) running windows. 
% This code is publicly available for demonstration and educational
% purposes, any commercial use without permission is strictly prohibited.  
%
% Please contact the author in case of any questions, comments, or Bug
% reports
%
% @CopyRight: Hamed Rezazadegan Tavakoli
% @Contact Email: hrezazad@ee.oulu.fi
% @date  : 2010
% @version: 0.1



[nrow, ncol, nChannel] = size(img);
if nChannel~=3 
    error('This works only on images of 3 channel preferably LAB color space needed')
end

fMat = reshape(double(img), [nrow*ncol, nChannel]); %Build feature matrix


%%
% get center coordinates
[yc xc] = meshgrid(1:ncol, 1:nrow);
xc = xc(:);
yc = yc(:);


%%
% process each sample
x=round(radius*cos((1:nSample)*2*pi/nSample));
y=round(-radius*sin((1:nSample)*2*pi/nSample));

LxcH0=zeros(numel(xc),1);
for i = 1:nSample
    fx = min(max(xc + x(i), 1), nrow);
    fy = min(max(yc + y(i), 1), ncol);
    Ind= fx+(fy-1)*nrow;
    fMatc = fMat(Ind, :);
    
    LxcH0=LxcH0+exp(-sum((fMatc-fMat).^2,2)/(2*(sigma0^2)));        
end

 ph0 = 1 - ph1;
 const=(ph0(:)*sigma1)./(ph1(:)*nSample*sigma0);

    
%% Compute final posterior
PH1xc=1./(1+const.*LxcH0);
saliency=reshape(PH1xc,[nrow,ncol]);
