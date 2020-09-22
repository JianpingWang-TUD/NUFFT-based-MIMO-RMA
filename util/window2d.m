function  filt2d = window2d(xNum,yNum,cutFreq, sigma)
% window2d generates a 2D window.
%
% Parameter:
%  xNum --- the size of the Gaussian filter window along the x axis
%  yNum --- the size of the Gaussian filter window along the y axis
%  cutFreq --- the cut frequency of the Gaussian filter window
%  sigma --- standard deviation of the Gaussian LPF
%  
% Author: J.Wang, Aug 16, 2018
%

[f1,f2] = freqspace([xNum yNum],'meshgrid');
Hd = ones([xNum yNum]); 
r = sqrt(f1.^2 + f2.^2);
Hd((r>cutFreq)) = 0;

win = fspecial('gaussian',[xNum yNum],sigma);
win = win ./ max(win(:));  % Normalize the maximum of the window be 1.

h = fwind2(Hd,win);
[filt2d, fv1 , fv2 ]= freqz2(h,[yNum xNum]);