function  filt2d = window2d(xNum,yNum,cutFreq, sigma)

[f1,f2] = freqspace([xNum yNum],'meshgrid');
Hd = ones([xNum yNum]); 
r = sqrt(f1.^2 + f2.^2);
Hd((r>cutFreq)) = 0;

win = fspecial('gaussian',[xNum yNum],sigma);
win = win ./ max(win(:));  % Make the maximum window value be 1.

h = fwind2(Hd,win);
[filt2d, fv1 , fv2 ]= freqz2(h,[yNum xNum]);