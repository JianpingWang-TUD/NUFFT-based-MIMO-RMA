function ImagData = Kdom2Space_IFFT2d(sigK3d, isWindw, KxNum_img, KyNum_img)
% Kdom2Space_IFFT2d transforms the K-space data array to the focused
% image by using 2D IFFT.
%
% Parameter:
%  sigK3d  --- signal data array in the wavenumber domain 
%  KxNum_img --- the point-number of IFFT along the Kx axis, i.e., the
%                pixel number of the spatial image in the x direction;
%  KyNum_img --- the point-number of IFFT along the Ky axis, i.e., the
%                pixel number of the spatial image in the y direction;
%  isWindw   --- the sign for windowing, true for windowing; otherwise false.
% Output:
%  ImagData  --- the focused image
%
% Author: J.Wang, Aug 15, 2018
%

[KzNum,KyNum,KxNum] = size(sigK3d);
if nargin<3
   KxNum_img = KxNum;
   KyNum_img = KyNum;
end

if KxNum_img<KxNum || KyNum_img<KyNum
    error(['Error. Input point-number of IFFT should be greater or equal to the point-number of the data set.'])
end

ImagData = zeros(KzNum, KyNum_img, KxNum_img);
%%
sigK3d = ifftshift( sigK3d, 1);
sigK3d_ifftshift_y = zeros(KzNum, KyNum_img, KxNum);
if mode(KyNum,2)==0
    sigK3d_ifftshift_y(:, 1:KyNum/2, :) = sigK3d(:, KyNum/2+1:end, :) ;
    sigK3d_ifftshift_y(:, end-KyNum/2+1:end, :) = sigK3d(:, 1:KyNum/2, :);
else
    sigK3d_ifftshift_y(:, 1:(KyNum+1)/2, :) = sigK3d(:, (KyNum+1)/2:end, :) ;
    sigK3d_ifftshift_y(:, end-(KyNum-1)/2+1:end, :) = sigK3d(:, 1:(KyNum-1)/2, :);
end
sigK3d_ifftshift_yx = zeros(KzNum, KyNum_img, KxNum_img);
if mode(KxNum,2)==0
    sigK3d_ifftshift_yx(:, :, 1:KxNum/2) = sigK3d_ifftshift_y(:, :, KxNum/2+1:end) ;
    sigK3d_ifftshift_yx(:, :, end-KxNum/2+1:end) = sigK3d_ifftshift_y(:, :, 1:KxNum/2);
else
    sigK3d_ifftshift_yx(:, :, 1:(KxNum+1)/2) = sigK3d_ifftshift_y(:, :, (KxNum+1)/2:end) ;
    sigK3d_ifftshift_yx(:, :, end-(KxNum-1)/2+1:end) = sigK3d_ifftshift_y(:, :, 1:(KxNum-1)/2);
end
%%
if isWindw
    filt2d = window2d(KyNum,KxNum,0.85, 2);
    filt2d_ifftshift_y = zeros(KyNum_img, KxNum);
    if mode(KyNum,2)==0
        filt2d_ifftshift_y(1:KyNum/2,:) = filt2d(KyNum/2+1:end,:);
        filt2d_ifftshift_y(end-KyNum/2+1:end, :) = filt2d(1:KyNum/2, :);
    else
        filt2d_ifftshift_y(1:(KyNum+1)/2, :) = filt2d((KyNum+1)/2:end, :);
        filt2d_ifftshift_y(end-(KyNum-1)/2+1:end, :) = filt2d(1:(KyNum-1)/2, :);
    end
    filt2d_ifftshift_yx = zeros(KyNum_img, KxNum_img);
    if mode(KxNum,2)==0
        filt2d_ifftshift_yx(:,1:KxNum/2) = filt2d_ifftshift_y(:, KxNum/2+1:end);
        filt2d_ifftshift_yx(end-KxNum/2+1:end) = filt2d_ifftshift_y(:, 1:KxNum/2);
    else
        filt2d_ifftshift_yx(:, 1:(KxNum+1)/2) = filt2d_ifftshift_y(:, (KxNum+1)/2:end);
        filt2d_ifftshift_yx(:, end-(KxNum-1)/2+1:end) = filt2d_ifftshift_y(:, 1:(KxNum-1)/2);
    end
    for kk = 1 : KzNum
        ImagData(kk,:,:) = fftshift( ifft2( squeeze(sigK3d_ifftshift_yx(kk,:,:))...
            .* filt2d_ifftshift_yx  ) );
    end
else
    for kk = 1 : KzNum
        ImagData(kk,:,:) = fftshift( ifft2( squeeze(sigK3d_ifftshift_yx(kk,:,:)) ) );
    end
end
%%


