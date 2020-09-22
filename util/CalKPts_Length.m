function [Kx,Ky] = CalKPts_Length(Lx,Ly,Kx_Num,Ky_Num)
% Function description:
%       This function is used to calculated the sampling points in the
%       wavenumber domain.
%
% Parameter?
%   Lx: Length of antenna aperture in the X direction
%   Ly: Length of antenna aperture in the Y direction
%   Kx_Num: the Number of sampling points in the Kx direction
%   Ky_Num: the Number of sampling points in the Ky direction
%
% Created by J.WANG, TU Delft
% Date: Aug 25,2015
%

if mod(Kx_Num,2)==1
    Kx = 2*pi/Lx * [-(Kx_Num - 1)/2 : (Kx_Num - 1)/2];
else
    Kx = 2*pi/Lx * [-Kx_Num/2 : (Kx_Num/2 - 1)];
end
if mod(Ky_Num,2)==1
    Ky = 2*pi/Ly * [-(Ky_Num - 1)/2 : (Ky_Num - 1)/2];
else
    Ky = 2*pi/Ly * [-Ky_Num/2 : (Ky_Num/2 - 1)];
end
