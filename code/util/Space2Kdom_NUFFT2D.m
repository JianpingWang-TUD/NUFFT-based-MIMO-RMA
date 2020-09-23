function Data_FK = Space2Kdom_NUFFT2D(echo, deltaSt, deltaSr,...
    TX_norm,TY_norm,RX_norm,RY_norm,Krx_Num,Kry_Num,Ktx_Num,Kty_Num)
% Space2Kdom_NUFFT2D transforms the echo from space domain to wavenumber
% domain by using 2D NUFFT over transmit and receive apertures.
%
% Parameter:
%   echo    ---  the echo data array, R_ant * T_ant * freq
%   deltaSt ---  the area of the Voronoi cell around each transmit antenna
%   deltaSr ---  the area of the Voronoi cell around each receive antenna
%   TX_norm ---  the normalized x coordinates of transmit antennas
%   TY_norm ---  the normalized y coordinates of transmit antennas
%   RX_norm ---  the normalized x coordinates of receive antennas
%   RY_norm ---  the normalized y coordinates of receive antennas
%   Krx_Num ---  number of Kx samples over the receive aperture
%   Kry_Num ---  number of Ky samples over the receive aperture
%   Ktx_Num ---  number of Kx samples over the transmit aperture
%   Kty_Num ---  number of Ky samples over the transmit aperture
%
% Output:
%   Data_FK --- the data in the F-K domain
%
% Author: J. Wang, Aug 25, 2015 
% 

freqNum = size(echo,3);
% Spatial weighting effect
deltaS = kron(deltaSr.',deltaSt);
deltaS = repmat(deltaS,[1,1, freqNum]);
echo = echo .* deltaS;

Data_FK = NUFFT2D_TRArray_MIMO(echo,TX_norm,TY_norm,RX_norm,RY_norm,...
    Krx_Num,Kry_Num,Ktx_Num,Kty_Num);

