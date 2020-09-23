%==========================================================================
% Range Migration Algorithm for 3-D imaging with arbitary planar MIMO array.
% In this algorithm, NUFFTs are used in two steps:
%       (1) 2D NUFFT is used to transform non-rectangular grids in space into
%       rectangular grids in the wavenumber domain;
%       (2) 1D NUFFT is used to replace the Stolt Interpolation
%
%   Ref of NUFFT:
%       L. Greengard and J.Y. Lee, Accelerating the Nonuniform Fast Fourier
%       Transform, SIAM Review, 46, 443-454(2004)
%
% Created by Jianping Wang
% Date: Feb 11, 2018
%==========================================================================
clear; close all; clc;
%% Directory of the NUFFT package 
addpath('.\nufftall_133')   % this should be modified based on your directory.
%%
addpath('.\util')
%% Load Data and Initialization
disp('Data Loading...')

load('..\Dataset\SynData_MIMO_E-shape_NumTx=70_NumRx=75.mat')

echo = EMData;   % the EM data array
freq = OperatingFreq;  % the frequency samples of the data set
%==========MIMO array===========
%transmiting array --- placed on the xz plane
TX = MIMOArray.Transmitter.TX;    % X coordinates
TY = MIMOArray.Transmitter.TY;    % Y coordinates
TZ = MIMOArray.Transmitter.TZ;    % Z coordinates
deltaSt = MIMOArray.Transmitter.vCellArea; % Area of Voronoi cell
Rd_t = 0.25;             % radius of the transmit antenna aperture
Num_Tx = length(TX);
%receiving array --- placed on the xz plane
RX = MIMOArray.Receiver.RX;       % X coordinates
RY = MIMOArray.Receiver.RY;       % Y coordinates
RZ = MIMOArray.Receiver.RZ;       % Z coordinates
deltaSr = MIMOArray.Receiver.vCellArea; %Area of Voronoi cell
Rd_r = 0.25;             %radius of the receive antenna aperture
Num_Rx = length(RX);

clear EMData
%% Some parameters/Constants
c = 299792458.000176; % speed of light

freqNum = length(freq); % number of frequency samples
Kr = 2*pi*freq/c;   % wavenumber related to each frequency
%% 2-D NUFFT over transmit and receive arrays
% The weight effects of the irregularlly spatial sampling should be taken
% into account.
disp('2-D NUFFT over transmit and receive arrays...')
%Normalizing the coordinates of transmit and receive arrays
TX_norm = TX/Rd_t*pi;   
TZ_norm = TZ/Rd_t*pi;

RX_norm = RX/Rd_r*pi;
RZ_norm = RZ/Rd_r*pi;

dKxt = 2*pi/2/Rd_t;
dKzt = 2*pi/2/Rd_t;
dKxr = 2*pi/2/Rd_r;
dKzr = 2*pi/2/Rd_r;

KxtMax = Kr;
KxtMin = 0;
Kxt_Num = floor((KxtMax - KxtMin)/dKxt)*2 + 1;
% Kxt_NumMax = max(Kxt_Num);
Kxt_NumMax = 30;
% Kxt = KxtMin + dKxt*[-(Kxt_NumMax-1)/2:(Kxt_NumMax-1)/2];
Kxt = Kxy_grid(KxtMin, Kxt_NumMax, dKxt);
Kxt_Num(Kxt_Num>Kxt_NumMax)=Kxt_NumMax;

KztMax = Kr;
KztMin = 0;
Kzt_Num = floor((KztMax - KztMin)/dKzt)*2 + 1;
% Kzt_NumMax = max(Kzt_Num);
Kzt_NumMax = 30;
% Kyt = KztMin + dKzt*[-(Kzt_NumMax-1)/2:(Kzt_NumMax-1)/2];
Kzt = Kxy_grid(KztMin, Kzt_NumMax, dKzt);
Kzt_Num(Kzt_Num>Kzt_NumMax)=Kzt_NumMax;

KxrMax = Kr;
KxrMin = 0;
Kxr_Num = floor((KxrMax - KxrMin)/dKxr)*2 + 1;
% Kxr_NumMax = max(Kxr_Num);
Kxr_NumMax = 30;
% Kxr = KxrMin + dKxr*[-(Kxr_NumMax-1)/2:(Kxr_NumMax-1)/2];
Kxr = Kxy_grid(KxrMin, Kxr_NumMax, dKxr);
Kxr_Num(Kxr_Num>Kxr_NumMax)=Kxr_NumMax;

KzrMax = Kr;
KzrMin = 0;
Kzr_Num = floor((KzrMax - KzrMin)/dKzr)*2 + 1;
% Kzr_NumMax = max(Kzr_Num);
Kzr_NumMax = 30;
% Kyr = KzrMin + dKzr*[-(Kzr_NumMax-1)/2:(Kzr_NumMax-1)/2];
Kzr = Kxy_grid(KzrMin, Kzr_NumMax, dKzr);
Kzr_Num(Kzr_Num>Kzr_NumMax)=Kzr_NumMax;

Data_FK = Space2Kdom_NUFFT2D(echo, deltaSt, deltaSr,...
    TX_norm,TZ_norm,RX_norm,RZ_norm,Kxr_Num,Kzr_Num,Kxt_Num,Kzt_Num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear echo
%% Interpolation & Data reformat
disp('Interpolation and data reformat from 5-D to 3-D')
[sig3dFK,Kx,Kz,dKy] = FK5DtoKdom3D(Kr,Kxt,Kzt,Kxr,Kzr,Data_FK);  % NUFFT
% [sig3dFK,Kx,Kz,dKy] = FK5DtoKdom3D_Interp(Kr,Kxt,Kzt,Kxr,Kzr,Data_FK); %Stolt interpolation
clear Data_FK
%% IFFT over the k-space
disp('2-D IFFT with respect to the wavenumber...')
KxNum_img = 100;
KzNum_img = 100;
KyNum = size(sig3dFK,1);

ImagData = Kdom2Space_IFFT2d(sig3dFK, 1, KxNum_img, KzNum_img);

%% Computation of Spatial sampling grid (i.e.,voxel grid) 
dKx = Kx(2)-Kx(1);
dKz = Kz(2)-Kz(1);

[X,Z] = CalKPts_Length(KxNum_img*dKx,KzNum_img*dKz,KxNum_img,KzNum_img);
Y = 2*pi/dKy/KyNum * [ 0:(KyNum-1) ];
%% Display
hf = figure;
Ftsz = 12;
maxvalue = max(abs(ImagData(:)));
ImagNorm = abs(ImagData)/maxvalue;
ImagNorm = permute(ImagNorm,[1 3 2]);
ImagNorm_dB = db(ImagNorm);

[meshx,meshy,meshz] = meshgrid(X,Y,Z);

ySliceIndex = 48:60;
hs = slice(meshx,meshy,meshz,ImagNorm_dB,[],Y(ySliceIndex),[]);
set(hs,'FaceColor','interp','EdgeColor','none');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3],'PaperSize',[4 3])
axis vis3d;
% axis([0.3 0.7 0.3 0.7  1.8 2.2])
% alpha('color')
for kk=1:length(ySliceIndex)
    set(hs(kk),'AlphaData',...
        squeeze(ImagNorm(kk+ySliceIndex(1)-1,:,:)),'FaceAlpha','interp');
end
colmp_type = 'jet';
colormap(colmp_type)
hc = colorbar
% view(12,14)
view(-16,11)
caxis([-25 0])
ylabel(hc, 'dB','fontsize',Ftsz)
daspect([1 1 1])
xlabel('X [m]','fontsize',Ftsz)
ylabel('Y [m]','fontsize',Ftsz)
zlabel('Z [m]','fontsize',Ftsz)
title(['NUFFT based RMA'],'fontsize',Ftsz)
%% Results output
print(hf, '-dpng','-r300',['..\results\Image_3D_NUFFT_RMA_NumTX=' num2str(Num_Tx) '_NumRX=' num2str(Num_Rx) '_' num2str(freq(1)/1e9) '_' num2str(freq(end)/1e9) 'GHz.png'])
print(hf, '-dpdf','-r300',['..\results\Image_3D_NUFFT_RMA_NumTX=' num2str(Num_Tx) '_NumRX=' num2str(Num_Rx) '_' num2str(freq(1)/1e9) '_' num2str(freq(end)/1e9) 'GHz.pdf'],'-opengl')
