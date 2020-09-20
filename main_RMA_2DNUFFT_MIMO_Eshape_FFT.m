%==========================================================================
% Range Migration Algorithm for MIMO array 3-D imaging.
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
addpath('D:\jianpingwang\Ref\PhD Topic\circular array imaging\RMA\RMA_NUFFT\nufftall_133')

% load('D:\jianpingwang\Ref\PhD Topic\circular array imaging\NUFFT_MIMO_Simulation\Elefield_MIMO_NumTx=9_NumRx=16.mat')
% load('.\MIMO_dense\raw_data\SiemensStar_script_MIMO3mm_EZ_antNumT_70_antNumR_75.mat')
load('D:\jianpingwang\Ref\PhD Topic\circular array imaging\NUFFT_MIMO_Simulation\NUFFT_MIMO Array_revise\Jan_2018\MIMO_dense\E-shape\E_shape_MIMO_antNumT_70_antNumR_75.mat')
load('D:\jianpingwang\Ref\PhD Topic\circular array imaging\NUFFT_MIMO_Simulation\NUFFT_MIMO Array_revise\Jan_2018\MIMO_dense\SiemenseStar\vArea_MIMO_SiemenseStar_Num_T=70_Num_R=75.mat')

deltaSt = vArea_T;
deltaSr = vArea_R;
%% Parameters
c = 299792458.000176;
ii = sqrt(-1);   % unit of the imagninary number

freql = 5e9; %low frequency
freqh = 15e9; %high frequency
freqstep = 100e6;   % frequency steps
freqNum = floor((freqh-freql)/freqstep)+1; %number of frequency samples
freq = freql + freqstep * (0:freqNum - 1);
Kr = 2*pi*freq/c;   % wavenumber related to each frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%transmitting antenna coordinates
xy_t = load('D:\jianpingwang\Ref\PhD Topic\circular array imaging\NUFFT_MIMO_Simulation\NUFFT_MIMO Array_revise\Jan_2018\MIMO_dense\Ant_transmit_Coord.dat');
xy_r = load('D:\jianpingwang\Ref\PhD Topic\circular array imaging\NUFFT_MIMO_Simulation\NUFFT_MIMO Array_revise\Jan_2018\MIMO_dense\Ant_receive_Coord.dat');
TX = xy_t(:,1).';
TY = xy_t(:,2).';

RX = xy_r(:,1).';
RY = xy_r(:,2).';

Rtx = 0.25;             % radius of the transmit antenna aperture
Rrx = 0.25;             %radius of the receive antenna aperture

%%the down-range coordinates of transmit and receive arrays
TantZ = 0;
RantZ = 0;

Num_Tx = length(TX);   % number of transmit antennas
Num_Rx = length(RX);   % number of receive antennas

%% Echo Loading
disp('Echo Loading...')

flInd = (freql-5e9)/freqstep+1;
fhInd = (freqh-5e9)/freqstep + 1;

% echo = zeros(Num_Rx,Num_Tx,fhInd-flInd+1);
echo = EMComponent.EZ(:,:,flInd:1:fhInd);

clear EX EY EZ
%% 2-D NUFFT over transmit and receive arrays
% The weight effects of the irregularlly spatial sampling should be taken
% into account.
disp('2-D NUFFT over transmit and receive arrays...')

%%Normalizing the coordinates of transmit and receive arrays
TX_norm = TX/Rtx*pi;   
TY_norm = TY/Rtx*pi;

RX_norm = RX/Rrx*pi;
RY_norm = RY/Rrx*pi;

dKxt = 2*pi/2/Rtx;
dKyt = 2*pi/2/Rtx;
dKxr = 2*pi/2/Rrx;
dKyr = 2*pi/2/Rrx;

KxtMax = Kr;
KxtMin = 0;
Kxt_Num = floor((KxtMax - KxtMin)/dKxt)*2 + 1;
% Kxt_NumMax = max(Kxt_Num);
Kxt_NumMax = 30;
% Kxt = KxtMin + dKxt*[-(Kxt_NumMax-1)/2:(Kxt_NumMax-1)/2];
Kxt = Kxy_grid(KxtMin, Kxt_NumMax, dKxt);
Kxt_Num(Kxt_Num>Kxt_NumMax)=Kxt_NumMax;

KytMax = Kr;
KytMin = 0;
Kyt_Num = floor((KytMax - KytMin)/dKyt)*2 + 1;
% Kyt_NumMax = max(Kyt_Num);
Kyt_NumMax = 30;
% Kyt = KytMin + dKyt*[-(Kyt_NumMax-1)/2:(Kyt_NumMax-1)/2];
Kyt = Kxy_grid(KytMin, Kyt_NumMax, dKyt);
Kyt_Num(Kyt_Num>Kyt_NumMax)=Kyt_NumMax;

KxrMax = Kr;
KxrMin = 0;
Kxr_Num = floor((KxrMax - KxrMin)/dKxr)*2 + 1;
% Kxr_NumMax = max(Kxr_Num);
Kxr_NumMax = 30;
% Kxr = KxrMin + dKxr*[-(Kxr_NumMax-1)/2:(Kxr_NumMax-1)/2];
Kxr = Kxy_grid(KxrMin, Kxr_NumMax, dKxr);
Kxr_Num(Kxr_Num>Kxr_NumMax)=Kxr_NumMax;

KyrMax = Kr;
KyrMin = 0;
Kyr_Num = floor((KyrMax - KyrMin)/dKyr)*2 + 1;
% Kyr_NumMax = max(Kyr_Num);
Kyr_NumMax = 30;
% Kyr = KyrMin + dKyr*[-(Kyr_NumMax-1)/2:(Kyr_NumMax-1)/2];
Kyr = Kxy_grid(KyrMin, Kyr_NumMax, dKyr);
Kyr_Num(Kyr_Num>Kyr_NumMax)=Kyr_NumMax;

Data_FK = Space2Kdom_NUFFT2D(echo, deltaSt, deltaSr,...
    TX_norm,TY_norm,RX_norm,RY_norm,Kxr_Num,Kyr_Num,Kxt_Num,Kyt_Num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear echo
%% Interpolation & Data reformat
disp('Interpolation and data reformat from 5-D to 3-D')
tic
% [sig3dFK,Kx,Ky,dKz] = Interp5Dto3D_2018_mtr(Kr,Kxt,Kyt,Kxr,Kyr,Data_FK);
% [sig3dFK,Kx,Ky,dKz] = Interp5Dto3D_20180410_mtr(Kr,Kxt,Kyt,Kxr,Kyr,Data_FK);
[sig3dFK,Kx,Ky,dKz] = FK5DtoKdom3D(Kr,Kxt,Kyt,Kxr,Kyr,Data_FK);
toc

clear Data_FK
%% IFFT over the k-space
disp('2-D IFFT with respect to the wavenumber...')
KyNum_img = 100;
KxNum_img = 100;

ImagData = Kdom2Space_IFFT2d(sig3dFK, KxNum_img, KyNum_img, 1);
toc
%%
dKx = Kx(2)-Kx(1);
dKy = Ky(2)-Ky(1);

% [X,Y] = CalKPts_Length(KxNum*dKx,KyNum*dKy,KxNum,KyNum);
[X,Y] = CalKPts_Length(KxNum_img*dKx,KyNum_img*dKy,KxNum_img,KyNum_img);
Z = 2*pi/dKz/KzNum * [ 0:(KzNum-1) ];
%% Display
hand_fig = figure;
Ftsz = 12;
maxvalue = max(abs(ImagData(:)));
ImagNorm = abs(ImagData)/maxvalue;
ImagNorm = permute(ImagNorm,[1 3 2]);
ImagNorm_dB = db(ImagNorm);

[meshx,meshz,meshy] = meshgrid(X,Z,Y);

zSliceIndex = 48:60;

hs = slice(meshx,meshz,meshy,ImagNorm_dB,[],Z(zSliceIndex),[]);
set(hs,'FaceColor','interp','EdgeColor','none');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3],'PaperSize',[4 3])
axis vis3d;
% axis([0.3 0.7 0.3 0.7  1.8 2.2])
% alpha('color')
for kk=1:length(zSliceIndex)
    set(hs(kk),'AlphaData',...
        squeeze(ImagNorm(kk+zSliceIndex(1)-1,:,:)),'FaceAlpha','interp');
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
%%