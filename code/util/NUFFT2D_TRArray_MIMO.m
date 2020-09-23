function Data_FK = NUFFT2D_TRArray_MIMO(echo,TX_norm,TY_norm,RX_norm,RY_norm,Kxr_Num,Kyr_Num,Kxt_Num,Kyt_Num)
% Function discription:
%    This function is used to implement the Fourier transform of the measured
%  data over non-uniform 2-D MIMO arrays
%
% Parameters:
% echo: the raw data of the measured signal
% TX_norm: X coordinates of the transmit antennas,normalized within [-pi,pi]
% TY_norm: Y coordinates of the transmit antennas,[-pi,pi]
% RX_norm: X coordinates of the receive antennas, [-pi,pi]
% RY_norm: Y coordinates of the receive antennas, [-pi,pi]
% TX_Num:     %the dimension of data in Kx direction after NUFFT
% TY_Num:     %the dimension of data in Ky direction after NUFFT
% RX_Num:     % the dimension of data in Kx direction after NUFFT
% RY_Num:     % the dimension of data in Ky direction after NUFFT
%
% Date: Aug 25, 2015   J.WANG @ TU Delft, the Netherlands
%

Kxt_NumMax = max(Kxt_Num);
Kyt_NumMax = max(Kyt_Num);
Kxr_NumMax = max(Kxr_Num);
Kyr_NumMax = max(Kyr_Num);

[Num_Rx,Num_Tx,freqNum] = size(echo);

% Data_AftNFFT_Tx = zeros(Num_Rx,Kxt_NumMax,Kyt_NumMax);
% Data_FK = zeros(Kxr_NumMax,Kyr_NumMax,Kxt_NumMax,Kyt_NumMax,freqNum);
%
% for kk = 1:freqNum
%     %%NUFFT over transmit array
%     iniX = ceil((Kxt_NumMax - Kxt_Num(kk))/2 + 1);
%     iniY = ceil((Kyt_NumMax - Kyt_Num(kk))/2 + 1);
%     for tt = 1:Num_Rx
%         [fk,ier]=nufft2d1(Num_Tx,TX_norm,TY_norm,echo(tt,:,kk), -1, 1e-6, Kxt_Num(kk), Kyt_Num(kk));
%         Data_AftNFFT_Tx(tt,iniX:(iniX+Kxt_Num(kk)-1),iniY:(iniY+Kyt_Num(kk)-1)) = fk;
%     end
%
%     %%NUFFT over receive array
%     iniXr = ceil((Kxr_NumMax - Kxr_Num(kk))/2 + 1);
%     iniYr = ceil((Kyr_NumMax - Kyr_Num(kk))/2 + 1);
%     for hh = 1:Kxt_Num(kk)
%         for jj = 1:Kyt_Num(kk)
%             [fk,ier]=nufft2d1(Num_Rx,RX_norm,RY_norm,Data_AftNFFT_Tx(:,hh,jj), -1, 1e-6, Kxr_Num(kk), Kyr_Num(kk));
%             Data_FK(iniXr:(iniXr+Kxr_Num(kk)-1),iniYr:(iniYr+Kyr_Num(kk)-1),hh,jj,kk) = fk;
%         end
%     end
% end

Data_FK = zeros(Kxr_NumMax,Kyr_NumMax,Kxt_NumMax,Kyt_NumMax,freqNum);

parfor kk = 1:freqNum
    disp(['parfor ' num2str(kk) '...'])
    Data_FK_temp = zeros(Kxr_NumMax,Kyr_NumMax,Kxt_NumMax,Kyt_NumMax);
    Data_AftNFFT_Tx = zeros(Num_Rx,Kxt_NumMax,Kyt_NumMax);
    %% NUFFT over transmit array
    iniX = ceil((Kxt_NumMax - Kxt_Num(kk))/2 + 1);
    iniY = ceil((Kyt_NumMax - Kyt_Num(kk))/2 + 1);
    for tt = 1:Num_Rx
        [fk,~]=nufft2d1(Num_Tx,TX_norm,TY_norm,echo(tt,:,kk), -1, 1e-6, Kxt_Num(kk), Kyt_Num(kk));
        Data_AftNFFT_Tx(tt,iniX:(iniX+Kxt_Num(kk)-1),iniY:(iniY+Kyt_Num(kk)-1)) = fk;
    end
    
    %% NUFFT over receive array
    iniXr = ceil((Kxr_NumMax - Kxr_Num(kk))/2 + 1);
    iniYr = ceil((Kyr_NumMax - Kyr_Num(kk))/2 + 1);
    for hh = 1:Kxt_NumMax
        for jj = 1:Kyt_NumMax
            [fk,~]=nufft2d1(Num_Rx,RX_norm,RY_norm,Data_AftNFFT_Tx(:,hh,jj), -1, 1e-6, Kxr_Num(kk), Kyr_Num(kk));
            Data_FK_temp(iniXr:(iniXr+Kxr_Num(kk)-1),iniYr:(iniYr+Kyr_Num(kk)-1),hh,jj) = fk;
        end
    end
    Data_FK(:,:,:,:,kk) = Data_FK_temp;
end

