function [Data,Kx,Ky,dKz]= FK5DtoKdom3D_Interp(Kr,Kxt,Kyt,Kxr,Kyr,Data_FK)
% FK5DtoKdom3D_Interp implements the 5-D F-K domain data to 3-D K-space
% mapping through Stolt interpolation.
% 
% Parameter:
%  Kr --- the wavenumber
%  Kxt --- the wavenumber along one of the axes(e.g., x-axis) of the
%          transmitting array
%  Kyt --- the wavenumber along the other of the axes(e.g., y-axis) of the
%          transmitting array
%  Kxr --- the wavenumber along the x-axis of the receiving array
%  Kyr --- the wavenumber along the y-axis of the receiving array
%  Data_FK --- the F-K domain data array
%
% Author: J.Wang, @MS3, TU Delft
%


[Kx_t,Ky_t] = meshgrid(Kxt,Kyt);
[Kx_r,Ky_r] = meshgrid(Kxr,Kyr);

Kxy_t_sq = Kx_t.^2 + Ky_t.^2;
Kxy_r_sq = Kx_r.^2 + Ky_r.^2;

[Kxr_Num,Kyr_Num,Kxt_Num,Kyt_Num,freqNum]=size(Data_FK);

KzMtr_t = zeros(freqNum,Kyt_Num,Kxt_Num);
KzMtr_r = zeros(freqNum,Kyr_Num,Kxr_Num);
for kk = 1:freqNum
    KzMtr_t(kk,:,:) = real(sqrt( Kr(kk)^2 - Kxy_t_sq ));
    KzMtr_r(kk,:,:) = real(sqrt( Kr(kk)^2 - Kxy_r_sq ));
end

Kzt_min = min(KzMtr_t(:));
Kzt_max = max(KzMtr_t(:));

Kzr_min = min(KzMtr_r(:));
Kzr_max = max(KzMtr_r(:));

Kzmin = Kzt_min + Kzr_min;
Kzmax = Kzt_max + Kzr_max;

dKz = 2*(Kr(2)-Kr(1));%(Kr(2)-Kr(1));

KzNum = floor((Kzmax - Kzmin)/dKz) + 1;   % could be chosen based on the system parameter for acceleration
Kz = (Kzmin + dKz * (0:KzNum-1)).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[KxtMtr, KxrMtr] = meshgrid(Kxt,Kxr);
[KytMtr, KyrMtr] = meshgrid(Kyt,Kyr);
KxMtr = KxtMtr + KxrMtr;
KyMtr = KytMtr + KyrMtr;

KxMtr = round(KxMtr/1e-5)*1e-5;
KyMtr = round(KyMtr/1e-5)*1e-5;

Kx = sort(unique(KxMtr));
Ky = sort(unique(KyMtr));

KxNum = length(Kx);
KyNum = length(Ky);

Data = zeros(KzNum,KyNum,KxNum); % the parameter order: Rx,Tx,Ry,Ty
%%==============================================
%%initialization
Data_FK = permute(Data_FK, [2 1 4 3 5]);
Data_FK = reshape(Data_FK, [Kxr_Num*Kyr_Num, Kxt_Num*Kyt_Num, freqNum]);

[KxtIndMtr, KytIndMtr] = meshgrid(1:Kxt_Num, 1:Kyt_Num);
[KxrIndMtr, KyrIndMtr] = meshgrid(1:Kxr_Num, 1:Kyr_Num);

for mm = 1 : Kxt_Num*Kyt_Num
    indX_t = KxtIndMtr(mm);
    indY_t = KytIndMtr(mm);
    
    Kz_t = squeeze(KzMtr_t(:, indY_t, indX_t));
    for nn = 1 : Kxr_Num*Kyr_Num
        indX_r = KxrIndMtr(nn);
        indY_r = KyrIndMtr(nn);
        
        Kz_r = squeeze(KzMtr_r(:, indY_r, indX_r));

        Kz_interp = Kz_t + Kz_r;
        
        I_1 = Kz_interp>0;
        if sum(I_1)>0
            I_2 = Kz <= max(Kz_interp(I_1)) & Kz >= min(Kz_interp(I_1));
            
            IndKx = find(Kx == KxMtr(indX_r,indX_t));%(indX_r + indX_t - 1); % can be accelerated for a specific data array.
            IndKy = find(Ky == KyMtr(indY_r,indY_t));%(indY_r + indY_t - 1); % 
            
            Data(I_2,IndKy,IndKx) = Data(I_2,IndKy,IndKx) + interp1(Kz_interp(I_1), squeeze(Data_FK(nn, mm, I_1)), Kz(I_2), 'linear');
        end
    end
end
%%========
