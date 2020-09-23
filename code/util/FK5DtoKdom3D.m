function [Data,Kx,Ky,dKz]= FK5DtoKdom3D(Kr,Kxt,Kyt,Kxr,Kyr,Data_FK)
% FK5DtoKdom3D performs the 5-D to 3-D mapping.
% 
% Parameter:
%  Kr  ---  the wavenumber
%  Kxt ---  the Kx samples over transmit aperture
%  Kyt ---  the Ky samples over transmit aperture
%  Kxr ---  the Kx samples over receive aperture
%  Kyr ---  the Ky samples over receive aperture
%  Data_FK ---- the F-K domain data array
%
% Output:
%   Data  ---  the data array after 5-D to 3-D mapping
%   Kx    ---  the Kx samples after mapping
%   Ky    ---  the Ky samples after mapping
%   dKz   ---  sampling interval of Kz
%
% Author: J.Wang, Aug 25, 2015
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

dKz = 2*(Kr(2)-Kr(1));% (Kr(2)-Kr(1));

KzNum = floor((Kzmax - Kzmin)/dKz) + 1;   % could be chosen based on the system parameter for acceleration 
Kz = Kzmin + dKz * (0:KzNum-1);

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
Data_FK = permute(Data_FK, [2 1 5 4 3]);
Data_FK = reshape(Data_FK, [Kxr_Num*Kyr_Num, freqNum, Kxt_Num*Kyt_Num]);

KzMtr_t = reshape(KzMtr_t, [freqNum Kxt_Num*Kyt_Num]);
KzMtr_r = reshape(KzMtr_r, [freqNum Kxr_Num*Kyr_Num]);

[KxtIndMtr, KytIndMtr] = meshgrid(1:Kxt_Num, 1:Kyt_Num);
[KxrIndMtr, KyrIndMtr] = meshgrid(1:Kxr_Num, 1:Kyr_Num);

for mm = 1 : Kxt_Num * Kyt_Num
    Kz_t_rep = repmat(KzMtr_t(:,mm),[1  Kxr_Num*Kyr_Num]);
    Kztemp = (Kz_t_rep + KzMtr_r)/dKz;
    Kztemp(Kztemp>KzNum/2) = Kztemp(Kztemp>KzNum/2) - KzNum;

    mm_t = KxtIndMtr(mm);
    nn_t = KytIndMtr(mm);
    
    Data_FK_temp = Data_FK(:, :, mm);
    for kk = 1 : Kxr_Num * Kyr_Num
        indKzMtr = (Kztemp(:,kk)~=0);
        indKzMtrNum = sum(indKzMtr);
        
        if indKzMtrNum>0
            [fk, ~] = nufft1d1(indKzMtrNum,Kztemp(indKzMtr,kk)/KzNum * 2 * pi, Data_FK_temp(kk,indKzMtr), 1, 1e-6, KzNum);
        
            hh = KxrIndMtr(kk);
            rr = KyrIndMtr(kk);
            
            IndKx = find(Kx == KxMtr(hh,mm_t));%(hh+mm_t-1); %
            IndKy = find(Ky == KyMtr(rr,nn_t));%(rr+nn_t-1); % 
                    
            Data(:,IndKy,IndKx) = Data(:,IndKy,IndKx) + fk;
        end
    end
end

%%========











