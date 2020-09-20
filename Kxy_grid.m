function Kxy = Kxy_grid(Kxy_0, Num_Kxy, dKxy)
%
%   Kxy_0 --- the zero or middle point of the Kxy samples;
%   Num_Kxy --- the number of Kxy samples;
%   dKxy --- the sampling interval of Kxy.
%
% $Created by J.Wang @ TUDelft, Dec 5, 2018
%

if mod(Num_Kxy,2)==0
    Kxy = Kxy_0 + dKxy * ( -Num_Kxy/2 : Num_Kxy/2-1 );
else
    Kxy = Kxy_0 + dKxy * ( -(Num_Kxy-1)/2 : (Num_Kxy-1)/2 );
end