function[out1,out2] = BHB(x)

% BHB comput the energy of a hydrogen bond (vr. 1.60) as a function of the
% distance between the hydrogen and the oxgen  atoms.
% 
%[dis,energy] = BHB(dis) also return dis.
% x (or dis) is the range of the distance. if no input x is given then x= 0:0.01:rMAx+1 
%
%example: plot(BHB)
%

si = 1.7;
epsi = 1;
rMax = 5.5;
Dmin = si*2^(1/6);

if nargin == 0, x = 0:0.01:rMax+1; end

for i = 1:length(x)
    if(x(i) > rMax)
        energy(i) = 0;
    else
       if(x(i) < Dmin)
          energy(i) = (0.01*epsi*(x(i) - Dmin)^2 - epsi);%0.9627);
        else
            energy(i) = (epsi*4*si^6*(si^6*x(i)^-12 - 1*x(i)^-6)) *((rMax - x(i))^2/((rMax - x(i))^2 + 0.00001));
        end
    end
end
if nargout == 2,
    out1 = x; out2 = energy;
else
    out1 = energy;
end