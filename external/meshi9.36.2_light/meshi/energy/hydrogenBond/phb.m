function[out1,out2,out3] = PHB(x,y)

%PHB computs the angle-punisg energy (vr 1.60) of an HB as a function of the
% distance between the hydrogen and the oxgen  atoms and the relavent angle
% 
%[dis,angle,energy] = BHB(dis) also return dis and angle.
% x (or dis) is the range of the distance. if no input x is given then x=
% 0:0.01:Dcut+0.5
% y (or angle) is the range of the angle. if no input y is given then
% y=90:PI/20:180
%
%example: plot(PHB)
%

Dcut = 4;
Acut = 150;
c1 = 10;
a1 = c1/Dcut^2;
b1 = -2*c1/Dcut;
y0Dis = 10/Dcut^2;
y0Ang = 40/Acut^2;


if (nargin == 0)
    x = 0:0.01:Dcut+0.5;
    y= 0:1:180;
end

for i = 1:length(x)
    for j=1:length(y)
        if(x(i) > Dcut | y(j) > Acut)
            energy(i,j) = 0;
        else
            %energy(i) = a1*x(i)^2 + b1*x(i) + c1;
            energy(i,j) = y0Dis*y0Ang*((x(i) - Dcut)^2)*((y(j) - Acut)^2);
        end
    end
end


if  nargout == 3,
    out1 = x; out2 = y; out3 = energy;
else
   if nargout == 2, 
        out1 = x;out2 = energy;
   else
       out1 = energy;
   end
end