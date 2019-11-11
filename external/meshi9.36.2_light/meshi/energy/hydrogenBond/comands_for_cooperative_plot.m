Y1 = BHB;
Y2 = BHB;
Z = Y1'*Y2;
Z1 = -1*Z;
mesh(X1,X1,Z1)
[X1,Y1] = BHB
xline = linespace(min(X1),max(X1),33)
xline = linspace(min(X1),max(X1),33)
[X,Y] = meshgrid(xline,xline)
Ans1 = griddata(X1,X1,Z,X,Y,'cubic');
mesh(X,Y,Ans1)
axis tight
view([133 45])