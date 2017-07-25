function [E1,E2] = minorthoF_mat(x1,x2,x3,y1,y2,y3)
% function [E1,E2] = minorthoF_mat(x1,x2,x3,y1,y2,y3)
%
% Estimates the orthographic essential matrix from three point correspondences
% Input x1,x2,x3 : 2xn matrices of points in view one, row 1 and 2 are x
% and y coordinate of the points respectively. The function handles multiple inputs so 
% columns correspond to different sets of 3-point correspondences.
% Input y1,y2,y3 : 2xn matrices of points in view two.
% Output E1 and E2 : The solution essential matrices, 5xn matrices (each 3-point
% correspondence gives two solutions).
% Rows are the values of the non-zero elements of the essential matrices,
% and columns correspond to the different sets of input 3-point
% correpondences.

xd1 = x2-x1;
yd1 = x3-x1;
xd2 = y2-y1;
yd2 = y3-y1;

denny = xd1(1,:).*yd1(2,:)-xd1(2,:).*yd1(1,:);

aac = (xd1(2,:).*yd2(1,:)-xd2(1,:).*yd1(2,:))./denny; 
aad = (xd1(2,:).*yd2(2,:)-xd2(2,:).*yd1(2,:))./denny;
bbc = (xd2(1,:).*yd1(1,:)-xd1(1,:).*yd2(1,:))./denny;
bbd = (xd2(2,:).*yd1(1,:)-xd1(1,:).*yd2(2,:))./denny;

dd_2 = - aac.^2 + aad.^2 - bbc.^2 + bbd.^2;
dd_1c = 2*aac.*aad + 2*bbc.*bbd;
dd_0 = aac.^2 + bbc.^2 - 1;

d4_4 = dd_1c.^2 + dd_2.^2;
d4_2 = - dd_1c.^2 + 2*dd_0.*dd_2;
d4_0 = dd_0.^2;


tmp = sqrt((d4_2.^2 - 4*d4_4.*d4_0));
d2sol_1 = -(d4_2+tmp)./d4_4/2;
d2sol_2 = -(d4_2-tmp)./d4_4/2;

dsol_1 = sqrt(d2sol_1);
dsol_2 = sqrt(d2sol_2);

csol_1 = -((-aac.^2+aad.^2-bbc.^2+bbd.^2).*dsol_1.^2+aac.^2+bbc.^2-1)./(2*aac.*aad.*dsol_1+2*bbc.*bbd.*dsol_1);
asol_1 = aac.*csol_1+aad.*dsol_1;
bsol_1 = bbc.*csol_1+bbd.*dsol_1;
esol_1 = -asol_1.*x1(1,:)-bsol_1.*x1(2,:)-csol_1.*y1(1,:)-dsol_1.*y1(2,:);

E1 = [asol_1; bsol_1; csol_1; dsol_1; esol_1];


csol_2 = -((-aac.^2+aad.^2-bbc.^2+bbd.^2).*dsol_2.^2+aac.^2+bbc.^2-1)./(2*aac.*aad.*dsol_2+2*bbc.*bbd.*dsol_2);
asol_2 = aac.*csol_2+aad.*dsol_2;
bsol_2 = bbc.*csol_2+bbd.*dsol_2;
esol_2 = -asol_2.*x1(1,:)-bsol_2.*x1(2,:)-csol_2.*y1(1,:)-dsol_2.*y1(2,:);

E2 = [asol_2; bsol_2; csol_2; dsol_2; esol_2];
















