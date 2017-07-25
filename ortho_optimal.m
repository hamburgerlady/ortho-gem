function [F,nrinliers,inliers]=ortho_optimal(x,y,bnd,nrinliers_break)
% function [F,nrinliers,inliers]=ortho_optimal(x,y,bnd,nrinliers_break)
%
% Finds the orthographic essential matrix that optimizes the number of correspondences
% among a given set of tentative correspondences. 
% Input:
% x,y: 2xn matrices of tentative correspondences, x in view one and y in
% view two. 
% bnd: the error bound
% nrinliers_break: set to a minimum number of found inliers if you want to
% exit earlier than exhaustive enumeration.
% Output:
% F: The optimal essential matrix (vector of the 5 non-zero elements).
% inliers: 1xn vector indicating the inlier set.
% nrinliers: the optimal number of inliers (=sum(inliers)).


if nargin<4,
    nrinliers_break = inf;
end
x = x(1:2,:);
y = y(1:2,:);
nn = size(x,2);
eppe=1e-7;
nrinliers = 0;
allX = [x(1:2,:);y(1:2,:);ones(1,nn)]';
for iii = 1:nn-2,
    for jjj = iii+1:nn-1
        for kkk = jjj+1:nn,
            
            Fsols = boundorthoF(x(:,[iii jjj kkk]),y(:,[iii jjj kkk]),bnd);
            Fsols = real(Fsols(:,abs(imag(Fsols(1,:)))<eppe));
            
            ressy = abs(allX*Fsols)<=bnd+2*eps;
            [maxi,id]=max(sum(ressy));
            
            if ~isempty(Fsols),
                if maxi>nrinliers && sum(Fsols(1:4,id))~=0,
                    nrinliers = maxi;
                    inliers = ressy(:,id)';
                    F = Fsols(:,id);
                end
            end
            if nrinliers>nrinliers_break,
                break;
            end
        end
    end
end



function Fsols = boundorthoF(x,y,bnd)

d1x = x(:,2)-x(:,1);
d2x = x(:,3)-x(:,1);
d1y = y(:,2)-y(:,1);
d2y = y(:,3)-y(:,1);
Fsols = zeros(5,32);
count = 1;
s1 = [-1 -1 -1 -1 1 1 1 1];
s2 = [-1 -1 1 1 -1 -1 1 1];
s3 = [-1 1 -1 1 -1 1 -1 1];

for iii = 1:8,
    bd1 = bnd*(s2(iii)-s1(iii));
    bd2 = bnd*(s3(iii)-s1(iii));
    cc1 = [d1x' d1y' -bd1];
    cc2 = [d2x' d2y' -bd2];
    sols = mindifforthoF(cc1,cc2);
    esol = -sols(1,:)*x(1,1)-sols(2,:)*x(2,1)-sols(3,:)*y(1,1)-sols(4,:)*y(2,1)+bnd*s1(iii);
    Fsols(:,count:count+3)=[sols;esol];
    count = count+4;
end


function sols = mindifforthoF(c1,c2)

% eqny = k1*c^2+ k2*c*d + k3*c + k4*d^2 + k5*d + k6
k1 = (c1(1)^2*c2(3)^2 - 2*c1(1)*c1(3)*c2(1)*c2(3) + c1(2)^2*c2(3)^2 - 2*c1(2)*c1(3)*c2(2)*c2(3) + c1(3)^2*c2(1)^2 + c1(3)^2*c2(2)^2)/(c1(1)*c2(2) - c1(2)*c2(1))^2;
k2 = (2*(c1(3)*c1(4)*c2(1)^2 + c1(3)*c1(4)*c2(2)^2 + c1(1)^2*c2(3)*c2(4) + c1(2)^2*c2(3)*c2(4) - c1(1)*c1(3)*c2(1)*c2(4) - c1(1)*c1(4)*c2(1)*c2(3) - c1(2)*c1(3)*c2(2)*c2(4) - c1(2)*c1(4)*c2(2)*c2(3)))/(c1(1)*c2(2) - c1(2)*c2(1))^2;
k3 = (2*(c1(3)*c1(5)*c2(1)^2 + c1(3)*c1(5)*c2(2)^2 + c1(1)^2*c2(3)*c2(5) + c1(2)^2*c2(3)*c2(5) - c1(1)*c1(3)*c2(1)*c2(5) - c1(1)*c1(5)*c2(1)*c2(3) - c1(2)*c1(3)*c2(2)*c2(5) - c1(2)*c1(5)*c2(2)*c2(3)))/(c1(1)*c2(2) - c1(2)*c2(1))^2;
k4 = (c1(1)^2*c2(4)^2 - 2*c1(1)*c1(4)*c2(1)*c2(4) + c1(2)^2*c2(4)^2 - 2*c1(2)*c1(4)*c2(2)*c2(4) + c1(4)^2*c2(1)^2 + c1(4)^2*c2(2)^2)/(c1(1)*c2(2) - c1(2)*c2(1))^2;
k5 = (2*(c1(4)*c1(5)*c2(1)^2 + c1(4)*c1(5)*c2(2)^2 + c1(1)^2*c2(4)*c2(5) + c1(2)^2*c2(4)*c2(5) - c1(1)*c1(4)*c2(1)*c2(5) - c1(1)*c1(5)*c2(1)*c2(4) - c1(2)*c1(4)*c2(2)*c2(5) - c1(2)*c1(5)*c2(2)*c2(4)))/(c1(1)*c2(2) - c1(2)*c2(1))^2;
k6 = (- c1(1)^2*c2(2)^2 + c1(1)^2*c2(5)^2 + 2*c1(1)*c1(2)*c2(1)*c2(2) - 2*c1(1)*c1(5)*c2(1)*c2(5) - c1(2)^2*c2(1)^2 + c1(2)^2*c2(5)^2 - 2*c1(2)*c1(5)*c2(2)*c2(5) + c1(5)^2*c2(1)^2 + c1(5)^2*c2(2)^2)/(c1(1)*c2(2) - c1(2)*c2(1))^2;
% p1 = k2*c*d + k3*c + (k4-k1)*d^2 + k5*d + k6+k1

q1 = (k1*k2^2 - k1*k3^2 + k3^2*k4 + k2^2*k6 - k2*k3*k5)/k2^2;
q2 = -(k1^2 - 2*k1*k4 + k2^2 + k4^2)/k2;
q3 = (k3*k1^2 + 2*k5*k1*k2 - 2*k3*k1*k4 - k3*k2^2 - 2*k5*k2*k4 + k3*k4^2)/k2^2;
q4 = (k1^2*k2 - k2*k5^2 + k2^3 - k1*k2*k4 + k1*k2*k6 - k1*k3*k5 - k2*k4*k6 + k3*k4*k5)/k2^2;
q5 = -(k1^2*k3 - k2^2*k3 + k1*k2*k5 - k1*k3*k4 + k1*k3*k6 + k2*k5*k6 - k3*k4*k6)/k2^2;
% p2 =  q1*c + q2*d^3 + q3*d^2 + q4*d + q5

% p0 = c^2+d^2-1
% basis b = [1 c d d^2]
% action d
% d*b = [d cd d^2 d^3] = [b(3) p1 b(4) p2]

M = [0 0 1 0;-k1/k2-k6/k2  -k3/k2 -k5/k2 k1/k2 - k4/k2;0 0 0 1;-q5/q2 -q1/q2 -q4/q2 -q3/q2];
if all(isfinite(M(:))),
    [sol_cd,~] = eig(M);
    sol_c = sol_cd(2,:)./sol_cd(1,:);
    sol_d = sol_cd(3,:)./sol_cd(1,:);
    
    sol_a = (c1(2)*c2(5) - c1(5)*c2(2) + sol_c*c1(2)*c2(3) - sol_c*c1(3)*c2(2) + c1(2)*c2(4)*sol_d - c1(4)*c2(2)*sol_d)/(c1(1)*c2(2) - c1(2)*c2(1));
    sol_b = -(c1(1)*c2(5) - c1(5)*c2(1) + sol_c*c1(1)*c2(3) - sol_c*c1(3)*c2(1) + c1(1)*c2(4)*sol_d - c1(4)*c2(1)*sol_d)/(c1(1)*c2(2) - c1(2)*c2(1));
    
    sols = [sol_a;sol_b;sol_c;sol_d];
else
    sols = zeros(4,4);
end

