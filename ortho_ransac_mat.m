function [F,nrinliers,inliers] = ortho_ransac_mat(x,y,bnd,iter)

if nargin<4,
    iter = 1000;
end
if nargin<3,
    bnd = 1e-3;
end

d = 3;
nn = size(x,2);

indy = zeros(d,iter);
for iii = 1:iter,
    indy(:,iii) = randperm(nn,d)';
end

xx = [x(1:2,:);y(1:2,:);ones(1,nn)];

x1 = x(:,indy(1,:));
x2 = x(:,indy(2,:));
x3 = x(:,indy(3,:));
y1 = y(:,indy(1,:));
y2 = y(:,indy(2,:));
y3 = y(:,indy(3,:));

[E1,E2]=minorthoF_mat(x1,x2,x3,y1,y2,y3);
%keyboard
E = [E1 E2];
ins = abs(E'*xx)<=bnd;
sumo = sum(ins,2);
[nrinliers,id]=max(sumo);
F = E(:,id);
inliers = ins(id,:);
