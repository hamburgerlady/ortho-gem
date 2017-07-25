% testscript for testing the rep ortho-gem

%% setup data
disp('Setting up data')
clear
load dino_withk
id1 = 10;
id2 = 14;
u1 = xx{id1};
u2 = xx{id2};

corrs = isfinite(u1(1,:)) & isfinite(u2(1,:));
u1 = u1(:,corrs);
u2 = u2(:,corrs);
% simulate gross outliers
nout = 20;
u1 = [u1 rand(2,nout)*500];
u2 = [u2 rand(2,nout)*500];
nn = size(u1,2);

% plot correspondences
figure(1)
clf
hold on
for iii = 1:nn,
    ll = plot([u1(1,iii) u2(1,iii)+500],[u1(2,iii) u2(2,iii)],'-');
    set(ll,'LineWidth',2);
end
title('Initial correspondences with simulated outliers');

% normalize
u1n = K\[u1;ones(1,nn)];
u2n = K\[u2;ones(1,nn)];
u1n = u1n(1:2,:);
u2n = u2n(1:2,:);



%% test minimal solver using RANSAC
disp('Testing minimal solver using RANSAC')

bnd = 1e-3;
nriter = 1000;
[F,nrinliers,inliers1] = ortho_ransac_mat(u1n,u2n,bnd,nriter);
disp(['Number of inliers after RANSAC: ' num2str(nrinliers)])

% plot inliers
figure(2)
clf
hold on
for iii = find(inliers1),
    ll = plot([u1(1,iii) u2(1,iii)+500],[u1(2,iii) u2(2,iii)],'-');
    set(ll,'LineWidth',2);
end
title('Inlier correspondences using RANSAC');



%% test optimal inlier solver
disp('Testing optimal inlier solver')

bnd = 1e-3;
[F,nrinliers,inliers2]=ortho_optimal(u1n,u2n,bnd);
disp(['Number of inliers after optimal estimation: ' num2str(nrinliers)])

% plot inliers
figure(3)
clf
hold on
for iii = find(inliers2),
    ll = plot([u1(1,iii) u2(1,iii)+500],[u1(2,iii) u2(2,iii)],'-');
    set(ll,'LineWidth',2);
end
title('Inlier correspondences maximizing number of inliers');



%% test least squares solver

disp('Testing least squares solver')

% using RANSAC inlier set 
[F1,rmin1] = mlorthoF_cm(u1n(:,inliers1),u2n(:,inliers1));
disp('Least squares Essential matrix from RANSAC:')
disp(F1)
disp('Gives residual:');
disp(rmin1)

% using RANSAC inlier set 
[F2,rmin2] = mlorthoF_cm(u1n(:,inliers2),u2n(:,inliers2));
disp('Least squares Essential matrix from optimal inlier set:')
disp(F2)
disp('Gives residual:');
disp(rmin2)




