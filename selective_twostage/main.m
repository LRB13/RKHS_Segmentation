

%% Image
load IM87

%% RKHS
nu=0.0001; %smoothing term (optional)
gamma = 1e-6; %regulariser
alpha = 1; %sparsity term
p = 4; %patchsize
o = 3; %overlap
[imNew,imT,imP] = model_rk(im,nu,gamma,alpha,p,o);

%% Geodesic distance
mask = roipoly(im); %region of interest
[gd,~] = geodistrkhs(im,imP,mask);

%% Segmentation
g = 1./(1+1000.*imP.^2);
mu=1; %length term in segmentation
lambda = 1; %intensity term in segmentation
theta = 80; %geodesic distance contsraint
flag = 2; %%%% flag determines which data fidelity term is used. flag = 1 is Chan-Vese, i.e. (z-c_1)^2 - (z-c_2)^2
        %%%% flag = 2 is a reformulated Chan-Vese  (usually performs better
        %%%% in general on medical data) https://drive.google.com/file/d/1Gz2gbhBgP89M0ca0kUP88s6b0MTRYMXW/view

 u = seg(im,lambda,theta,mu,mask,gd,g,flag);
 
 

 