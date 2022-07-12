function [imNew,imT,imP] = model_rk(im,nu,lambda,alpha,patchsize,overlap)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[n,m] = size(im);

%patchsize = 4; overlap = 2;

%%%
sig=12;
K = construct_Gaussian_kernel(patchsize,patchsize,1,sig,0);
[P,~] = consP_2d_arbDim(patchsize,patchsize,patchsize,patchsize,1e-4,12);

 %%% Parameters
 rho1=2; rho2=2; rho3=2;
 %alpha=1; 
 %lambda = 1e-6; 
 mu = 1; %fitting parameter
 %nu = 1e1; %parameter in new smoothing term
 iota = 1000; %parameter in g, edge detector
 
param.r1 = rho1;
param.r2 = rho2;
param.r3 = rho3;
param.alpha = alpha;
param.lambda = lambda;
param.mu = mu;
param.nu = nu;
param.iota = iota;

 [mm2, nn2] = size(P);

%%%%
Matrix.C1 = rho2*(K'*K) + 2*lambda*K;
Matrix.C2 = rho2*K'*P;
Matrix.P1 = rho2*P'*K;
Matrix.P2 = rho2*(P'*P) + rho1*eye(nn2); 
Matrix.P2I = inv(Matrix.P2);
Matrix.CA = inv(Matrix.C1-Matrix.C2*inv(Matrix.P2)*Matrix.P1);


h=n; w=m;
%%% grid for patch 
gridy = 1:patchsize - overlap: w-(mod(w,patchsize-overlap)+1+patchsize-overlap);
gridy = setdiff(gridy, [(w-patchsize+1):w]);  % delete rest elements
gridx = 1:patchsize - overlap: h-(mod(h,patchsize-overlap)+1+patchsize-overlap);
gridx = setdiff(gridx, [(h-patchsize+1):h]); % delete rest elements
gridx = [gridx h-patchsize+1]; gridy = [gridy w-patchsize+1];

counter = 0;

A = zeros(h,w);
B = zeros(h,w);
Pb = zeros(h,w);
Kd = zeros(h,w);


tic
for i=1:length(gridx)
for j=1:length(gridy)
    counter = counter+1;
    xx = gridx(i); yy = gridy(j);
    imPatch = im(xx:xx+patchsize-1,yy:yy+patchsize-1);
    %krPatch = K(xx:xx+patchsize-1,yy:yy+patchsize-1,:); %old kernel
    if norm(imPatch,2) > 0
    %    keyboard
    end
    %%%
    
    if norm(imPatch,2) == 0 %dont compute if patch is all zeroes
        c = zeros(size(K,1),1);
        beta = zeros(nn2,1);
    else
    [c,beta] = iterate(imPatch,K,P,param,Matrix);
    end
    a1 = K*c;
    a2 = P*beta;
    
    a1 = reshape(a1,[patchsize,patchsize]);
    a2 = reshape(a2,[patchsize,patchsize]);
    newPatch = a1 + a2;
    
    A(xx:xx+patchsize-1, yy:yy+patchsize-1) = A(xx:xx+patchsize-1, yy:yy+patchsize-1) + newPatch;
    B(xx:xx+patchsize-1, yy:yy+patchsize-1) = B(xx:xx+patchsize-1, yy:yy+patchsize-1) + 1; %OVERLAP
    
    %%%%highres testing
%     XX = Gridx(i); YY = Gridy(j);
%     B_H(XX:XX+size_h-1, YY:YY+size_h-1) = B_H(XX:XX+size_h-1, YY:YY+size_h-1) + 1;
%     Z(XX:XX+size_h-1, YY:YY+size_h-1) = Z(XX:XX+size_h-1, YY:YY+size_h-1) + b2;
    
    Kd(xx:xx+patchsize-1, yy:yy+patchsize-1) = Kd(xx:xx+patchsize-1, yy:yy+patchsize-1) + a1;
    Pb(xx:xx+patchsize-1, yy:yy+patchsize-1) = Pb(xx:xx+patchsize-1, yy:yy+patchsize-1) + a2;
    %if j==length(gridy)
    if i==50 && j==50
      % keyboard
    end
end
end
toc

B(B==0) = 1;
imNew = A./B;

imT = Kd./B;
imP = Pb./B;


end

function [c,beta] = iterate(im,K,P,param,Matrix)

[n,m] = size(im);
%im = im(:);

rho1 = param.r1;
rho2 = param.r2;
rho3 = param.r3;
alpha = param.alpha;
lambda = param.lambda;
mu = param.mu;
nu = param.nu;
iota = param.iota;

z = zeros(n,m);
b2 = z;
v1 = zeros(n,m);
v2 = zeros(n,m);
b3x = v1;
b3y = v2;

[m3,n3] = size(P);
beta = zeros(n3,1);
b1 = zeros(n3,1); 
u = b1;

C1 = Matrix.C1;
C2 = Matrix.C2;
P1 = Matrix.P1;
P2 = Matrix.P2;
P2I = Matrix.P2I;
CA = Matrix.CA;

   lambda1 = psf2otf([1,-1],[n,m]); lambda2 = psf2otf([1;-1],[n,m]);
   lambda1Conj = conj(lambda1); lambda2Conj = conj(lambda2);

trackc=[];
trackbeta=[];
E=[];

for t = 1:10
    
    zvec=z(:); b2vec = b2(:);
    %%% c beta problem
   CR = rho2*K'*(zvec+b2vec);
   PR = rho2*P'*(zvec+b2vec) + rho1*(u+b1);
   
   c = CA*(CR - C2*P2I*PR);
   beta = P2I*(PR-P1*c);
   Kc = K*c; Kc = reshape(Kc,[n,m]);
   Pb = P*beta; Pb = reshape(Pb,[n,m]);
   g = 1./(1+iota.*(Pb).^2); g = g.*nu;
   
   %%% u problem
   p1 = sign(beta - b1);
   p2 = abs(beta - b1);
   u = p1.*max(p2 - alpha/rho1, 0);
   b1 = b1 + (u - beta);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
   %%% v problem
   [v10,v20] = gradient(z);
   v1a = sign(v10-b3x);
   v2a = sign(v20-b3y);
   v1b = abs(v10-b3x);
   v2b = abs(v20-b3y);
   v1 = v1a.*max(v1b - g./rho3,0);
   v2 = v2a.*max(v2b - g./rho3,0);
   b3x = b3x + v1 - v10;
   b3y = b3y + v2 - v20;
   %%%%%%%%%%%

   
   %%% z problem   
   z0a = rho3*(lambda1.*lambda1Conj + lambda2.*lambda2Conj);
   z0b = 2*mu+rho2; zL = z0a+z0b;
    
   
   b2 = b2 + z - (Kc + Pb);
   
   zR0 = fft2(2*mu*im + rho2*(Kc+Pb-b2)) + rho3.*(lambda1Conj.*fft2(v1+b3x) + lambda2Conj.*fft2(v2+b3y));
   z0 = zR0./zL;
   z = real(ifft2(z0));
   
   trackc=[trackc;c(1)];
   trackbeta=[trackbeta;beta];
   
   
e1 = sqrt(v1.^2 + v2.^2); e1 = sum(e1(:));
e2 = mu.*(im-z).^2; e2 = sum(e2(:));
e3 = lambda*c'*K*c;
e4 = alpha*norm(u,1);
e5 = (rho1/2)*norm(u-beta+b1,2); e5 = e5.^2;
e6 = (rho2/2)*norm(z-(Kc+Pb)+b2,2); e6 = e6.^2;
[zx,zy] = gradient(z);
e7 = (rho3/2)*norm(v1-zx+b3x,2); e7 = e7.^2;
e8 = (rho3/2)*norm(v2-zy+b3y,2); e8 = e8.^2;
energy = e1 + e2 + e3 + e4 + e5 + e6 + e7 + e8;
E = [E;energy];
end




%keyboard

end