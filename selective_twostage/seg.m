    %%% take in image, g
function u = seg(Im,lambda,theta,mu,u0,gd,g,flag)
%%%% This function gives u in [0,1]
%%%% Model: mu \int g | \grad u | + \lambda \int ( (z-c1)^2 - (z-c2)^2 ) v +
%%%% 1/2\theta \int (u-v)^2
%%%% BRESSON CHAMBOLLE 2007

if nargin < 6
    [gd,tempu0] = geodist(Im);
    if nargin < 4
        u0 = tempu0;
    end
end
    

if nargin < 7
    ims = imgaussfilt(Im,1.5);
    [grad1,grad2] = gradient(ims);
    grad = grad1.^2 + grad2.^2;
    g = 1./(1+1000.*grad);
end

if nargin < 8
    flag = 1;
end

z = Im;

[h,w] = size(Im);

%%%% Optional

%Params
%lambda = 1; %fitting term parameter
rho = 1; %term on 1/{2 \theta} \int (u-v)^2
maxit = 250;
dt = 1/8;
%stop = 1e-4;
stop=0.001;
%c1 = 0.8;
%c2 = 0.2;

%%%%%%%
sigma = 5;%size of kernel
Ksigma=fspecial('gaussian',round(2*sigma)*2 + 1,sigma); %  kernel
%%%%%%%%


Pd = theta.*gd;

p1 = zeros(h,w); p2 = zeros(h,w);
u = zeros(h,w); v = zeros(h,w); u=u0; v=u0;
res0 = [];
for k=1:maxit

uold = u; vold = v;

c1 = sum(sum(v.*z))/sum(sum(v));
c2 = sum(sum((1-v).*z))/sum(sum((1-v)));

r0 = (Im - c1).^2 - (Im - c2).^2;

if flag == 1
r1 = lambda.*r0;
elseif flag == 2
r1 = RS_Fitting(z,lambda,c1);
end

r = r1 + Pd;

[p1,p2,divP] = Pstuff(p1,p2,v,g,rho*mu,dt);

    
u = v - rho.*mu.*divP;

v = max(u-rho*r,0);
v = min(v,1);

if mod(k,10)==0
imagesc(Im); colormap gray; title("k = " + k); hold on; contour(u,[0.5,0.5],'r','LineWidth',2); drawnow
end
resu = norm(u - uold); resv = norm(v - vold);
Res = max([resu; resv]); res0 = [res0, Res];
if Res < stop; break; end
%%% END Chambolle %%%

if k==100
    %keyboard
end
    
end
%figure; plot(res0); title("Residual");

end

function [newP1,newP2,divP] = Pstuff(p1,p2,v,g,theta,dt)

divP = divp(p1,p2);
n0 = divP - v./theta; 
[n1,n2] = dualgrad(n0); %grad of divp - v/theta
D0 = sqrt(n1.^2 + n2.^2); %rhs of denominator
D = 1 + (dt./g).*D0; % 1+ dt/g | \nabla ( divp - v/theta) | i.e. the denominator
 
a1 = p1 + dt.*n1;
a2 = p2 + dt.*n2;

newP1 = a1./D;
newP2 = a2./D;
end

function [gradp1,gradp2] = dualgrad(a)
[h,w] = size(a);
gradp1 = zeros(h,w); gradp2 = zeros(h,w);

gradp1(1:h-1,:) = a(2:h,:) - a(1:h-1,:);
gradp2(:,1:w-1) = a(:,2:w) - a(:,1:w-1);
end

function divP = divp(p1,p2)
[h,w] = size(p1);
divP = zeros(h,w);
divP(1:h-1,:) = divP(1:h-1,:) + p1(1:h-1,:);
divP(2:h,:) = divP(2:h,:) - p1(1:h-1,:);
divP(:,1:w-1) = divP(:,1:w-1) + p2(:,1:w-1);
divP(:,2:w) = divP(:,2:w) - p2(:,1:w-1);
end

