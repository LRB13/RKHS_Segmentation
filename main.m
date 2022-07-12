%%
load im2_40


%%
gamma = 1e-6; %reg of
iota = 1000;
lambda = 1e-6;
mu = 1e-3;
rho2 = 1e-9;


dProxConst = 1e-9;
betaProxConst = 10;
uProxConst = 4.0e-6;
c1ProxConst = 1e-9;
c2ProxConst = 1e-9;

[imNew,imT,imP,u] = rk_seg_convergencefinal_slower(im,mask,mu,lambda,gamma,rho2,dProxConst,betaProxConst,uProxConst,c1ProxConst,c2ProxConst);

%% figure and save
threshold = 0.5;
FigH = figure('Position', get(0, 'Screensize'));
imagesc(imNew); colormap gray; axis off; axis image; %title("uPc = " + uProxConst + ", lambda = " + lambda);
hold on; contour(u,[threshold,threshold],'r','LineWidth',2);
saveas(gcf,'output.png');