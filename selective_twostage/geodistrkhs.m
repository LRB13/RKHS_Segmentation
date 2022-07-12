function [gd,mask,cols,rows] = geodistrkhs(Im,imP,seg_result,cols,rows)
%Inputs:
%%% Im - image
%%% imP - imP from RKHS
%%% seg_result - u result from an initial segmentation
%%% cols, rows - coords to input a mask


%%

[n,m] = size(Im);

% seg_result(seg_result>0.5) = 1;
% seg_result(seg_result<=0.5)=0;
% seg_result = extendregion(seg_result,2);
if nargin > 2
    seg_result(seg_result>.5) = 1;
    seg_result(seg_result<=.5) = 0;
end

if nargin == 3 
    mask = seg_result;
       [cols,rows] = find_boundary_coords(mask);
end
if nargin ==5 || nargin ==4
    mask = roipoly(Im,cols,rows);
    mask = mask+seg_result;
    mask(mask>1)=1;
end
    
    
 if nargin < 3
     figure; set(gcf,'units','normalized','outerposition',[0 0 1 1]);
 [mask,cols,rows] = roipoly(Im);
% else
%     mask = roipoly(Im,cols,rows);
 end

if size(cols)==1
    col = cols; row = rows;
    cols = []; rows = [];
    cols = [col col col+1];
    rows = [row row+1 row];
    mask = roipoly(Im,cols,rows);
end

%%% get gradient
ims = imgaussfilt(Im,2);
[gx,gy] = gradient(ims);
grad = gx.^2 + gy.^2;
%%%

%eucdist = timesweep(ones(n,m),mask);

beta = 1;
epsi = 1e-3;
theta= 0.05;
theta=0;

%imPs = imtgvsmooth(imP,0.01,0.01,10);
imPs = imP; %imPsc = imP;
imPsc = imP./max(abs(imP(:)));
f2 = epsi + 1000.*abs(imPsc.^2); %+ theta.*eucdist; %this is great!
tic
D = timesweep(f2,mask);
toc

gd=D./max(D(:));
end



function newu = extendregion(u,n)
%given an image u (heaviside fn essentialy), move it up down left and right
%note i'm imagining coordinate axis are as follows
% 5|
% 4|
%y |
% 2|
% 1|
% 0 - - - - - - - - -
%  0 1 2 3 4 5 6 ...
%            x
% so up, down, left, right aren't exactly what I imagine, but it doesn't
% matter
%%%%%5%% output is a slightly extended heaviside fn
for k=1:n
invu = 1-u;
[n,m] = size(u);
up = ones(n,m); down = up; left = up; right = up;

up(:,2:end) = invu(:,1:end-1);
down(:,1:end-1) = invu(:,2:end);
left(1:end-1,:) = invu(2:end,:);
right(2:end,:) = invu(1:end-1,:);

newinv = up.*down.*right.*left;
newu = 1-newinv;
u=newu;
end
end