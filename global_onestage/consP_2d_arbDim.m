function [H, HH] = consP_2d_arbDim(n_gridx, n_gridy, N_gridx, N_gridy, ksi, m)

   % ksi: control the amount of smooth
   % larger ksi, more smoothing
   

 
   %----------coarse grid--------%
    tx = [0:(n_gridx-1)]'/(n_gridx-1);
    ty = [0:(n_gridy-1)]'/(n_gridy-1);
%     tx = [1:n_gridx]'/n_gridx;
%     ty = [1:n_gridy]'/n_gridy;
    
    X = repmat(tx,n_gridy,1);
    Y = [];

    for i = 1:n_gridy
        e = repmat(ty(i), n_gridx,1);
        Y = [Y; e];
    end

    %------fine grid--------------
    ttx = [0:(N_gridx-1)]'/(N_gridx-1);
    tty = [0:(N_gridy-1)]'/(N_gridy-1);
%     ttx = [1:N_gridx]'/N_gridx;
%     tty = [1:N_gridy]'/N_gridy;
    
    XX = repmat(ttx,N_gridy,1);
    YY = [];

    for i = 1:N_gridy
        ee = repmat(tty(i), N_gridx,1);
        YY = [YY; ee];
    end
    
    %====== function: H=1/2+1/pi*arctan(z/ksi) ===============%
    %====== z = cos(theta)*x + sin(theta)*y + c==============%
    %====== c = tx or ttx  =====================%
     n = n_gridx*n_gridy;
     N = N_gridx*N_gridy;
     %ksi = 1e-4;
    %m =4;
    theta = [0: 2*pi/m: 2*pi];
    %theta = [0: pi/m: pi];
    %====== construct H ===========%
%     grid = 19;
%     c = [0:grid]/grid; 
%     c = X;
%     grid = size(c,1);
    %grid = 25;
    grid = n;
    c = [0:(grid-1)]'/(grid-1);
    iter = 1;
    for j = 1:m
       for i = 1:grid
         z = cos(theta(j))*X + sin(theta(j))*Y + repmat(c(i),n,1);
         z = z/ksi;
         r = atan(z);
         H(:,iter) = repmat(1/2, n, 1) + 1/pi*r;
         iter = iter + 1;
       end
    end
    
   %====== construct HH =========%
    iter = 1;
    d = ttx;
      for j = 1:m
       for i = 1:grid
         z = cos(theta(j))*XX + sin(theta(j))*YY + repmat(c(i),N,1);
         z = z/ksi;
         r = atan(z);
         HH(:,iter) = repmat(1/2, N, 1) + 1/pi*r;
         iter = iter + 1;
       end
    end 
    
    
    
    
    
    
    
