function K = construct_Gaussian_kernel(n_gridx, n_gridy, cls, sig, coeff1)
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
    
    
    type=9; %3,8,9 all right
    %note most of these cases contain incorrect code from a long time ago,
    %ignore all below except case 9.
switch type
    case 1
        %%%%%% This case takes \xi in [0,1], partioned depending on cls,
        %%%%%% and exvaluates gaussian kernel exp( - | x - \xi_i )
     %====== function: K= a * exp{ - | x - \xi_i|^2  / 2sig^2} ===============%
    % ==== a = 1/(sqrt(2*pi)*sig)
    K = zeros(n_gridx*n_gridy,cls);
a = 1/(sqrt(2*pi)*sig);
xi = 0:1/(cls-1):1; % clsX1 vec between 0 and 1
for L = 1:cls
    hhh = exp( -( X - xi(L) ).^2/(2*sig^2) ).*exp( -( Y - xi(L) ).^2/(2*sig^2) ); %xi_i [0,1]
    
    K(:,L) = a*hhh;  
end
    case 2
        %%%% This evaluates thekernel exp{ -(X.^2 +Y.^2)/2sig^2 } on the
        %%%% discretized [0,1]x[0,1] grid
        n = n_gridx*n_gridy;
        
        
        repX = repmat(X,1,n);
        repY = repmat(Y,1,n);
        [repX,repY] = meshgrid(tx,ty);

        a = exp(-(repX.^2)/(2*sig^2));
        b = exp(-(repY.^2)/(2*sig^2));
        mEXP = a.*b;    
        
        coeff = 1/(sqrt(2*pi)*sig);
        K = coeff^2*mEXP;        
        
    case 3
        %%%% This evaluates thekernel exp{ -((X-X').^2 + (Y-Y').^2)/2sig^2 } on the
        %%%% discretized [0,1]x[0,1] grid
        n = n_gridx*n_gridy;
        
        [repX,repY] = meshgrid(tx,ty);
        repX = repmat(X,1,n);
        repY = repmat(Y,1,n);
        myX = repX - repX';
        myY = repY - repY';
        a = exp(-(myX.^2)/(2*sig^2));
        b = exp(-(myY.^2)/(2*sig^2));
        mEXP = a.*b;
        
        coeff = 1/(sqrt(2*pi)*sig);
        K = coeff^2*mEXP;  
        
    case 4
        %%%This evaluates exp( - (x^2 + y^2)) and repeats the 36x1 vec no
        %%%times
        no = 10;
        K = zeros(n_gridx*n_gridy,no);
        for L=1:no
            K(:,L) = exp(-(X.^2 + Y.^2)/(2*sig^2));
        end
        coeff = 1/(sqrt(2*pi)*sig);
        K = coeff*K;
        
    case 5
        %Consider only 6x6 grid
        X6 = reshape(X,[n_gridx,n_gridy]);
        Y6 = reshape(Y,[n_gridx,n_gridy]);
        
        a = exp(-(X6.^2)/(2*sig^2));
        b = exp(-(Y6.^2)/(2*sig^2));
        mEXP = a.*b;
        
        coeff = 1/(sqrt(2*pi)*sig);
        K = coeff*mEXP; 
        
    case 6
        %We look at T_ij = K(x_i,y_j) ~ exp( (x_i - y_j)^2 / 2sigma^2 );
        %Use our X and Y 36x1 vecs
        n= n_gridx*n_gridy;
        repX = repmat(X,1,n);
        repY = repmat(Y,1,n);
        
        mExp = exp( ( - (repX-repY).^2) / (2*sig^2) );
        coeff = 1/(sqrt(2*pi)*sig);
        
        K = coeff*mExp;
        
    case 7
        %We look at T_ij = K(x_i,y_j) ~ exp( (x_i - x_j)^2 / 2sigma^2 ) exp( (y_i - y_j)^2 / 2sigma^2 )
        n= n_gridx*n_gridy;
        %K will be 36x36
        for L=1:n
        Kx(L,:) = exp( - (( X(L) - X(:)).^2)/(2*sig^2));
        Ky(L,:) = exp( - (( Y(L) - Y(:)).^2)/(2*sig^2));
        end
        
        mExp = Kx.*Ky;
        coeff = 1/(sqrt(2*pi)*sig);
        coeff = coeff^2;
        K = coeff.*mExp;
        
    case 8
        %We look at now considering advice after 36x2
        XX = [X Y];
        YY = [X Y];
        
        n= n_gridx*n_gridy;
        %K will be 36x36
        %K(i,j) = exp( - ( (XX(i)-YY(i)).^2)/(2*sig^2)) * j part
        %T_12 = exp( - ( (XX(1,1)-YY(2,1
        
        for i=1:n
        for j=1:n
           a = XX(i,1) - YY(j,1);
           b = XX(i,2) - YY(j,2);
           tao = a.^2 + b.^2;
           K(i,j) = exp(-tao./(2*sig^2));
        end
        end
        
%         for i=1:n
%         for j=1:n
%         Kx(i,j) = exp( - (( XX(i,1) - YY(i,1)).^2)/(2*sig^2));
%         Ky(i,j) = exp( - (( XX(j,2) - YY(j,2)).^2)/(2*sig^2));
%         end
%         end
        
        %mExp = Kx.*Ky;
        coeff = 1/(sqrt(2*pi)*sig);
        coeff = coeff^2;
        K = coeff.*K;
        
    case 9
        n = n_gridx*n_gridy;
        
        [repX,repY] = meshgrid(tx,ty);
        repX = repmat(X,1,n);
        repY = repmat(Y,1,n);
        a = repX - repX';
        b = repY - repY';
        tao = sqrt(a.^2 + b.^2);
        tao(find(tao==0)) = eps;
        K = exp(-tao.^2./(2*sig^2));
        if coeff1 == 1
            coeff = 1;
        else
        coeff = 1/(sqrt(2*pi)*sig);
        end
        coeff = coeff^2;
        K = coeff.*K;
end
    

end