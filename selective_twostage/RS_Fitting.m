function [ f1 ] = RS_Fitting( z, lambda3, c1 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

        %%% normalise c1 and z
        c1 = (c1-min(z(:)))/(max(z(:))-min(z(:))); 
        z = (z-min(z(:)))/(max(z(:))-min(z(:)));
        
        N = 3;
        K = [0,multithresh(z(:),N-1),1];
        
        % LOWER THRESHOLD
        L_vect = max(c1 - K,0);
        L_vect2 = L_vect;
        TF_vect = 1:size(L_vect,2);
        TF_vect(L_vect==0) = [];
        L_vect2(L_vect==0) = [];
        TF2 = find(L_vect2==min(L_vect2));
        L = K(TF_vect(TF2(1)));
        
        % UPPER THRESHOLD
        H_vect = max(K - c1,0);
        H_vect2 = H_vect;
        TF_vect = 1:size(L_vect,2);
        TF_vect(H_vect==0) = [];
        H_vect2(H_vect==0) = [];
        if isempty(H_vect2)
            H = K(end);
        else
            TF2 = find(H_vect2==min(H_vect2));
            H = K(TF_vect(TF2(1)));
        end
        
        gamma1 = c1 - L;
        gamma2 = H - c1;
        
        TF3 = (z>=c1-gamma1) & (z<=c1);
        TF4 = (z<=c1+gamma2) & (z>c1);
        
        f3 = ( 1+((z-c1)/gamma1) ).*TF3 +...
            ( 1-((z-c1)/gamma2) ).*TF4 ;
        
        f1 = ((z-c1).^2) - lambda3 * f3;

end

