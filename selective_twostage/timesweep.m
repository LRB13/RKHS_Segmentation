function D = timesweep(f,mask)

[n,m] = size(f);

%maxit = 1200;

inf = 1e6;
T = 1000;
u = inf*(1-mask);

unew = zeros(n,m);

res0 = []; stop = .00005; R = 10*stop;

count = 0;
xx = 1:n; yy = 1:m;
orderx = [xx;fliplr(xx);fliplr(xx);xx];
ordery = [yy;yy;fliplr(yy);fliplr(yy)];

h=0.005;

while(R>stop)
    count = count+1;
    order = mod(count-1,4)+1;
    x = orderx(order,:);
    y = ordery(order,:);
    
    oldu = u;
    for i=x
        for j=y
            
            %%% set a and b
            if i==1
                a = u(2,j);
            elseif i==n
                a = u(n-1,j);
            else
                a = min( u(i-1,j), u(i+1,j));
            end
            
            if j==1
                b = u(i,2);
            elseif j==m
                b = u(i,m-1);
            else
                b = min( u(i,j-1), u(i,j+1));
            end
            
            if min(a,b) < T
                cond = abs(a-b);
                if cond >= f(i,j)*h
                    ubar = min(a,b) + f(i,j)*h;
                else
                    fh = f(i,j)*h;
                    ubar = (1/2)*(a+b+ sqrt( (2*fh^2) - (a-b)^2 ));
                end
                u(i,j) = min(u(i,j),ubar);
            end
            
            
        end
    end
       % imagesc(unew); drawnow
        R = norm(oldu-u,2);
    res0 = [res0 R];
    %u = unew;
    %imagesc(u); drawnow
    
end


D = (u-min(u(:)))./(max(u(:))-min(u(:)));
disp("iter = " + count);

end