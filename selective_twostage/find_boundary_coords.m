function [cols,rows] = find_boundary_coords(u0)
[n,m] = size(u0);

u = zeros(n,m);
u(u0>0.5) = 1;

[gx,gy] = gradient(u);
grad = gx.^2 + gy.^2;

cols=[]; rows=[];
for i=1:n
    for j=1:m
        if grad(i,j) > 0
            cols = [cols;j];
            rows = [rows;i];
        end
    end
end



end