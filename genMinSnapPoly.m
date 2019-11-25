function coeff = genMinSnapPoly(w, n, k)
m = size(w);
m = m(1) - 1;
t0 = w(1,2);
tm = w(end,2);

dim = (n+1)*m;
Aeq = zeros(k*(m+1) + 2*m, dim);
beq = zeros(k*(m+1) + 2*m, 1);

% set boundary conditions for differentials less than equal to k (=0)
for z = 1:m
    for x = 1:k
        for y = 1:n-x+1
            if x == k
                continue;
            else
                Aeq(k*(z-1)+x, (z-1)*(n+1) + y) =  (factorial(n-y+1)/factorial(n-x-y+1))*(w(z,2)^(n-y-x+1));
                Aeq(k*(z-1)+k+x, (z-1)*(n+1) + y) =  -(factorial(n-y+1)/factorial(n-x-y+1))*(w(z+1,2)^(n-y-x+1));
            end
        end
    end
end

% set boundary condition for each waypoint
for x = 1:m
    for y = 1:n+1
        Aeq(k*(m+1) + 2*x - 1, y + (x-1)*(n+1)) = w(x, 2)^(n-y+1);
        Aeq(k*(m+1) + 2*x, y + (x-1)*(n+1)) = w(x+1, 2)^(n-y+1);
    end
    beq(k*(m+1) + 2*x-1, 1) = w(x,1);
    beq(k*(m+1) + 2*x, 1) = w(x+1,1);
end
% 
Aeq
beq

% H = zeros(dim, dim);
% for i = 1:m
%     H()
H = zeros(dim,dim);
f = zeros(1, dim);

for i = 1:m
    H((i-1)*(n+1)+1,(i-1)*(n+1)+1) = 1;

options = optimoptions('quadprog','Display','iter','MaxIterations',2000);
x0 = zeros(dim,1); 
coeff = quadprog(H,f,[],[],Aeq,beq,[],[],x0,options);
   
end