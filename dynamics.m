function dx = dynamics(t,y,xcoeff,ycoeff,zcoeff,psicoeff,time,params)

for i = 1:params.m
    if t>=time(i) && t <=time(i+1)
        index = i;
        break
    end
end

syms T real
x_des = 0;
y_des = 0;
z_des = 0;
psi_des = 0;
for j = 1:params.nx+1
    x_des = x_des + xcoeff((index-1)*(params.nx+1) + j)*(T^(params.nx+1-j));
end

for j = 1:params.ny+1
    y_des = y_des + ycoeff((index-1)*(params.ny+1) + j)*(T^(params.ny+1-j));
end

for j = 1:params.nz+1
    z_des = z_des + zcoeff((index-1)*(params.nz+1) + j)*(T^(params.nz+1-j));
end

for j = 1:params.npsi+1
    psi_des = psi_des + psicoeff((index-1)*(params.npsi+1) + j)*(T^(params.npsi+1-j));
end

sigma(T) = [x_des,y_des,z_des,psi_des];
dsigma(T) = diff(sigma, T);
d2sigma(T) = diff(dsigma, T);
d3sigma(T) = diff(d2sigma, T);

phi = y(4);
theta = y(5);
psi = y(6);

r = reshape(y(1:3),[3,1]);
dr = reshape(y(7:9),[3,1]);
w_b = reshape(y(10:12),[3,1]);

xW = [1;0;0];
yW = [0;1;0];
zW = [0;0;1];

RW_c = [cos(psi), -sin(psi), 0;...
        sin(psi), cos(psi), 0;...
        0, 0, 1];

Rc_x = [1, 0, 0;...
        0, cos(theta), -sin(theta);...
        0, sin(theta), cos(theta)];
    
Rx_b = [ cos(phi), 0, sin(phi);...
        0, 1, 0;...
        -sin(phi), 0, cos(phi)];
    
Rc_b = Rc_x * Rx_b;
RW_b = RW_c * Rc_b;

w_bW = RW_b*w_b;

u = control(sigma(t), dsigma(t), d2sigma(t), d3sigma(t), r, dr, RW_b, w_b, params);

d2r = -params.g*zW + u(1)*RW_b(:,3);
dw_b = params.J\(cross(-w_b, params.J*w_b) + u(2:4));

dx = double([dr;w_b;d2r;dw_b]);

end