syms phi theta psi real  % roll pitch yaw - (Z-X-Y) transform from world to body
syms x y z real          % position - world frame
syms dx dy dz real       % linear velocity - world frame
syms p q r real          % angular velocity - body frame

syms m g real

syms u1 u2 u3 u4 real    % u1-net thrust, u2,u3,u4-moment(x,y,z)

syms d2x d2y d2z d2psi dpsi real
syms d3x d3y d3z d3psi real
syms d4x d4y d4z d4psi real

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

sigma  = [x; y; z; psi];
dsigma  = [dx; dy; dz; dpsi];
d2sigma  = [d2x; d2y; d2z; d2psi];
d3sigma  = [d3x; d3y; d3z; d3psi];   % a_dot

r = sigma(1:3);
dr = dsigma(1:3);
d2r = d2sigma(1:3);

zb = (d2sigma(1:3) + m*g*[0;0;1])/norm(d2sigma(1:3) + m*g*[0;0;1]);

xc = [cos(sigma(4));sin(sigma(4));0];

yb = cross(zb, xc)/norm(cross(zb, xc));

xb = cross(yb, zb);

h_w = (m/u1)*(d3sigma(1:3) - dot(zb, d3sigma(1:3))*zb);

w_bW = [dot(-h_w, yb);...
        dot(h_w, xb);...
        dot(dsigma(4)*zW, zb)];
    



