function u = control(sigma, dsigma, d2sigma, d3sigma, r, dr, R, w_b, params)

r_des = reshape(sigma(1:3),[3,1]);
dr_des = reshape(dsigma(1:3),[3,1]);
d2r_des = reshape(d2sigma(1:3),[3,1]);

temp = reshape(d2sigma(1:3),[3,1]) + params.M*params.g*[0;0;1];
zb_des = temp/norm(temp);

xc_des = [cos(sigma(4));sin(sigma(4));0];

yb_des = cross(zb_des, xc_des)/norm(cross(zb_des, xc_des));

xb_des = cross(yb_des, zb_des);

h_w_des = (1/norm(temp))*(reshape(d3sigma(1:3),[3,1]) - dot(zb_des, reshape(d3sigma(1:3),[3,1]))*zb_des);

w_b_des = [dot(-h_w_des, yb_des);...
        dot(h_w_des, xb_des);...
        dot(dsigma(4)*[0;0;1], zb_des)];
    
RW_b_des = [xb_des,yb_des,zb_des];

ep = r - r_des;
ev = dr - dr_des;

F_des = params.Kp*ep + params.Kv*ev + params.M*([0;0;params.g] + d2r_des);

u1 = dot(F_des,R(:,3));

er = 0.5*vee(RW_b_des'*R - R'*RW_b_des);
ew = w_b - w_b_des;

u234 = - params.Kr*er - params.Kw*ew;

u = double([u1;u234]);

end