clear all;
digits(4)
params.nx = 9;
params.kx = 4;
params.ny = 9;
params.ky = 4;
params.nz = 9;
params.kz = 4;
params.npsi = 6;
params.kpsi = 2;

waypoints = [0,0,1,0,0;
             1,0,1,pi/4,1;
             1,1,1,3*pi/4,2;
             0,1,1,5*pi/4,3;
             0,1,2,2*pi,4];
         
params.m = size(waypoints,1) - 1;

params.M = 1;
params.J = eye(3);
params.Kp = eye(3);
params.Kv = eye(3);
params.Kr = eye(3);
params.Kw = eye(3);
params.g = 9.8;

[xcoeff,ycoeff,zcoeff,psicoeff] = genWaypointPolys(waypoints, params)

plotTraj(xcoeff,ycoeff,zcoeff,psicoeff,params,waypoints(:,5))

dyn = @(t,y)dynamics(t,y,xcoeff,ycoeff,zcoeff,psicoeff,waypoints(:,5),params);

% y0 = [waypoints(1,1),waypoints(1,2),waypoints(1,3),0,0,waypoints(1,4),0,0,0,0,0,0]
% 
% [t,y] = ode45(dyn,[waypoints(1,5),waypoints(end,5)],y0);
% 
% data = [t,y];
% save('case1_data.mat','data');
% 
% 
