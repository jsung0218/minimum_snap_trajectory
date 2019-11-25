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
             0,1,1,5*pi/4,3];
%              0,1,2,2*pi,4];
         
params.m = size(waypoints,1) - 1;

[xcoeff,ycoeff,zcoeff,psicoeff] = genWaypointPolys(waypoints, params)

plotTraj(xcoeff,ycoeff,zcoeff,psicoeff,params,waypoints(:,5))

