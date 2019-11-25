function [xcoeff,ycoeff,zcoeff,psicoeff] = genWaypointPolys(waypoints, params)

% t0 = waypoints(1,5);
% tm = waypoints(end,5);
% 
% alpha = tm - t0;
% 
% % Temporal Scaling
% waypoints(:,5) = (waypoints(:,5) - t0)/alpha;

xt = waypoints(:,1);
xt = [xt,waypoints(:,5)];
yt = waypoints(:,2);
yt = [yt,waypoints(:,5)];
zt = waypoints(:,3);
zt = [zt,waypoints(:,5)];
psit = waypoints(:,4);
psit = [psit,waypoints(:,5)];

xcoeff = genMinSnapPoly(xt, params.nx, params.kx);
ycoeff = genMinSnapPoly(yt, params.ny, params.ky);
zcoeff = genMinSnapPoly(zt, params.nz, params.kz);
psicoeff = genMinSnapPoly(psit, params.npsi, params.kpsi);

end