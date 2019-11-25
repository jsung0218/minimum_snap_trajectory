function plotTraj(xcoeff,ycoeff,zcoeff,psicoeff,params,time)
syms t

xpoly = [];
ypoly = [];
zpoly = [];
psipoly = [];

for i = 1:params.m
    xpol = 0;
    ypol = 0;
    zpol = 0;
    psipol = 0;
    for j = 1:params.nx+1
        xpol = xpol + xcoeff((i-1)*(params.nx+1) + j)*(t^(params.nx+1-j));
    end
    xpoly = [xpoly;xpol];
    
    for j = 1:params.ny+1
        ypol = ypol + ycoeff((i-1)*(params.ny+1) + j)*(t^(params.ny+1-j));
    end
    ypoly = [ypoly;ypol];
    
    for j = 1:params.nz+1
        zpol = zpol + zcoeff((i-1)*(params.nz+1) + j)*(t^(params.nz+1-j));
    end
    zpoly = [zpoly;zpol];
    
    for j = 1:params.npsi+1
        psipol = psipol + psicoeff((i-1)*(params.npsi+1) + j)*(t^(params.npsi+1-j));
    end
    psipoly = [psipoly;psipol];
end

% xpoly = vpa(xpoly)
% ypoly = vpa(ypoly)
% zpoly = vpa(zpoly)
% 
% subs(xpoly(1),t,1)
% subs(xpoly(2),t,1)

figure(1)
for i = 1:params.m
    fplot3(xpoly(i),ypoly(i),zpoly(i),[time(i),time(i+1)]);
    hold on;
    plot3(subs(xpoly(i),t,time(i)),subs(ypoly(i),t,time(i)),subs(zpoly(i),t,time(i)),'x');
    hold on;
    plot3(subs(xpoly(i),t,time(i+1)),subs(ypoly(i),t,time(i+1)),subs(zpoly(i),t,time(i+1)),'x');
    hold on;
    grid on;
    axis([-2 2 -2 2 0 2]);
    axis on;
end
drawnow;


end