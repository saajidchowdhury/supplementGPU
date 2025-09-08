clear all;
setenv('TZ', 'America/New_York');
fclose('all');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',50); %get(groot,'factory')
set(groot,'defaultAxesLineWidth',3);
set(groot,'defaultLineLineWidth',3);
set(groot,'defaultLineMarkerSize',50);
set(groot,'defaultErrorbarLineWidth',3);
set(groot,'defaultErrorbarMarkerSize',50);
set(groot,'defaultErrorbarCapSize',20);
set(groot,'defaultAxesView',[0,90]);
set(groot,'defaultAxesBox','on');
set(groot,'defaultTextFontSize',50);
set(groot,'defaultConstantlineLineWidth',3);
set(groot,'defaultConstantlineAlpha',1);
set(groot,'defaultAxesLabelFontSizeMultiplier',1);
set(groot,'defaultFigurePosition',[790 1 1267 1173]);
main = string(datetime('now','Format','M-d-y@HH.mm.ss'))+"code4olderversion";
mkdir(main);

a0 = 5.29177210544e-11; %meters
hbar = 1.054571817e-34; %Joules*seconds
Eh = 4.3597447222060e-18; %Joules
tau = 2.4188843265864e-17; %seconds. tau = hbar/Eh
vau = 2.18769126216e6; %meters per second
me = 9.1093837139e-31; %kilograms
u = 1.66053906892e-27; %kilograms
amu = u/me; %1822.9 electron masses, au
kB = 1.380649e-23; %Joules per Kelvin
EhK = 3.1577502480398e5; %Kelvin

tmax = 1000e-6 / tau; %seconds -> au.
r0 = 5000; %au
umin = 0; umax = 1/2; ntheta = 50;
phimin = 0; phimax = pi/2; nphi = 50;
us = linspace(umin,umax,ntheta);
thetas = acos(1-2*us);
phis = linspace(phimin,phimax,nphi);
T = 1e-6; %Kelvin
mion = 170.936323 * amu; %u -> au
matom = 86.9091805310 * amu; %u -> au
De = 1 / EhK; %Kelvin -> Hartree
C4 = 319.091/2; %au
C8 = C4^2/(4*De); %au
v0 = sqrt(3*kB*T/Eh/matom); %au
Re = (2*C8/C4)^(1/4); %au
rbounce = Re; %au
rtol = 1e-10;
atol = 1e-20;
options = odeset('Stats','off','RelTol',rtol,'AbsTol',atol,'OutputFcn',@halt);

ax = -2.982e-4;
ay = ax;
az = -2*ax;
qx = 0.219;
qy = -qx;
qz = 0;
fRF = 2.5e6 * tau; %au
OmegaRF = 2*pi*fRF; %au

%gpu scan
fprintf("\nRunning GPU scan.\n");
angleGPU = zeros(ntheta,nphi);
tendGPU = zeros(ntheta,nphi);
for i = 1:ntheta
    for j = 1:nphi
        tstartGPU = tic;
        fprintf("GPU i=%d/%d, j=%d/%d.\n",i,ntheta,j,nphi);
        theta0 = thetas(i);
        phi0 = phis(j);
        x0 = r0*sin(theta0)*cos(phi0);
        y0 = r0*sin(theta0)*sin(phi0);
        z0 = r0*cos(theta0);
        vx0 = -v0*sin(theta0)*cos(phi0);
        vy0 = -v0*sin(theta0)*sin(phi0);
        vz0 = -v0*cos(theta0);
        [tend,yend1,yend2,yend3,yend4,yend5,yend6,yend7,yend8,yend9,yend10,yend11,yend12] = ...
            ode45gpu(C4,C8,mion,matom,ax,ay,az,qx,qy,qz,OmegaRF, ...
            0,tmax,0,0,0,0,0,0,x0,y0,z0,vx0,vy0,vz0,rtol,atol);
        angleGPU(i,j) = acos((vx0*yend10+vy0*yend11+vz0*yend12)/sqrt(vx0^2+vy0^2+vz0^2)/sqrt(yend10^2+yend11^2+yend12^2));
        tendGPU(i,j) = toc(tstartGPU);
    end
end

%cpu scan
fprintf("\nRunning CPU scan.\n");
angleCPU = zeros(ntheta,nphi);
tendCPU = zeros(ntheta,nphi);
for i = 1:ntheta
    for j = 1:nphi
        tstartCPU = tic;
        fprintf("CPU i=%d/%d, j=%d/%d.\n",i,ntheta,j,nphi);
        theta0 = thetas(i);
        phi0 = phis(j);
        x0 = r0*sin(theta0)*cos(phi0);
        y0 = r0*sin(theta0)*sin(phi0);
        z0 = r0*cos(theta0);
        vx0 = -v0*sin(theta0)*cos(phi0);
        vy0 = -v0*sin(theta0)*sin(phi0);
        vz0 = -v0*cos(theta0);
        y00 = [0;0;0;0;0;0;x0;y0;z0;vx0;vy0;vz0];
        [~,y] = ode45(@(t,y)f(C4,C8,mion,matom,ax,ay,az,qx,qy,qz,OmegaRF,t,y), [0,tmax], y00, options);
        angleCPU(i,j) = acos((vx0*y(end,10)+vy0*y(end,11)+vz0*y(end,12))/sqrt(vx0^2+vy0^2+vz0^2)/sqrt(y(end,10)^2+y(end,11)^2+y(end,12)^2));
        tendCPU(i,j) = toc(tstartCPU);
    end
end

fprintf("\nThe CPU took %.3f seconds.\n",sum(tendCPU(:)));
fprintf("The GPU took %.3f seconds.\n",sum(tendGPU(:)));
fprintf("The maximum difference in CPU and GPU scattering angles is %.2e radians.\n",max(angleCPU(:)-angleGPU(:)));

figure; hold on;
tendCPUplot = tendCPU';
tendGPUplot = tendGPU';
plot(1:ntheta*nphi, cumsum(tendCPUplot(:))/3600);
plot(1:ntheta*nphi, cumsum(tendGPUplot(:))/3600);
legend(["\verb+ode45+","older version"],"Location","Northwest");
xlabel("Number of collisions, $n$");
ylabel("Cumulative runtime, $t$ (hours)");
print(gcf,'-vector','-dsvg',main+"/time.svg");
saveas(gcf,main+"/time.png");
hold off;

figure; hold on;
imagesc(us,phis,angleCPU');
axis([umin,umax,phimin,phimax]);
set(gca,'YDir','normal');
set(gcf,'Position', [790 1 1267 1173]); %get(0, 'Screensize')
xlabel("$\theta$ (radians)");
ylabel("$\phi$ (radians)");
xticks(0.5*(1-cos([0 0.5 1 1.5])));
xticklabels({'0','0.5','1','1.5'});
colormap(jet);
c = colorbar;
c.Label.Interpreter = "Latex";
c.Label.String = "Scattering angle (radians)";
c.TickLabelInterpreter = "Latex";
print(gcf,'-vector','-dsvg',main+"/angle.svg");
saveas(gcf,main+'/angle.png');
hold off;

clear y;
save(main+"/fig4olderversion.mat");

function status = halt(t,~,flag)
    status = 0;
    persistent nstep
    switch flag
        case 'init'
            nstep = 0;
        case []
            nstep = nstep + length(t);
            if nstep >= 4000000
                status = 1;
            end
    end
end

function f = f(C4,C8,mion,matom,ax,ay,az,qx,qy,qz,OmegaRF,t,y)
    f = zeros(12,1);
    f(1:3) = y(4:6);
    f(7:9) = y(10:12);
    r = sqrt((y(1)-y(7))^2+(y(2)-y(8))^2+(y(3)-y(9))^2);
    dVdr = 4*C4/r^5 - 8*C8/r^9;
    f(4:6) = -1/mion*dVdr*(y(1:3)-y(7:9))/r;
    f(10:12) = -1/matom*dVdr*(y(7:9)-y(1:3))/r;
    f(4) = f(4) - (ax+2*qx*cos(OmegaRF*t))*OmegaRF^2/4*y(1);
    f(5) = f(5) - (ay+2*qy*cos(OmegaRF*t))*OmegaRF^2/4*y(2);
    f(6) = f(6) - (az+2*qz*cos(OmegaRF*t))*OmegaRF^2/4*y(3);
end

function [tend,yend1,yend2,yend3,yend4,yend5,yend6,yend7,yend8,yend9,yend10,yend11,yend12] = ...
    ode45gpu(C4,C8,mion,matom,ax,ay,az,qx,qy,qz,OmegaRF, ...
    t0,tf,y01,y02,y03,y04,y05,y06,y07,y08,y09,y010,y011,y012,rtol,atol)
    
    y0 = [y01;y02;y03;y04;y05;y06;y07;y08;y09;y010;y011;y012];
    nsteps = 0;
    htspan = tf - t0;
    %ode start
    f0 = zeros(12,1);
    f0(1:3) = y0(4:6);
    f0(7:9) = y0(10:12);
    r = sqrt((y0(1)-y0(7))^2+(y0(2)-y0(8))^2+(y0(3)-y0(9))^2);
    dVdr = 4*C4/r^5 - 8*C8/r^9;
    f0(4:6) = -1/mion*dVdr*(y0(1:3)-y0(7:9))/r;
    f0(10:12) = -1/matom*dVdr*(y0(7:9)-y0(1:3))/r;
    f0(4) = f0(4) - (ax+2*qx*cos(OmegaRF*t0))*OmegaRF^2/4*y0(1);
    f0(5) = f0(5) - (ay+2*qy*cos(OmegaRF*t0))*OmegaRF^2/4*y0(2);
    f0(6) = f0(6) - (az+2*qz*cos(OmegaRF*t0))*OmegaRF^2/4*y0(3);
    %ode end
    threshold = atol ./ rtol;
    hmax = min(tf-t0,max(0.1*(tf-t0),16.0*eps*tf));
    t = t0;
    y = y0;
    pow = 1/5;
    a2=1/5;
    a3=3/10;
    a4=4/5;
    a5=8/9;
    b11=1/5;
    b21=3/40;
    b31=44/45;
    b41=19372/6561;
    b51=9017/3168;
    b61=35/384;
    b22=9/40;
    b32=-56/15;
    b42=-25360/2187;
    b52=-355/33;
    b33=32/9;
    b43=64448/6561;
    b53=46732/5247;
    b63=500/1113;
    b44=-212/729;
    b54=49/176;
    b64=125/192;
    b55=-5103/18656;
    b65=-2187/6784;
    b66=11/84;
    e1=71/57600;
    e3=-71/16695;
    e4=71/1920;
    e5=-17253/339200;
    e6=22/525;
    e7=-1/40;
    hmin = 16*eps(t);
    absh = min(hmax, htspan);
    rh = norm(f0 ./ max(abs(y),threshold),inf) / (0.8 * rtol^pow);
    if absh * rh > 1
        absh = 1 / rh;
    end
    absh = max(absh, hmin);
    f1 = f0;
    
    
    
    done = false;
    while ~done
        hmin = 16*eps(t);
        absh = min(hmax, max(hmin, absh));
        h = absh;
        if 1.1*absh >= abs(tf - t)
            h = tf - t;
            absh = abs(h);
            done = true;
        end
        nofailed = true;
        while true
            y2 = y + h .* (b11.*f1 );
            t2 = t + h .* a2;
            %ode start
            f2 = zeros(12,1);
            f2(1:3) = y2(4:6);
            f2(7:9) = y2(10:12);
            r = sqrt((y2(1)-y2(7))^2+(y2(2)-y2(8))^2+(y2(3)-y2(9))^2);
            dVdr = 4*C4/r^5 - 8*C8/r^9;
            f2(4:6) = -1/mion*dVdr*(y2(1:3)-y2(7:9))/r;
            f2(10:12) = -1/matom*dVdr*(y2(7:9)-y2(1:3))/r;
            f2(4) = f2(4) - (ax+2*qx*cos(OmegaRF*t2))*OmegaRF^2/4*y2(1);
            f2(5) = f2(5) - (ay+2*qy*cos(OmegaRF*t2))*OmegaRF^2/4*y2(2);
            f2(6) = f2(6) - (az+2*qz*cos(OmegaRF*t2))*OmegaRF^2/4*y2(3);
            %ode end
            y3 = y + h .* (b21.*f1 + b22.*f2 );
            t3 = t + h .* a3;
            %ode start
            f3 = zeros(12,1);
            f3(1:3) = y3(4:6);
            f3(7:9) = y3(10:12);
            r = sqrt((y3(1)-y3(7))^2+(y3(2)-y3(8))^2+(y3(3)-y3(9))^2);
            dVdr = 4*C4/r^5 - 8*C8/r^9;
            f3(4:6) = -1/mion*dVdr*(y3(1:3)-y3(7:9))/r;
            f3(10:12) = -1/matom*dVdr*(y3(7:9)-y3(1:3))/r;
            f3(4) = f3(4) - (ax+2*qx*cos(OmegaRF*t3))*OmegaRF^2/4*y3(1);
            f3(5) = f3(5) - (ay+2*qy*cos(OmegaRF*t3))*OmegaRF^2/4*y3(2);
            f3(6) = f3(6) - (az+2*qz*cos(OmegaRF*t3))*OmegaRF^2/4*y3(3);
            %ode end
            y4 = y + h .* (b31.*f1 + b32.*f2 + b33.*f3 );
            t4 = t + h .* a4;
            %ode start
            f4 = zeros(12,1);
            f4(1:3) = y4(4:6);
            f4(7:9) = y4(10:12);
            r = sqrt((y4(1)-y4(7))^2+(y4(2)-y4(8))^2+(y4(3)-y4(9))^2);
            dVdr = 4*C4/r^5 - 8*C8/r^9;
            f4(4:6) = -1/mion*dVdr*(y4(1:3)-y4(7:9))/r;
            f4(10:12) = -1/matom*dVdr*(y4(7:9)-y4(1:3))/r;
            f4(4) = f4(4) - (ax+2*qx*cos(OmegaRF*t4))*OmegaRF^2/4*y4(1);
            f4(5) = f4(5) - (ay+2*qy*cos(OmegaRF*t4))*OmegaRF^2/4*y4(2);
            f4(6) = f4(6) - (az+2*qz*cos(OmegaRF*t4))*OmegaRF^2/4*y4(3);
            %ode end
            y5 = y + h .* (b41.*f1 + b42.*f2 + b43.*f3 + b44.*f4 );
            t5 = t + h .* a5;
            %ode start
            f5 = zeros(12,1);
            f5(1:3) = y5(4:6);
            f5(7:9) = y5(10:12);
            r = sqrt((y5(1)-y5(7))^2+(y5(2)-y5(8))^2+(y5(3)-y5(9))^2);
            dVdr = 4*C4/r^5 - 8*C8/r^9;
            f5(4:6) = -1/mion*dVdr*(y5(1:3)-y5(7:9))/r;
            f5(10:12) = -1/matom*dVdr*(y5(7:9)-y5(1:3))/r;
            f5(4) = f5(4) - (ax+2*qx*cos(OmegaRF*t5))*OmegaRF^2/4*y5(1);
            f5(5) = f5(5) - (ay+2*qy*cos(OmegaRF*t5))*OmegaRF^2/4*y5(2);
            f5(6) = f5(6) - (az+2*qz*cos(OmegaRF*t5))*OmegaRF^2/4*y5(3);
            %ode end
            y6 = y + h .* (b51.*f1 + b52.*f2 + b53.*f3 + b54.*f4 + b55.*f5 );
            t6 = t + h;
            %ode start
            f6 = zeros(12,1);
            f6(1:3) = y6(4:6);
            f6(7:9) = y6(10:12);
            r = sqrt((y6(1)-y6(7))^2+(y6(2)-y6(8))^2+(y6(3)-y6(9))^2);
            dVdr = 4*C4/r^5 - 8*C8/r^9;
            f6(4:6) = -1/mion*dVdr*(y6(1:3)-y6(7:9))/r;
            f6(10:12) = -1/matom*dVdr*(y6(7:9)-y6(1:3))/r;
            f6(4) = f6(4) - (ax+2*qx*cos(OmegaRF*t6))*OmegaRF^2/4*y6(1);
            f6(5) = f6(5) - (ay+2*qy*cos(OmegaRF*t6))*OmegaRF^2/4*y6(2);
            f6(6) = f6(6) - (az+2*qz*cos(OmegaRF*t6))*OmegaRF^2/4*y6(3);
            %ode end
            tnew = t + h;
            if done
                tnew = tf;
            end
            h = tnew - t;
            ynew = y + h.* ( b61.*f1 + b63.*f3 + b64.*f4 + b65.*f5 + b66.*f6 );
            %ode start
            f7 = zeros(12,1);
            f7(1:3) = ynew(4:6);
            f7(7:9) = ynew(10:12);
            r = sqrt((ynew(1)-ynew(7))^2+(ynew(2)-ynew(8))^2+(ynew(3)-ynew(9))^2);
            dVdr = 4*C4/r^5 - 8*C8/r^9;
            f7(4:6) = -1/mion*dVdr*(ynew(1:3)-ynew(7:9))/r;
            f7(10:12) = -1/matom*dVdr*(ynew(7:9)-ynew(1:3))/r;
            f7(4) = f7(4) - (ax+2*qx*cos(OmegaRF*tnew))*OmegaRF^2/4*ynew(1);
            f7(5) = f7(5) - (ay+2*qy*cos(OmegaRF*tnew))*OmegaRF^2/4*ynew(2);
            f7(6) = f7(6) - (az+2*qz*cos(OmegaRF*tnew))*OmegaRF^2/4*ynew(3);
            %ode end
            fE = f1*e1 + f3*e3 + f4*e4 + f5*e5 + f6*e6 + f7*e7;
            err = absh * norm((fE) ./ max(max(abs(y),abs(ynew)),threshold),inf);
            if err > rtol
                if nofailed
                    nofailed = false;
                    absh = max(hmin, absh * max(0.1, 0.8*(rtol/err)^pow));
                else
                    absh = max(hmin, 0.5 * absh);
                end
                h = absh;
                done = false;
            else
                break;
            end
        end
        nsteps = nsteps + 1;
        if nsteps >= 1000000
            done = true; %this if may be wrong / off by one
        end
        if done
            break
        end
        if nofailed
            temp = 1.25*(err/rtol)^pow;
            if temp > 0.2
                absh = absh / temp;
            else
                absh = 5.0*absh;
            end
        end
        t = tnew;
        y = ynew;
        f1 = f7;  % Already have f(tnew,ynew)
    end % end of while loop
    tend = tnew;
    yend1 = ynew(1); yend2 = ynew(2); yend3 = ynew(3);
    yend4 = ynew(4); yend5 = ynew(5); yend6 = ynew(6);
    yend7 = ynew(7); yend8 = ynew(8); yend9 = ynew(9);
    yend10 = ynew(10); yend11 = ynew(11); yend12 = ynew(12);
end
