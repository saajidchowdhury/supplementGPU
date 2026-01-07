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
main = string(datetime('now','Format','M-d-y@HH.mm.ss'))+"code3";
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
fprintf("\nRunning CPU-improved scan.\n");
angleGPU = zeros(ntheta,nphi);
tendGPU = zeros(ntheta,nphi);
for i = 1:ntheta
    for j = 1:nphi
        tstartGPU = tic;
        fprintf("CPU-improved i=%d/%d, j=%d/%d.\n",i,ntheta,j,nphi);
        theta0 = thetas(i);
        phi0 = phis(j);
        x0 = r0*sin(theta0)*cos(phi0);
        y0 = r0*sin(theta0)*sin(phi0);
        z0 = r0*cos(theta0);
        vx0 = -v0*sin(theta0)*cos(phi0);
        vy0 = -v0*sin(theta0)*sin(phi0);
        vz0 = -v0*cos(theta0);
        y00 = [0;0;0;0;0;0; x0;y0;z0;vx0;vy0;vz0];
        [~,y] = ode45(@(t,y)g(C4,C8,mion,matom,ax,ay,az,qx,qy,qz,OmegaRF,t,y), [0,tmax], y00, options);
        angleGPU(i,j) = acos((vx0*y(end,10)+vy0*y(end,11)+vz0*y(end,12))/sqrt(vx0^2+vy0^2+vz0^2)/sqrt(y(end,10)^2+y(end,11)^2+y(end,12)^2));
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
        y00 = [0;0;0;0;0;0; x0;y0;z0;vx0;vy0;vz0];
        [~,y] = ode45(@(t,y)f(C4,C8,mion,matom,ax,ay,az,qx,qy,qz,OmegaRF,t,y), [0,tmax], y00, options);
        angleCPU(i,j) = acos((vx0*y(end,10)+vy0*y(end,11)+vz0*y(end,12))/sqrt(vx0^2+vy0^2+vz0^2)/sqrt(y(end,10)^2+y(end,11)^2+y(end,12)^2));
        tendCPU(i,j) = toc(tstartCPU);
    end
end

fprintf("\nThe CPU took %.3f seconds.\n",sum(tendCPU(:)));
fprintf("The CPU-improved took %.3f seconds.\n",sum(tendGPU(:)));
fprintf("The maximum difference in CPU and CPU-improved scattering angles is %.2e radians.\n",max(angleCPU(:)-angleGPU(:)));

figure; hold on;
tendCPUplot = tendCPU';
tendGPUplot = tendGPU';
plot(1:ntheta*nphi, cumsum(tendCPUplot(:))/3600);
plot(1:ntheta*nphi, cumsum(tendGPUplot(:))/3600);
legend(["\verb+ode45+","\verb+ode45+, improved $f(t,y)$"],"Location","Northwest");
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
save(main+"/fig3.mat");

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

function f = g(C4,C8,mion,matom,ax,ay,az,qx,qy,qz,OmegaRF,t,y)
    f = zeros(12,1);
    f(1) = y(4); f(2) = y(5); f(3) = y(6);
    f(7) = y(10); f(8) = y(11); f(9) = y(12);
    r = sqrt((y(1)-y(7))^2+(y(2)-y(8))^2+(y(3)-y(9))^2);
    dVdr = 4*C4/r^5 - 8*C8/r^9;
    f(4) = -1/mion*dVdr*(y(1)-y(7))/r;
    f(5) = -1/mion*dVdr*(y(2)-y(8))/r;
    f(6) = -1/mion*dVdr*(y(3)-y(9))/r;
    f(10) = -1/matom*dVdr*(y(7)-y(1))/r;
    f(11) = -1/matom*dVdr*(y(8)-y(2))/r;
    f(12) = -1/matom*dVdr*(y(9)-y(3))/r;
    f(4) = f(4) - (ax+2*qx*cos(OmegaRF*t))*OmegaRF^2/4*y(1);
    f(5) = f(5) - (ay+2*qy*cos(OmegaRF*t))*OmegaRF^2/4*y(2);
    f(6) = f(6) - (az+2*qz*cos(OmegaRF*t))*OmegaRF^2/4*y(3);
end
