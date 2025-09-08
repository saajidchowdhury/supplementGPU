setenv('TZ', 'America/New_York');
fclose('all');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',50); %get(groot,'factory')
set(groot,'defaultAxesLineWidth',4);
set(groot,'defaultLineLineWidth',4);
set(groot,'defaultLineMarkerSize',50);
set(groot,'defaultErrorbarLineWidth',4);
set(groot,'defaultErrorbarMarkerSize',50);
set(groot,'defaultErrorbarCapSize',20);
set(groot,'defaultAxesView',[0,90]);
set(groot,'defaultAxesBox','on');
set(groot,'defaultTextFontSize',50);
set(groot,'defaultConstantlineLineWidth',4);
set(groot,'defaultConstantlineAlpha',1);
set(groot,'defaultAxesLabelFontSizeMultiplier',1);
set(groot,'defaultFigurePosition',[1317 541 740 665]);
main = string(datetime('now','Format','M-d-y@HH.mm.ss'))+"fig1";
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
b0 = 0; %au
r0 = 5000; %au
theta0 = 1.01921414770684326534;
phi0 = 0.22205754524899390390;
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

ax = -2.982e-4;
ay = ax;
az = -2*ax;
qx = 0.219;
qy = -qx;
qz = 0;
fRF = 2.5e6 * tau; %au
OmegaRF = 2*pi*fRF; %au

fprintf("\nRunning ode45 with initial\ntheta0 = %.20f radians,\nphi0 = %.20f radians.\n",theta0,phi0);
tstart = tic;
x0 = r0*sin(theta0)*cos(phi0) + b0*sin(theta0+pi/2)*cos(phi0);
y00 = r0*sin(theta0)*sin(phi0) + b0*sin(theta0+pi/2)*sin(phi0);
z0 = r0*cos(theta0) + b0*cos(theta0+pi/2);
vx0 = -v0*sin(theta0)*cos(phi0);
vy0 = -v0*sin(theta0)*sin(phi0);
vz0 = -v0*cos(theta0);
tspan = [0 tmax];
y0 = [0;0;0;0;0;0; x0;y00;z0;vx0;vy0;vz0];
options = odeset('Stats','on','RelTol',rtol,'AbsTol',atol,'OutputFcn',@halt);
[t,y] = ode45(@(t,y)f(C4,C8,mion,matom,ax,ay,az,qx,qy,qz,OmegaRF,t,y), tspan, y0, options);
r = sqrt((y(:,1)-y(:,7)).^2+(y(:,2)-y(:,8)).^2+(y(:,3)-y(:,9)).^2); %meters
KE = 1/2 * mion * (y(end,4)^2+y(end,5)^2+y(end,6)^2);
fprintf("With ode45, scattering angle is %.20e radians.\n", ...
    acos((vx0*y(end,10)+vy0*y(end,11)+vz0*y(end,12))/sqrt(vx0^2+vy0^2+vz0^2)/sqrt(y(end,10)^2+y(end,11)^2+y(end,12)^2)));
tend = toc(tstart);
fprintf("It took %.3f seconds.\n",tend);

fprintf("\nRunning gpu-optimized function with initial\ntheta0 = %.20f radians,\nphi0 = %.20f radians.\n",theta0,phi0);
tstart = tic;
[~,yend1,yend2,yend3,yend4,yend5,yend6,yend7,yend8,yend9,yend10,yend11,yend12] = ...
    g(C4,C8,mion,matom,ax,ay,az,qx,qy,qz,OmegaRF, ...
    0,tmax,0,0,0,0,0,0,x0,y00,z0,vx0,vy0,vz0,rtol,atol);
fprintf("With gpu-optimized function, scattering angle is %.20e radians.\n", ...
    acos((vx0*yend10+vy0*yend11+vz0*yend12)/sqrt(vx0^2+vy0^2+vz0^2)/sqrt(yend10^2+yend11^2+yend12^2)));
tend = toc(tstart);
fprintf("It took %.3f seconds.\n",tend);

l = length(y(:,1));
fraction = 1;
figure; hold on;
plot(y(1:floor(fraction*l),2),y(1:floor(fraction*l),1)); %transposed!!
plot(y(1:floor(fraction*l),8),y(1:floor(fraction*l),7)); %transposed!!
xlim([-1285,2100]);
ylim([-1035,1060]);
xlabel("$x$ (Bohr)");
ylabel("$y$ (Bohr)");
legend("ion","atom");
print(gcf,'-vector','-dsvg',main+"/fig1.svg");
hold off;

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
    %y(1:3) = position of ion
    %y(4:6) = velocity of ion
    %y(7:9) = position of atom
    %y(10:12) = velocity of atom
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

%wherever it says ode, compute the f(t,y) function
function [tend,yend1,yend2,yend3,yend4,yend5,yend6,yend7,yend8,yend9,yend10,yend11,yend12] = ...
    g(C4,C8,mion,matom,ax,ay,az,qx,qy,qz,OmegaRF,t0,tfinal, ...
    y01,y02,y03,y04,y05,y06,y07,y08,y09,y010,y011,y012,rtol,atol)

    nsteps  = 0;
    f01 = y04; f02 = y05; f03 = y06; %ode start
    f07 = y010; f08 = y011; f09 = y012;
    r = sqrt((y07-y01)^2+(y08-y02)^2+(y09-y03)^2);
    dVdr = 4*C4/r^5 - 8*C8/r^9;
    f04 = -1/mion*dVdr*(y01-y07)/r;
    f05 = -1/mion*dVdr*(y02-y08)/r;
    f06 = -1/mion*dVdr*(y03-y09)/r;
    f010 = -1/matom*dVdr*(y07-y01)/r;
    f011 = -1/matom*dVdr*(y08-y02)/r;
    f012 = -1/matom*dVdr*(y09-y03)/r;
    f04 = f04 - (ax+2*qx*cos(OmegaRF*t0))*OmegaRF^2/4*y01;
    f05 = f05 - (ay+2*qy*cos(OmegaRF*t0))*OmegaRF^2/4*y02;
    f06 = f06 - (az+2*qz*cos(OmegaRF*t0))*OmegaRF^2/4*y03; %odeends
    threshold = atol / rtol;
    hmax = min(tfinal-t0,max(0.1*(tfinal-t0),16.0*eps*(tfinal-t0)));
    t = t0;
    y1 = y01; y2 = y02; y3 = y03; y4 = y04; y5 = y05; y6 = y06;
    y7 = y07; y8 = y08; y9 = y09; y10 = y010; y11 = y011; y12 = y012;
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
    absh = min(hmax, tfinal-t0);
    maxabs = abs(f01 / threshold);
    maxabs = max(maxabs, abs(f02 / threshold));
    maxabs = max(maxabs, abs(f03 / threshold));
    maxabs = max(maxabs, abs(f04 / threshold));
    maxabs = max(maxabs, abs(f05 / threshold));
    maxabs = max(maxabs, abs(f06 / threshold));
    maxabs = max(maxabs, abs(f07 / max(abs(y07),threshold)));
    maxabs = max(maxabs, abs(f08 / max(abs(y08),threshold)));
    maxabs = max(maxabs, abs(f09 / max(abs(y09),threshold)));
    maxabs = max(maxabs, abs(f010 / max(abs(y010),threshold)));
    maxabs = max(maxabs, abs(f011 / max(abs(y011),threshold)));
    maxabs = max(maxabs, abs(f012 / max(abs(y012),threshold)));
    rh = maxabs / (0.8 * rtol^pow);
    if absh * rh > 1
        absh = 1 / rh;
    end
    absh = max(absh, hmin);
    f11 = f01; f12 = f02; f13 = f03; f14 = f04; f15 = f05; f16 = f06;
    f17 = f07; f18 = f08; f19 = f09; f110 = f010; f111 = f011; f112 = f012;
    done = false;
    err = eps; %uninit
    tnew = t; %uninit
    ynew1 = y1; ynew2 = y2; ynew3 = y3; ynew4 = y4; ynew5 = y5; ynew6 = y6; %uninit
    ynew7 = y7; ynew8 = y8; ynew9 = y9; ynew10 = y10; ynew11 = y11; ynew12 = y12; %uninit
    f71 = f11; f72 = f12; f73 = f13; f74 = f14; f75 = f15; f76 = f16; %uninit
    f77 = f17; f78 = f18; f79 = f19; f710 = f110; f711 = f111; f712 = f112; %uninit
    while ~done
        hmin = 16*eps(t);
        absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
        h = absh;
        if 1.1*absh >= abs(tfinal - t)
            h = tfinal - t;
            absh = abs(h);
            done = true;
        end
        nofailed = true;                      % no failed attempts
        while true
            y21 = y1 + h * (b11*f11 );
            y22 = y2 + h * (b11*f12 );
            y23 = y3 + h * (b11*f13 );
            y24 = y4 + h * (b11*f14 );
            y25 = y5 + h * (b11*f15 );
            y26 = y6 + h * (b11*f16 );
            y27 = y7 + h * (b11*f17 );
            y28 = y8 + h * (b11*f18 );
            y29 = y9 + h * (b11*f19 );
            y210 = y10 + h * (b11*f110 );
            y211 = y11 + h * (b11*f111 );
            y212 = y12 + h * (b11*f112 );
            t2 = t + h * a2;
            f21 = y24; f22 = y25; f23 = y26; %ode start
            f27 = y210; f28 = y211; f29 = y212;
            r = sqrt((y27-y21)^2+(y28-y22)^2+(y29-y23)^2);
            dVdr = 4*C4/r^5 - 8*C8/r^9;
            f24 = -1/mion*dVdr*(y21-y27)/r;
            f25 = -1/mion*dVdr*(y22-y28)/r;
            f26 = -1/mion*dVdr*(y23-y29)/r;
            f210 = -1/matom*dVdr*(y27-y21)/r;
            f211 = -1/matom*dVdr*(y28-y22)/r;
            f212 = -1/matom*dVdr*(y29-y23)/r;
            f24 = f24 - (ax+2*qx*cos(OmegaRF*t2))*OmegaRF^2/4*y21;
            f25 = f25 - (ay+2*qy*cos(OmegaRF*t2))*OmegaRF^2/4*y22;
            f26 = f26 - (az+2*qz*cos(OmegaRF*t2))*OmegaRF^2/4*y23; %odeend
            y21 = y1 + h * (b21*f11 + b22*f21 );
            y22 = y2 + h * (b21*f12 + b22*f22 );
            y23 = y3 + h * (b21*f13 + b22*f23 );
            y24 = y4 + h * (b21*f14 + b22*f24 );
            y25 = y5 + h * (b21*f15 + b22*f25 );
            y26 = y6 + h * (b21*f16 + b22*f26 );
            y27 = y7 + h * (b21*f17 + b22*f27 );
            y28 = y8 + h * (b21*f18 + b22*f28 );
            y29 = y9 + h * (b21*f19 + b22*f29 );
            y210 = y10 + h * (b21*f110 + b22*f210 );
            y211 = y11 + h * (b21*f111 + b22*f211 );
            y212 = y12 + h * (b21*f112 + b22*f212 );
            t2 = t + h * a3;
            f31 = y24; f32 = y25; f33 = y26; %ode start
            f37 = y210; f38 = y211; f39 = y212;
            r = sqrt((y27-y21)^2+(y28-y22)^2+(y29-y23)^2);
            dVdr = 4*C4/r^5 - 8*C8/r^9;
            f34 = -1/mion*dVdr*(y21-y27)/r;
            f35 = -1/mion*dVdr*(y22-y28)/r;
            f36 = -1/mion*dVdr*(y23-y29)/r;
            f310 = -1/matom*dVdr*(y27-y21)/r;
            f311 = -1/matom*dVdr*(y28-y22)/r;
            f312 = -1/matom*dVdr*(y29-y23)/r;
            f34 = f34 - (ax+2*qx*cos(OmegaRF*t2))*OmegaRF^2/4*y21;
            f35 = f35 - (ay+2*qy*cos(OmegaRF*t2))*OmegaRF^2/4*y22;
            f36 = f36 - (az+2*qz*cos(OmegaRF*t2))*OmegaRF^2/4*y23; %odeend
            y21 = y1 + h * (b31*f11 + b32*f21 + b33*f31 );
            y22 = y2 + h * (b31*f12 + b32*f22 + b33*f32 );
            y23 = y3 + h * (b31*f13 + b32*f23 + b33*f33 );
            y24 = y4 + h * (b31*f14 + b32*f24 + b33*f34 );
            y25 = y5 + h * (b31*f15 + b32*f25 + b33*f35 );
            y26 = y6 + h * (b31*f16 + b32*f26 + b33*f36 );
            y27 = y7 + h * (b31*f17 + b32*f27 + b33*f37 );
            y28 = y8 + h * (b31*f18 + b32*f28 + b33*f38 );
            y29 = y9 + h * (b31*f19 + b32*f29 + b33*f39 );
            y210 = y10 + h * (b31*f110 + b32*f210 + b33*f310 );
            y211 = y11 + h * (b31*f111 + b32*f211 + b33*f311 );
            y212 = y12 + h * (b31*f112 + b32*f212 + b33*f312 );
            t2 = t + h * a4;
            f41 = y24; f42 = y25; f43 = y26; %ode start
            f47 = y210; f48 = y211; f49 = y212;
            r = sqrt((y27-y21)^2+(y28-y22)^2+(y29-y23)^2);
            dVdr = 4*C4/r^5 - 8*C8/r^9;
            f44 = -1/mion*dVdr*(y21-y27)/r;
            f45 = -1/mion*dVdr*(y22-y28)/r;
            f46 = -1/mion*dVdr*(y23-y29)/r;
            f410 = -1/matom*dVdr*(y27-y21)/r;
            f411 = -1/matom*dVdr*(y28-y22)/r;
            f412 = -1/matom*dVdr*(y29-y23)/r;
            f44 = f44 - (ax+2*qx*cos(OmegaRF*t2))*OmegaRF^2/4*y21;
            f45 = f45 - (ay+2*qy*cos(OmegaRF*t2))*OmegaRF^2/4*y22;
            f46 = f46 - (az+2*qz*cos(OmegaRF*t2))*OmegaRF^2/4*y23; %odeend
            y21 = y1 + h * (b41*f11 + b42*f21 + b43*f31 + b44*f41 );
            y22 = y2 + h * (b41*f12 + b42*f22 + b43*f32 + b44*f42 );
            y23 = y3 + h * (b41*f13 + b42*f23 + b43*f33 + b44*f43 );
            y24 = y4 + h * (b41*f14 + b42*f24 + b43*f34 + b44*f44 );
            y25 = y5 + h * (b41*f15 + b42*f25 + b43*f35 + b44*f45 );
            y26 = y6 + h * (b41*f16 + b42*f26 + b43*f36 + b44*f46 );
            y27 = y7 + h * (b41*f17 + b42*f27 + b43*f37 + b44*f47 );
            y28 = y8 + h * (b41*f18 + b42*f28 + b43*f38 + b44*f48 );
            y29 = y9 + h * (b41*f19 + b42*f29 + b43*f39 + b44*f49 );
            y210 = y10 + h * (b41*f110 + b42*f210 + b43*f310 + b44*f410 );
            y211 = y11 + h * (b41*f111 + b42*f211 + b43*f311 + b44*f411 );
            y212 = y12 + h * (b41*f112 + b42*f212 + b43*f312 + b44*f412 );
            t2 = t + h * a5;
            f51 = y24; f52 = y25; f53 = y26; %ode start
            f57 = y210; f58 = y211; f59 = y212;
            r = sqrt((y27-y21)^2+(y28-y22)^2+(y29-y23)^2);
            dVdr = 4*C4/r^5 - 8*C8/r^9;
            f54 = -1/mion*dVdr*(y21-y27)/r;
            f55 = -1/mion*dVdr*(y22-y28)/r;
            f56 = -1/mion*dVdr*(y23-y29)/r;
            f510 = -1/matom*dVdr*(y27-y21)/r;
            f511 = -1/matom*dVdr*(y28-y22)/r;
            f512 = -1/matom*dVdr*(y29-y23)/r;
            f54 = f54 - (ax+2*qx*cos(OmegaRF*t2))*OmegaRF^2/4*y21;
            f55 = f55 - (ay+2*qy*cos(OmegaRF*t2))*OmegaRF^2/4*y22;
            f56 = f56 - (az+2*qz*cos(OmegaRF*t2))*OmegaRF^2/4*y23; %odeend
            y21 = y1 + h * (b51*f11 + b52*f21 + b53*f31 + b54*f41 + b55*f51 );
            y22 = y2 + h * (b51*f12 + b52*f22 + b53*f32 + b54*f42 + b55*f52 );
            y23 = y3 + h * (b51*f13 + b52*f23 + b53*f33 + b54*f43 + b55*f53 );
            y24 = y4 + h * (b51*f14 + b52*f24 + b53*f34 + b54*f44 + b55*f54 );
            y25 = y5 + h * (b51*f15 + b52*f25 + b53*f35 + b54*f45 + b55*f55 );
            y26 = y6 + h * (b51*f16 + b52*f26 + b53*f36 + b54*f46 + b55*f56 );
            y27 = y7 + h * (b51*f17 + b52*f27 + b53*f37 + b54*f47 + b55*f57 );
            y28 = y8 + h * (b51*f18 + b52*f28 + b53*f38 + b54*f48 + b55*f58 );
            y29 = y9 + h * (b51*f19 + b52*f29 + b53*f39 + b54*f49 + b55*f59 );
            y210 = y10 + h * (b51*f110 + b52*f210 + b53*f310 + b54*f410 + b55*f510 );
            y211 = y11 + h * (b51*f111 + b52*f211 + b53*f311 + b54*f411 + b55*f511 );
            y212 = y12 + h * (b51*f112 + b52*f212 + b53*f312 + b54*f412 + b55*f512 );
            t2 = t + h;
            f61 = y24; f62 = y25; f63 = y26; %ode start
            f67 = y210; f68 = y211; f69 = y212;
            r = sqrt((y27-y21)^2+(y28-y22)^2+(y29-y23)^2);
            dVdr = 4*C4/r^5 - 8*C8/r^9;
            f64 = -1/mion*dVdr*(y21-y27)/r;
            f65 = -1/mion*dVdr*(y22-y28)/r;
            f66 = -1/mion*dVdr*(y23-y29)/r;
            f610 = -1/matom*dVdr*(y27-y21)/r;
            f611 = -1/matom*dVdr*(y28-y22)/r;
            f612 = -1/matom*dVdr*(y29-y23)/r;
            f64 = f64 - (ax+2*qx*cos(OmegaRF*t2))*OmegaRF^2/4*y21;
            f65 = f65 - (ay+2*qy*cos(OmegaRF*t2))*OmegaRF^2/4*y22;
            f66 = f66 - (az+2*qz*cos(OmegaRF*t2))*OmegaRF^2/4*y23; %odeend
            tnew = t + h;
            if done
                tnew = tfinal;   % Hit end point exactly.
            end
            h = tnew - t;      % Purify h.
            ynew1 = y1 + h* ( b61*f11 + b63*f31 + b64*f41 + b65*f51 + b66*f61 );
            ynew2 = y2 + h* ( b61*f12 + b63*f32 + b64*f42 + b65*f52 + b66*f62 );
            ynew3 = y3 + h* ( b61*f13 + b63*f33 + b64*f43 + b65*f53 + b66*f63 );
            ynew4 = y4 + h* ( b61*f14 + b63*f34 + b64*f44 + b65*f54 + b66*f64 );
            ynew5 = y5 + h* ( b61*f15 + b63*f35 + b64*f45 + b65*f55 + b66*f65 );
            ynew6 = y6 + h* ( b61*f16 + b63*f36 + b64*f46 + b65*f56 + b66*f66 );
            ynew7 = y7 + h* ( b61*f17 + b63*f37 + b64*f47 + b65*f57 + b66*f67 );
            ynew8 = y8 + h* ( b61*f18 + b63*f38 + b64*f48 + b65*f58 + b66*f68 );
            ynew9 = y9 + h* ( b61*f19 + b63*f39 + b64*f49 + b65*f59 + b66*f69 );
            ynew10 = y10 + h* ( b61*f110 + b63*f310 + b64*f410 + b65*f510 + b66*f610 );
            ynew11 = y11 + h* ( b61*f111 + b63*f311 + b64*f411 + b65*f511 + b66*f611 );
            ynew12 = y12 + h* ( b61*f112 + b63*f312 + b64*f412 + b65*f512 + b66*f612 );
            f71 = ynew4; f72 = ynew5; f73 = ynew6; %ode start
            f77 = ynew10; f78 = ynew11; f79 = ynew12;
            r = sqrt((ynew7-ynew1)^2+(ynew8-ynew2)^2+(ynew9-ynew3)^2);
            dVdr = 4*C4/r^5 - 8*C8/r^9;
            f74 = -1/mion*dVdr*(ynew1-ynew7)/r;
            f75 = -1/mion*dVdr*(ynew2-ynew8)/r;
            f76 = -1/mion*dVdr*(ynew3-ynew9)/r;
            f710 = -1/matom*dVdr*(ynew7-ynew1)/r;
            f711 = -1/matom*dVdr*(ynew8-ynew2)/r;
            f712 = -1/matom*dVdr*(ynew9-ynew3)/r;
            f74 = f74 - (ax+2*qx*cos(OmegaRF*tnew))*OmegaRF^2/4*ynew1;
            f75 = f75 - (ay+2*qy*cos(OmegaRF*tnew))*OmegaRF^2/4*ynew2;
            f76 = f76 - (az+2*qz*cos(OmegaRF*tnew))*OmegaRF^2/4*ynew3; %odeend
            fE1 = f11*e1 + f31*e3 + f41*e4 + f51*e5 + f61*e6 + f71*e7;
            fE2 = f12*e1 + f32*e3 + f42*e4 + f52*e5 + f62*e6 + f72*e7;
            fE3 = f13*e1 + f33*e3 + f43*e4 + f53*e5 + f63*e6 + f73*e7;
            fE4 = f14*e1 + f34*e3 + f44*e4 + f54*e5 + f64*e6 + f74*e7;
            fE5 = f15*e1 + f35*e3 + f45*e4 + f55*e5 + f65*e6 + f75*e7;
            fE6 = f16*e1 + f36*e3 + f46*e4 + f56*e5 + f66*e6 + f76*e7;
            fE7 = f17*e1 + f37*e3 + f47*e4 + f57*e5 + f67*e6 + f77*e7;
            fE8 = f18*e1 + f38*e3 + f48*e4 + f58*e5 + f68*e6 + f78*e7;
            fE9 = f19*e1 + f39*e3 + f49*e4 + f59*e5 + f69*e6 + f79*e7;
            fE10 = f110*e1 + f310*e3 + f410*e4 + f510*e5 + f610*e6 + f710*e7;
            fE11 = f111*e1 + f311*e3 + f411*e4 + f511*e5 + f611*e6 + f711*e7;
            fE12 = f112*e1 + f312*e3 + f412*e4 + f512*e5 + f612*e6 + f712*e7;
            maxabs =             abs(fE1 / max(max(abs(y1),abs(ynew1)),threshold));
            maxabs = max(maxabs, abs(fE2 / max(max(abs(y2),abs(ynew2)),threshold)));
            maxabs = max(maxabs, abs(fE3 / max(max(abs(y3),abs(ynew3)),threshold)));
            maxabs = max(maxabs, abs(fE4 / max(max(abs(y4),abs(ynew4)),threshold)));
            maxabs = max(maxabs, abs(fE5 / max(max(abs(y5),abs(ynew5)),threshold)));
            maxabs = max(maxabs, abs(fE6 / max(max(abs(y6),abs(ynew6)),threshold)));
            maxabs = max(maxabs, abs(fE7 / max(max(abs(y7),abs(ynew7)),threshold)));
            maxabs = max(maxabs, abs(fE8 / max(max(abs(y8),abs(ynew8)),threshold)));
            maxabs = max(maxabs, abs(fE9 / max(max(abs(y9),abs(ynew9)),threshold)));
            maxabs = max(maxabs, abs(fE10 / max(max(abs(y10),abs(ynew10)),threshold)));
            maxabs = max(maxabs, abs(fE11 / max(max(abs(y11),abs(ynew11)),threshold)));
            maxabs = max(maxabs, abs(fE12 / max(max(abs(y12),abs(ynew12)),threshold)));
            err = absh * maxabs;
            if err > rtol                       % Failed step
                if nofailed
                    nofailed = false;
                    absh = max(hmin, absh * max(0.1, 0.8*(rtol/err)^pow));
                else
                    absh = max(hmin, 0.5 * absh);
                end
                h = absh;
                done = false;
            else                                % Successful step
                break;
            end
        end
        nsteps = nsteps + 1;
        if nsteps >= 1000000
            done = true;
        end
        if done
            break;
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
        y1 = ynew1; y2 = ynew2; y3 = ynew3; y4 = ynew4; y5 = ynew5; y6 = ynew6;
        y7 = ynew7; y8 = ynew8; y9 = ynew9; y10 = ynew10; y11 = ynew11; y12 = ynew12;
        f11 = f71; f12 = f72; f13 = f73; f14 = f74; f15 = f75; f16 = f76; % Already have f(tnew,ynew)
        f17 = f77; f18 = f78; f19 = f79; f110 = f710; f111 = f711; f112 = f712;
    end % end of while loop
    tend = tnew;
    yend1 = ynew1; yend2 = ynew2; yend3 = ynew3; yend4 = ynew4; yend5 = ynew5; yend6 = ynew6;
    yend7 = ynew7; yend8 = ynew8; yend9 = ynew9; yend10 = ynew10; yend11 = ynew11; yend12 = ynew12;
end