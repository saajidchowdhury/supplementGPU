clear all;

% Input start (units are a.u., unless otherwise specified):

mion          = 170.936323 * 1822.88848628276014097;
matom         = 86.9091805310 * 1822.88848628276014097;
Tion          = 10.0e-6;   % Kelvin
Tatom         = 1.0e-6;    % Kelvin
potential     = "CnCm";    % "CnCm" or "DeRe" or "CnDe"
n             = 4;         % unitless
m             = 8;         % unitless
Cn            = 159.5455;
Cm            = 2.00949488627505397797e+09;
De            = 1.0 * 3.16681156345647751469e-06;
Re            = 70.8448174250339519631;
ax            = -2.982e-4;
ay            = ax;
az            = -2*ax;
qx            = 0.219;
qy            = -qx;
qz            = 0;
OmegaRF       = 2*pi * 2.5e6 * 2.4188843265864e-17;

tmax          = 1000e-6 / 2.4188843265864e-17;
r0            = 5000;
ntrajectories = 100;       % unitless
nanglesauto   = "no";      % "yes" or "no"
ntheta        = 10;        % unitless
nphi          = 10;        % unitless
thetamin      = 0;         % radians
thetamax      = pi;        % radians
phimin        = 0;         % radians
phimax        = 2*pi;      % radians
rtol          = 1e-10;     % unitless
atol          = 1e-20;     % unitless
maxsteps      = 1000000;   % unitless
ncores        = 1;         % unitless
positions     = "uniform"; % "uniform" or "random"
collisiontype = "head-on"; % "head-on" or "thermal"
processor     = "CPU";     % "GPU" or "CPU"

longestlived  = "yes";     % "yes" or "no"
custom        = "yes";     % "yes" or "no"
t0custom      = 0;
tfcustom      = tmax;
y0custom      = [0,0,0,0,0,0,4154,938,2620,...
                -6.4e-9,-1.45e-9,-4.1e-9];

saveworkspace = "yes";     % "yes" or "no"
savecsvs      = "yes";     % "yes" or "no"
onlyonecsv    = "angle";   % "no" or "angle" or ...

% Input end










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
main = string(datetime('now','Format','M-d-y@HH.mm.ss'))+"atomiongpu";
mkdir(main);

a0 = 5.29177210544e-11; % meters
hbar = 1.054571817e-34; % Joules*seconds
Eh = 4.3597447222060e-18; % Joules
tau = 2.4188843265864e-17; % seconds. tau = hbar/Eh
vau = 2.18769126216e6; % meters per second
me = 9.1093837139e-31; % kilograms
u = 1.66053906892e-27; % kilograms
amu = u/me; % 1822.9 electron masses, au
kB = 1.380649e-23; % Joules per Kelvin
EhK = 3.1577502480398e5; % Kelvin

% m6Li = 6.0151228874 * amu; C4Li = 82.0563;
% m23Na = 22.989767 * amu; C4Na = 81.3512;
% m39K = 38.9637064864 * amu; C4K = 146.4389;
% m87Rb = 86.9091805310 * amu; C4Rb = 319.091/2;

if nanglesauto == "yes"
    ntheta = floor(sqrt(ntrajectories));
    nphi = ceil(ntrajectories / ntheta);
    ntrajectories = ntheta * nphi;
elseif nanglesauto == "no"
    if ntheta * nphi ~= ntrajectories
        error("Error, ntheta * nphi needs to equal ntrajectories.");
    end
else
    error("Error, nanglesauto must be ""yes"" or ""no"".");
end

if collisiontype == "head-on"
    s1 = "colliding head-on at fixed speed";
    s3 = "rest at origin";
elseif collisiontype == "thermal"
    s1 = "in thermal bath colliding";
    s3 = sprintf("temp %.1e Kelvin near origin", Tion);
end
if qx == 0 && qy == 0 && qz == 0
    s2 = "harmonic";
else
    s2 = "Paul";
end
fprintf("\nSimulating %d trajectories of atom (mass %.1f amu, temp %.1e Kelvin) starting %.1f Bohr away from origin, \n" + ...
    "%s with %s-trapped ion (mass %.1f amu) starting at %s,\n" + ...
    "for simulated time %.1f microseconds using %s with %d cores.\n\n", ...
    ntrajectories, matom/amu, Tatom, r0, s1, s2, mion/amu, s3, tmax*tau*1e6, processor, ncores);

if saveworkspace == "yes"
    memoryworkspace = 280*ntrajectories + 34000;
elseif saveworkspace == "no"
    memoryworkspace = 0;
else
    error("Error, saveworkspace must be ""yes"" or ""no"".");
end
if savecsvs == "yes"
    memorycsvs = 64*ntrajectories * 2.5;
elseif savecsvs == "no"
    if any(onlyonecsv == ["angle","bounces","lifetime","position","transfer","dist","nsteps","KE"])
        memorycsvs = 8*ntrajectories * 2.5;
    elseif onlyonecsv == "no"
        memorycsvs = 0;
    else
        error("Error, onlyonecsv must be ""no"" or ""angle"" or ""bounces"" or ""lifetime"" or ""position"" or ""transfer"" or ""dist"" or ""nsteps"" or ""KE"".");
    end
else
    error("Error, savecsvs must be ""yes"" or ""no"".");
end
memoryfigures = 16*1048576;
memorytotal = memoryworkspace + memorycsvs + memoryfigures;
fprintf("We will use approximately %.3f megabytes (MB) of storage (%.3f MB for workspace, %.3f MB for csvs, %.3f MB for figures).\n\n", ...
    memorytotal/1048576, memoryworkspace/1048576, memorycsvs/1048576, memoryfigures/1048576);

if potential == "CnCm"
    Re = (m/n*Cm/Cn)^(1/(m-n));
    De = Cn/Re^n - Cm/Re^m;
elseif potential == "DeRe"
    Cn = De*Re^n * (n/m)^(n/(m-n)) / ((n/m)^(n/(m-n))-(n/m)^(m/(m-n)));
    Cm = De*Re^m * (n/m)^(m/(m-n)) / ((n/m)^(n/(m-n))-(n/m)^(m/(m-n)));
elseif potential == "CnDe"
    Re = (Cn/De)^(1/n) * (((n/m)^(n/(m-n))-(n/m)^(m/(m-n))) / (n/m)^(n/(m-n)))^(1/n);
    Cm = (Cn/Re^n-De) * Re^m;
else
    error("Error, potential must be ""CnCm"" or ""DeRe"".");
end
fprintf("The atom-ion potential is V(r) = -%.3f/r^%d + %.1e/r^%d, with Re = %.1f Bohr and De = %.1e Kelvin.\n\n", ...
    Cn,n,Cm,m,Re,De/3.16681156345647751469e-06);
rbounce = Re;

probability = 0;

tends       = zeros(ntheta,nphi);
yend1s      = zeros(ntheta,nphi);
yend2s      = zeros(ntheta,nphi);
yend3s      = zeros(ntheta,nphi);
yend4s      = zeros(ntheta,nphi);
yend5s      = zeros(ntheta,nphi);
yend6s      = zeros(ntheta,nphi);
yend7s      = zeros(ntheta,nphi);
yend8s      = zeros(ntheta,nphi);
yend9s      = zeros(ntheta,nphi);
yend10s     = zeros(ntheta,nphi);
yend11s     = zeros(ntheta,nphi);
yend12s     = zeros(ntheta,nphi);

angle       = zeros(ntheta,nphi);
bounces     = zeros(ntheta,nphi);
lifetime    = zeros(ntheta,nphi);
position    = zeros(ntheta,nphi);
transfer    = zeros(ntheta,nphi);
dist        = zeros(ntheta,nphi);
nsteps      = zeros(ntheta,nphi);
KE          = zeros(ntheta,nphi);

mions = ones(ntheta,nphi) * mion;
matoms = ones(ntheta,nphi) * matom;
ns = ones(ntheta,nphi) * n;
ms = ones(ntheta,nphi) * m;
Cns = ones(ntheta,nphi) * Cn;
Cms = ones(ntheta,nphi) * Cm;
axs = ones(ntheta,nphi) * ax;
ays = ones(ntheta,nphi) * ay;
azs = ones(ntheta,nphi) * az;
qxs = ones(ntheta,nphi) * qx;
qys = ones(ntheta,nphi) * qy;
qzs = ones(ntheta,nphi) * qz;
OmegaRFs = ones(ntheta,nphi) * OmegaRF;
rbounces = ones(ntheta,nphi) * rbounce;
tmaxs = ones(ntheta,nphi) * tmax;
rtols = ones(ntheta,nphi) * rtol;
atols = ones(ntheta,nphi) * atol;
maxstepss = ones(ntheta,nphi) * maxsteps;

umin = (1-cos(thetamin))/2; % 0 -> 0
umax = (1-cos(thetamax))/2; % pi/2 -> 1/2
us = linspace(umin,umax,ntheta);
thetas = acos(1-2*us);
phis = linspace(phimin,phimax,nphi);
if positions == "uniform"
    rthetas = repmat(thetas',1,nphi);
    rphis = repmat(phis,ntheta,1);
elseif positions == "random"
    rthetas = acos(1-2*(umin+rand(ntheta,nphi)*(umax-umin)));
    rphis = phimin + rand(ntheta,nphi) * (phimax-phimin);
else
    error("Error, positions must be ""uniform"" or ""random"".");
end

if collisiontype == "head-on"
    v0 = sqrt(3*kB*Tatom/Eh/matom);
    t0s = zeros(ntheta,nphi);
    y01s = zeros(ntheta,nphi);
    y02s = zeros(ntheta,nphi);
    y03s = zeros(ntheta,nphi);
    y04s = zeros(ntheta,nphi);
    y05s = zeros(ntheta,nphi);
    y06s = zeros(ntheta,nphi);
    y07s = r0*sin(rthetas).*cos(rphis);
    y08s = r0*sin(rthetas).*sin(rphis);
    y09s = r0*cos(rthetas);
    y010s = -v0*sin(rthetas).*cos(rphis);
    y011s = -v0*sin(rthetas).*sin(rphis);
    y012s = -v0*cos(rthetas);
elseif collisiontype == "thermal"
    w = 1/2 * OmegaRF * sqrt(abs(ax+1/2*qx^2));
    t0s = zeros(ntheta,nphi);
    y01s = randn(ntheta,nphi)*sqrt(kB*Tion/Eh/mion/w^2);
    y02s = randn(ntheta,nphi)*sqrt(kB*Tion/Eh/mion/w^2);
    y03s = randn(ntheta,nphi)*sqrt(kB*Tion/Eh/mion/w^2);
    y04s = randn(ntheta,nphi)*sqrt(kB*Tion/Eh/mion);
    y05s = randn(ntheta,nphi)*sqrt(kB*Tion/Eh/mion);
    y06s = randn(ntheta,nphi)*sqrt(kB*Tion/Eh/mion);
    y07s = r0*sin(rthetas).*cos(rphis);
    y08s = r0*sin(rthetas).*sin(rphis);
    y09s = r0*cos(rthetas);
    y010s = randn(ntheta,nphi)*sqrt(kB*Tatom/Eh/matom);
    y011s = randn(ntheta,nphi)*sqrt(kB*Tatom/Eh/matom);
    y012s = randn(ntheta,nphi)*sqrt(kB*Tatom/Eh/matom);
    for i = 1:ntheta*nphi
        if y07s(i)*y010s(i) + y08s(i)*y011s(i) + y09s(i)*y012s(i) > 0
            y010s(i) = -y010s(i);
            y011s(i) = -y011s(i);
            y012s(i) = -y012s(i);
        end
    end
else
    error("Error, collisiontype must be ""head-on"" or ""thermal"".");
end

if ncores >= 2
    p = gcp('nocreate');
    if isempty(p)
        p = parpool(ncores);
    elseif p.NumWorkers ~= ncores
        delete(p);
        p = parpool(ncores);
    end
end

tendsp = cell(ncores,1);
yend1sp = cell(ncores,1);
yend2sp = cell(ncores,1);
yend3sp = cell(ncores,1);
yend4sp = cell(ncores,1);
yend5sp = cell(ncores,1);
yend6sp = cell(ncores,1);
yend7sp = cell(ncores,1);
yend8sp = cell(ncores,1);
yend9sp = cell(ncores,1);
yend10sp = cell(ncores,1);
yend11sp = cell(ncores,1);
yend12sp = cell(ncores,1);

anglep = cell(ncores,1);
bouncesp = cell(ncores,1);
lifetimep = cell(ncores,1);
positionp = cell(ncores,1);
transferp = cell(ncores,1);
distp = cell(ncores,1);
nstepsp = cell(ncores,1);
KEp = cell(ncores,1);

mionsp = cell(ncores,1);
matomsp = cell(ncores,1);
nsp = cell(ncores,1);
msp = cell(ncores,1);
Cnsp = cell(ncores,1);
Cmsp = cell(ncores,1);
axsp = cell(ncores,1);
aysp = cell(ncores,1);
azsp = cell(ncores,1);
qxsp = cell(ncores,1);
qysp = cell(ncores,1);
qzsp = cell(ncores,1);
OmegaRFsp = cell(ncores,1);
rbouncesp = cell(ncores,1);
tmaxsp = cell(ncores,1);
rtolsp = cell(ncores,1);
atolsp = cell(ncores,1);
maxstepssp = cell(ncores,1);

t0sp = cell(ncores,1);
y01sp = cell(ncores,1);
y02sp = cell(ncores,1);
y03sp = cell(ncores,1);
y04sp = cell(ncores,1);
y05sp = cell(ncores,1);
y06sp = cell(ncores,1);
y07sp = cell(ncores,1);
y08sp = cell(ncores,1);
y09sp = cell(ncores,1);
y010sp = cell(ncores,1);
y011sp = cell(ncores,1);
y012sp = cell(ncores,1);

if processor == "CPU"
    g = @(x)x;
elseif processor == "GPU"
    g = @(x)gpuArray(x);
else
    error("Error, processor must be ""CPU"" or ""GPU"".");
end

for i = 1:ncores
    mionsp{i} = g(mions(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    matomsp{i} = g(matoms(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    nsp{i} = g(ns(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    msp{i} = g(ms(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    Cnsp{i} = g(Cns(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    Cmsp{i} = g(Cms(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    axsp{i} = g(axs(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    aysp{i} = g(ays(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    azsp{i} = g(azs(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    qxsp{i} = g(qxs(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    qysp{i} = g(qys(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    qzsp{i} = g(qzs(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    OmegaRFsp{i} = g(OmegaRFs(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    rbouncesp{i} = g(rbounces(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    tmaxsp{i} = g(tmaxs(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    rtolsp{i} = g(rtols(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    atolsp{i} = g(atols(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    maxstepssp{i} = g(maxstepss(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    t0sp{i} = g(t0s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    y01sp{i} = g(y01s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    y02sp{i} = g(y02s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    y03sp{i} = g(y03s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    y04sp{i} = g(y04s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    y05sp{i} = g(y05s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    y06sp{i} = g(y06s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    y07sp{i} = g(y07s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    y08sp{i} = g(y08s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    y09sp{i} = g(y09s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    y010sp{i} = g(y010s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    y011sp{i} = g(y011s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
    y012sp{i} = g(y012s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)));
end

fprintf("\nRunning calculations ...\n\n");
tstart = tic;
parfor i = 1:ncores
    [tendsp{i}, yend1sp{i}, yend2sp{i}, yend3sp{i}, yend4sp{i}, yend5sp{i}, yend6sp{i}, ...
        yend7sp{i}, yend8sp{i}, yend9sp{i}, yend10sp{i}, yend11sp{i}, yend12sp{i}, ...
        anglep{i}, bouncesp{i}, lifetimep{i}, positionp{i}, transferp{i}, distp{i}, nstepsp{i}, KEp{i}] = ...
        arrayfun(@ode45gpu, mionsp{i}, matomsp{i}, nsp{i}, msp{i}, Cnsp{i}, Cmsp{i}, ...
        axsp{i}, aysp{i}, azsp{i}, qxsp{i}, qysp{i}, qzsp{i}, OmegaRFsp{i}, rbouncesp{i}, ...
        t0sp{i}, tmaxsp{i}, y01sp{i}, y02sp{i}, y03sp{i}, y04sp{i}, y05sp{i}, y06sp{i}, ...
        y07sp{i}, y08sp{i}, y09sp{i}, y010sp{i}, y011sp{i}, y012sp{i}, rtolsp{i}, atolsp{i}, maxstepssp{i});
end
tend = toc(tstart);
fprintf("\nFinished calculating. It took %.3f seconds.\n\n",tend);

for i = 1:ncores
    tends(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(tendsp{i});
    yend1s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(yend1sp{i});
    yend2s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(yend2sp{i});
    yend3s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(yend3sp{i});
    yend4s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(yend4sp{i});
    yend5s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(yend5sp{i});
    yend6s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(yend6sp{i});
    yend7s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(yend7sp{i});
    yend8s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(yend8sp{i});
    yend9s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(yend9sp{i});
    yend10s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(yend10sp{i});
    yend11s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(yend11sp{i});
    yend12s(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(yend12sp{i});
    angle(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(anglep{i});
    bounces(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(bouncesp{i});
    lifetime(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(lifetimep{i});
    position(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(positionp{i});
    transfer(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(transferp{i});
    dist(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(distp{i});
    nsteps(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(nstepsp{i});
    KE(floor(ntrajectories/ncores*(i-1))+1 : floor(ntrajectories/ncores*i)) = gather(KEp{i});
end

probability = length(bounces(bounces>=2))/length(bounces(:));
fprintf("The probability of forming a complex is %.2f%%.\n\n",100*probability);
file = fopen(main+"/probability.txt",'w');
fprintf(file,"The probability of forming a complex is %.20f.\n\n",probability);
fclose(file);

if longestlived == "yes"
    [~,imax] = max(lifetime(:));
    fprintf("According to ode45gpu on the %s, the longest-lived trajectory had bounces = %d and lifetime = %.1f microseconds.\n", ...
        processor, bounces(imax), lifetime(imax)*tau*1e6);
    if processor == "GPU"
        fprintf("Calculating this trajectory using ode45gpu on the CPU ...\n");
        [tlend,ylend1,ylend2,ylend3,ylend4,ylend5,ylend6,ylend7,ylend8,ylend9,ylend10,ylend11,ylend12,...
            anglel,bouncesl,lifetimel,positionl,transferl,distl,nstepsl,KEl] = ...
            ode45gpu(mion,matom,n,m,Cn,Cm,ax,ay,az,qx,qy,qz,OmegaRF,rbounce,...
            0,tmax,y01s(imax),y02s(imax),y03s(imax),y04s(imax),y05s(imax),y06s(imax), ...
            y07s(imax),y08s(imax),y09s(imax),y010s(imax),y011s(imax),y012s(imax),rtol,atol,maxsteps);
        fprintf("Done. According to ode45gpu on the CPU, bounces = %d and lifetime = %.1f microseconds.\n", bouncesl, lifetimel*tau*1e6);
    end
    fprintf("Calculating this trajectory using ode45 ...\n");
    y0l = [y01s(imax),y02s(imax),y03s(imax),y04s(imax),y05s(imax),y06s(imax),...
        y07s(imax),y08s(imax),y09s(imax),y010s(imax),y011s(imax),y012s(imax)];
    options = odeset('Stats','off','RelTol',rtol,'AbsTol',atol,'OutputFcn',@(t,ruichan,flag)halt(t,ruichan,flag,4*maxsteps));
    [tl, yl] = ode45(@(t,y)f(mion,matom,n,m,Cn,Cm,ax,ay,az,qx,qy,qz,OmegaRF,t,y), [0,tmax], y0l, options);
    rl = sqrt((yl(:,1)-yl(:,7)).^2+(yl(:,2)-yl(:,8)).^2+(yl(:,3)-yl(:,9)).^2);
    if processor == "CPU"
        compare = [yend1s(imax),yend2s(imax),yend3s(imax),yend4s(imax),yend5s(imax),yend6s(imax),...
            yend7s(imax),yend8s(imax),yend9s(imax),yend10s(imax),yend11s(imax),yend12s(imax)];
    elseif processor == "GPU"
        compare = [ylend1,ylend2,ylend3,ylend4,ylend5,ylend6,ylend7,ylend8,ylend9,ylend10,ylend11,ylend12];
    end
    if compare == yl(end,:)
        fprintf("Done. Success, ode45 and ode45gpu match on the CPU!\n\n");
    else
        fprintf("Done. Failure, ode45 and ode45gpu do not match on the CPU.\n\n");
    end
    figure; hold on;
    plot(tl*tau*1e6,rl);
    if ~isempty(tl(rl<rl(1)+100))
        timestamp = tl(rl<rl(1)+100)*tau*1e6;
        timestamp = timestamp(end);
    else
        timestamp = tl(end)*tau*1e6;
    end
    xlim([0,timestamp]);
    ylim([0,2.5*rl(1)]);
    title("Longest-lived trajectory");
    xlabel("Time $t$ ($\mu$s)");
    ylabel("Atom-ion distance $r$ (a.u.)");
    legend("ode45 atom-ion distance $r$");
    saveas(gcf,main+"/trajectorylongest-r.png");
    hold off;
    figure; hold on;
    plot(yl(:,1),yl(:,2));
    plot(yl(:,7),yl(:,8));
    xlim([-2000,2000]);
    ylim([-2000,2000]);
    title("Longest-lived trajectory");
    xlabel("$x$ (Bohr)");
    ylabel("$y$ (Bohr)");
    legend("ion","atom");
    saveas(gcf,main+"/trajectorylongest-xy.png");
    hold off;
    figure; hold on;
    plot(yl(:,1),yl(:,3));
    plot(yl(:,7),yl(:,9));
    xlim([-2000,2000]);
    ylim([-2000,2000]);
    title("Longest-lived trajectory");
    xlabel("$x$ (Bohr)");
    ylabel("$z$ (Bohr)");
    legend("ion","atom");
    saveas(gcf,main+"/trajectorylongest-xz.png");
    hold off;
    figure; hold on;
    plot(yl(:,2),yl(:,3));
    plot(yl(:,8),yl(:,9));
    xlim([-2000,2000]);
    ylim([-2000,2000]);
    title("Longest-lived trajectory");
    xlabel("$y$ (Bohr)");
    ylabel("$z$ (Bohr)");
    legend("ion","atom");
    saveas(gcf,main+"/trajectorylongest-yz.png");
    hold off;
elseif longestlived == "no"
    % do nothing
else
    error("Error, longestlived must be ""yes"" or ""no"".");
end

if custom == "yes"
    fprintf("\nCalculating custom trajectory using ode45gpu ...\n");
    y0c = y0custom;
    [tcend,ycend1,ycend2,ycend3,ycend4,ycend5,ycend6,ycend7,ycend8,ycend9,ycend10,ycend11,ycend12,...
        anglec,bouncesc,lifetimec,positionc,transferc,distc,nstepsc,KEc] = ...
        ode45gpu(mion,matom,n,m,Cn,Cm,ax,ay,az,qx,qy,qz,OmegaRF,rbounce,...
        t0custom,tfcustom,y0c(1),y0c(2),y0c(3),y0c(4),y0c(5),y0c(6),y0c(7),y0c(8),y0c(9),y0c(10),y0c(11),y0c(12),rtol,atol,maxsteps);
    fprintf("Done. According to ode45gpu, bounces = %d and lifetime = %.1f microseconds.\n", bouncesc, lifetimec*tau*1e6);
    fprintf("Calculating custom trajectory using ode45 ...\n");
    options = odeset('Stats','off','RelTol',rtol,'AbsTol',atol,'OutputFcn',@(t,ruichan,flag)halt(t,ruichan,flag,4*maxsteps));
    [tc, yc] = ode45(@(t,y)f(mion,matom,n,m,Cn,Cm,ax,ay,az,qx,qy,qz,OmegaRF,t,y), [t0custom,tfcustom], y0custom, options);
    rc = sqrt((yc(:,1)-yc(:,7)).^2+(yc(:,2)-yc(:,8)).^2+(yc(:,3)-yc(:,9)).^2);
    if [ycend1,ycend2,ycend3,ycend4,ycend5,ycend6,ycend7,ycend8,ycend9,ycend10,ycend11,ycend12] == yc(end,:)
        fprintf("Done. Success, ode45 and ode45gpu match on the CPU!\n\n");
    else
        fprintf("Done. Failure, ode45 and ode45gpu do not match on the CPU.\n\n");
    end
    figure; hold on;
    plot(tc*tau*1e6,rc);
    if ~isempty(tc(rc<rc(1)+100))
        timestamp = tc(rc<rc(1)+100)*tau*1e6;
        timestamp = timestamp(end);
    else
        timestamp = tc(end)*tau*1e6;
    end
    xlim([0,timestamp]);
    ylim([0,2.5*rc(1)]);
    title("Custom trajectory");
    xlabel("Time $t$ ($\mu$s)");
    ylabel("Atom-ion distance $r$ (a.u.)");
    legend("ode45 atom-ion distance $r$");
    saveas(gcf,main+"/trajectorycustom-r.png");
    hold off;
    figure; hold on;
    plot(yc(:,1),yc(:,2));
    plot(yc(:,7),yc(:,8));
    xlim([-2000,2000]);
    ylim([-2000,2000]);
    title("Custom trajectory");
    xlabel("$x$ (Bohr)");
    ylabel("$y$ (Bohr)");
    legend("ion","atom");
    saveas(gcf,main+"/trajectorycustom-xy.png");
    hold off;
    figure; hold on;
    plot(yc(:,1),yc(:,3));
    plot(yc(:,7),yc(:,9));
    xlim([-2000,2000]);
    ylim([-2000,2000]);
    title("Custom trajectory");
    xlabel("$x$ (Bohr)");
    ylabel("$z$ (Bohr)");
    legend("ion","atom");
    saveas(gcf,main+"/trajectorycustom-xz.png");
    hold off;
    figure; hold on;
    plot(yc(:,2),yc(:,3));
    plot(yc(:,8),yc(:,9));
    xlim([-2000,2000]);
    ylim([-2000,2000]);
    title("Custom trajectory");
    xlabel("$y$ (Bohr)");
    ylabel("$z$ (Bohr)");
    legend("ion","atom");
    saveas(gcf,main+"/trajectorycustom-yz.png");
    hold off;
elseif custom == "no"
    % do nothing
else
    error("Error, custom must be ""yes"" or ""no"".");
end

file = fopen(main+"/outliers.txt",'w');
[~,i] = min(angle(:));
fprintf(file,"The trajectory with the smallest scattering angle had angle %.20e radians, starting from\n" + ...
    "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
    "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",angle(i), ...
    y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
[~,i] = max(angle(:));
fprintf(file,"The trajectory with the largest scattering angle had angle %.20e radians, starting from\n" + ...
    "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
    "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",angle(i), ...
    y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
[~,i] = min(bounces(:));
fprintf(file,"The trajectory with the least bounces had %d bounces, starting from\n" + ...
    "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
    "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",bounces(i), ...
    y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
[~,i] = max(bounces(:));
fprintf(file,"The trajectory with the most bounces had %d bounces, starting from\n" + ...
    "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
    "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",bounces(i), ...
    y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
if max(lifetime(:)) > 0
    lt = lifetime;
    for i = 1:ntrajectories
        if lt(i) == 0
            lt(i) = 1e300;
        end
    end
    [~,i] = min(lt(:));
    fprintf(file,"The trajectory with the smallest complex lifetime had lifetime %.20f microseconds, starting from\n" + ...
        "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
        "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",lifetime(i)*tau*1e6, ...
        y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
    [~,i] = max(lifetime(:));
    fprintf(file,"The trajectory with the largest complex lifetime had lifetime %.20f microseconds, starting from\n" + ...
        "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
        "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",lifetime(i)*tau*1e6, ...
        y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
else
    fprintf(file,"No trajectories formed a complex.\n");
end
[~,i] = min(position(:));
fprintf(file,"The trajectory with the smallest ion position at first bounce had position %.20f Bohr, starting from\n" + ...
    "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
    "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",position(i), ...
    y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
[~,i] = max(position(:));
fprintf(file,"The trajectory with the largest ion position at first bounce had position %.20f Bohr, starting from\n" + ...
    "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
    "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",position(i), ...
    y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
[~,i] = min(transfer(:));
fprintf(file,"The trajectory with the smallest momentum transfer had transfer %.20e a.u., starting from\n" + ...
    "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
    "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",transfer(i), ...
    y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
[~,i] = max(transfer(:));
fprintf(file,"The trajectory with the largest momentum transfer had transfer %.20e a.u., starting from\n" + ...
    "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
    "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",transfer(i), ...
    y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
[~,i] = min(dist(:));
fprintf(file,"The trajectory with the smallest final atom-ion distance had distance %.20f Bohr, starting from\n" + ...
    "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
    "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",dist(i), ...
    y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
[~,i] = max(dist(:));
fprintf(file,"The trajectory with the largest final atom-ion distance had distance %.20f Bohr, starting from\n" + ...
    "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
    "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",dist(i), ...
    y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
[~,i] = min(nsteps(:));
fprintf(file,"The trajectory with the smallest number of timesteps had %d timesteps, starting from\n" + ...
    "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
    "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",nsteps(i), ...
    y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
[~,i] = max(nsteps(:));
fprintf(file,"The trajectory with the largest number of timesteps had %d timesteps, starting from\n" + ...
    "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
    "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",nsteps(i), ...
    y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
[~,i] = min(KE(:));
fprintf(file,"The trajectory with the smallest final ion kinetic energy had KE %.20e Hartree, starting from\n" + ...
    "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
    "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",KE(i), ...
    y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
[~,i] = max(KE(:));
fprintf(file,"The trajectory with the largest final ion kinetic energy had KE %.20e Hartree, starting from\n" + ...
    "y0custom = [%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n" + ...
    "%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e,...\n%.20e];\n\n",KE(i), ...
    y01s(i),y02s(i),y03s(i),y04s(i),y05s(i),y06s(i),y07s(i),y08s(i),y09s(i),y010s(i),y011s(i),y012s(i));
fclose(file);

if savecsvs == "yes"
    writematrix(rthetas,main+"/athetas.csv");
    writematrix(rphis,main+"/aphis.csv");
end

color = jet(1024);
v = linspace(thetamin,thetamax,4);

%angle
figure; hold on;
imagesc(us,phis,angle');
axis([umin,umax,phimin,phimax]);
set(gca,'YDir','normal');
xlabel("$\theta$ (radians)");
ylabel("$\phi$ (radians)");
xticks(0.5*(1-cos(v)));
xticklabels({sprintf("%.1f",v(1)),sprintf("%.1f",v(2)),sprintf("%.1f",v(3)),sprintf("%.1f",v(4))});
colormap(jet);
c = colorbar;
c.Label.Interpreter = "Latex";
c.Label.String = "Scattering angle (radians)";
c.TickLabelInterpreter = "Latex";
saveas(gcf,main+'/heatmapangle.png');
if savecsvs == "yes" || onlyonecsv == "angle"
    writematrix(angle,main+"/dataangle.csv");
end
hold off;

%bounces
figure; hold on;
imagesc(us,phis,bounces');
axis([umin,umax,phimin,phimax]);
set(gca,'YDir','normal');
xlabel("$\theta$ (radians)");
ylabel("$\phi$ (radians)");
xticks(0.5*(1-cos(v)));
xticklabels({sprintf("%.1f",v(1)),sprintf("%.1f",v(2)),sprintf("%.1f",v(3)),sprintf("%.1f",v(4))});
colormap(jet);
c = colorbar;
c.Label.Interpreter = "Latex";
c.Label.String = "Number of bounces";
c.TickLabelInterpreter = "Latex";
saveas(gcf,main+'/heatmapbounces.png');
if savecsvs == "yes" || onlyonecsv == "bounces"
    writematrix(bounces,main+"/databounces.csv");
end
hold off;

%lifetime
figure; hold on;
imagesc(us,phis,lifetime');
axis([umin,umax,phimin,phimax]);
set(gca,'YDir','normal');
xlabel("$\theta$ (radians)");
ylabel("$\phi$ (radians)");
xticks(0.5*(1-cos(v)));
xticklabels({sprintf("%.1f",v(1)),sprintf("%.1f",v(2)),sprintf("%.1f",v(3)),sprintf("%.1f",v(4))});
colormap(jet);
c = colorbar;
c.Label.Interpreter = "Latex";
c.Label.String = "Complex lifetime (a.u.)";
c.TickLabelInterpreter = "Latex";
saveas(gcf,main+'/heatmaplifetime.png');
if savecsvs == "yes" || onlyonecsv == "lifetime"
    writematrix(lifetime,main+"/datalifetime.csv");
end
hold off;

%position
figure; hold on;
imagesc(us,phis,position');
axis([umin,umax,phimin,phimax]);
set(gca,'YDir','normal');
xlabel("$\theta$ (radians)");
ylabel("$\phi$ (radians)");
xticks(0.5*(1-cos(v)));
xticklabels({sprintf("%.1f",v(1)),sprintf("%.1f",v(2)),sprintf("%.1f",v(3)),sprintf("%.1f",v(4))});
colormap(jet);
c = colorbar;
c.Label.Interpreter = "Latex";
c.Label.String = "Ion position at first bounce (Bohr)";
c.TickLabelInterpreter = "Latex";
saveas(gcf,main+'/heatmapposition.png');
if savecsvs == "yes" || onlyonecsv == "position"
    writematrix(position,main+"/dataposition.csv");
end
hold off;

%transfer
figure; hold on;
imagesc(us,phis,transfer');
axis([umin,umax,phimin,phimax]);
set(gca,'YDir','normal');
xlabel("$\theta$ (radians)");
ylabel("$\phi$ (radians)");
xticks(0.5*(1-cos(v)));
xticklabels({sprintf("%.1f",v(1)),sprintf("%.1f",v(2)),sprintf("%.1f",v(3)),sprintf("%.1f",v(4))});
colormap(jet);
c = colorbar;
c.Label.Interpreter = "Latex";
c.Label.String = "Momentum transfer (a.u.)";
c.TickLabelInterpreter = "Latex";
saveas(gcf,main+'/heatmaptransfer.png');
if savecsvs == "yes" || onlyonecsv == "transfer"
    writematrix(transfer,main+"/datatransfer.csv");
end
hold off;

%dist
figure; hold on;
imagesc(us,phis,dist');
axis([umin,umax,phimin,phimax]);
set(gca,'YDir','normal');
xlabel("$\theta$ (radians)");
ylabel("$\phi$ (radians)");
xticks(0.5*(1-cos(v)));
xticklabels({sprintf("%.1f",v(1)),sprintf("%.1f",v(2)),sprintf("%.1f",v(3)),sprintf("%.1f",v(4))});
colormap(jet);
c = colorbar;
c.Label.Interpreter = "Latex";
c.Label.String = "Final atom-ion distance (Bohr)";
c.TickLabelInterpreter = "Latex";
saveas(gcf,main+'/heatmapdist.png');
if savecsvs == "yes" || onlyonecsv == "dist"
    writematrix(dist,main+"/datadist.csv");
end
hold off;

%nsteps
figure; hold on;
imagesc(us,phis,nsteps');
axis([umin,umax,phimin,phimax]);
set(gca,'YDir','normal');
xlabel("$\theta$ (radians)");
ylabel("$\phi$ (radians)");
xticks(0.5*(1-cos(v)));
xticklabels({sprintf("%.1f",v(1)),sprintf("%.1f",v(2)),sprintf("%.1f",v(3)),sprintf("%.1f",v(4))});
colormap(jet);
c = colorbar;
c.Label.Interpreter = "Latex";
c.Label.String = "Number of timesteps";
c.TickLabelInterpreter = "Latex";
saveas(gcf,main+'/heatmapnsteps.png');
if savecsvs == "yes" || onlyonecsv == "nsteps"
    writematrix(nsteps,main+"/datansteps.csv");
end
hold off;

%KE
figure; hold on;
imagesc(us,phis,KE');
axis([umin,umax,phimin,phimax]);
set(gca,'YDir','normal');
xlabel("$\theta$ (radians)");
ylabel("$\phi$ (radians)");
xticks(0.5*(1-cos(v)));
xticklabels({sprintf("%.1f",v(1)),sprintf("%.1f",v(2)),sprintf("%.1f",v(3)),sprintf("%.1f",v(4))});
colormap(jet);
c = colorbar;
c.Label.Interpreter = "Latex";
c.Label.String = "Final ion kinetic energy (Hartree)";
c.TickLabelInterpreter = "Latex";
saveas(gcf,main+'/heatmapKE.png');
if savecsvs == "yes" || onlyonecsv == "KE"
    writematrix(KE,main+"/dataKE.csv");
end
hold off;

clear tendsp yend1sp yend2sp yend3sp yend4sp yend5sp yend6sp yend7sp yend8sp yend9sp yend10sp yend11sp yend12sp ...
    anglep bouncesp lifetimep positionp transferp distp nstepsp KEp rthetasp rphisp ...
    atolsp axsp aysp azsp Cmsp Cnsp matomsp maxstepssp mionsp msp nsp OmegaRFsp qxsp qysp qzsp rbouncesp rtolsp t0sp tmaxsp ...
    y01sp y02sp y03sp y04sp y05sp y06sp y07sp y08sp y09sp y010sp y011sp y012sp ...
    atols axs ays azs Cms Cns matoms maxstepss mions ms ns OmegaRFs qxs qys qzs rbounces rtols t0s tmaxs ...
    tl yl rl tc yc rc
if saveworkspace == "yes"
    save(main+"/work.mat");
end
fprintf("The results are in the folder %s.\n\n", main);



function status = halt(t,~,flag,maxsteps)
    status = 0;
    persistent nstep
    switch flag
        case 'init'
            nstep = 0;
        case []
            nstep = nstep + length(t);
            if nstep >= maxsteps
                status = 1;
            end
    end
end

function f = f(mion,matom,n,m,Cn,Cm,ax,ay,az,qx,qy,qz,OmegaRF,t,y)
    f = zeros(12,1);
    f(1:3) = y(4:6);
    f(7:9) = y(10:12);
    r = sqrt((y(1)-y(7))^2+(y(2)-y(8))^2+(y(3)-y(9))^2);
    dVdr = n*Cn/r^(n+1) - m*Cm/r^(m+1);
    f(4:6) = -1/mion*dVdr*(y(1:3)-y(7:9))/r;
    f(10:12) = -1/matom*dVdr*(y(7:9)-y(1:3))/r;
    f(4) = f(4) - (ax+2*qx*cos(OmegaRF*t))*OmegaRF^2/4*y(1);
    f(5) = f(5) - (ay+2*qy*cos(OmegaRF*t))*OmegaRF^2/4*y(2);
    f(6) = f(6) - (az+2*qz*cos(OmegaRF*t))*OmegaRF^2/4*y(3);
end

function [tend,yend1,yend2,yend3,yend4,yend5,yend6,yend7,yend8,yend9,yend10,yend11,yend12,...
    angle,bounces,lifetime,position,transfer,dist,nsteps,KE] = ...
    ode45gpu(mion,matom,n,m,Cn,Cm,ax,ay,az,qx,qy,qz,OmegaRF,rbounce,...
    t0,tf,y01,y02,y03,y04,y05,y06,y07,y08,y09,y010,y011,y012,rtol,atol,maxsteps)
    
    bounces = 0;
    bouncing = false;
    tfirst = 0;
    tlast = 0;
    tmin = 0;
    rmin = 1e300;
    lifetime = 0;
    pmin = 0;
    position = 0;
    nsteps = 0;
    f01 = y04; f02 = y05; f03 = y06;                       %ode start
    f07 = y010; f08 = y011; f09 = y012;
    r = sqrt((y01-y07)^2+(y02-y08)^2+(y03-y09)^2);
    dVdr = n*Cn/r^(n+1) - m*Cm/r^(m+1);
    f04 = -1/mion*dVdr*(y01-y07)/r;
    f05 = -1/mion*dVdr*(y02-y08)/r;
    f06 = -1/mion*dVdr*(y03-y09)/r;
    f010 = -1/matom*dVdr*(y07-y01)/r;
    f011 = -1/matom*dVdr*(y08-y02)/r;
    f012 = -1/matom*dVdr*(y09-y03)/r;
    f04 = f04 - (ax+2*qx*cos(OmegaRF*t0))*OmegaRF^2/4*y01;
    f05 = f05 - (ay+2*qy*cos(OmegaRF*t0))*OmegaRF^2/4*y02;
    f06 = f06 - (az+2*qz*cos(OmegaRF*t0))*OmegaRF^2/4*y03; %ode end
    threshold = atol ./ rtol;
    hmax = min(tf-t0,max(0.1*(tf-t0),16.0*eps*tf));
    t = t0;
    y11 = y01; y12 = y02; y13 = y03; y14 = y04; y15 = y05; y16 = y06;
    y17 = y07; y18 = y08; y19 = y09; y110 = y010; y111 = y011; y112 = y012;
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
    absh = min(hmax, tf - t0);
    rh =         abs(f01 / max(abs(y01),threshold));
    rh = max(rh, abs(f02 / max(abs(y02),threshold)));
    rh = max(rh, abs(f03 / max(abs(y03),threshold)));
    rh = max(rh, abs(f04 / max(abs(y04),threshold)));
    rh = max(rh, abs(f05 / max(abs(y05),threshold)));
    rh = max(rh, abs(f06 / max(abs(y06),threshold)));
    rh = max(rh, abs(f07 / max(abs(y07),threshold)));
    rh = max(rh, abs(f08 / max(abs(y08),threshold)));
    rh = max(rh, abs(f09 / max(abs(y09),threshold)));
    rh = max(rh, abs(f010 / max(abs(y010),threshold)));
    rh = max(rh, abs(f011 / max(abs(y011),threshold)));
    rh = max(rh, abs(f012 / max(abs(y012),threshold)));
    rh = rh / (0.8 * rtol^pow);
    if absh * rh > 1
        absh = 1 / rh;
    end
    absh = max(absh, hmin);
    f11 = f01; f12 = f02; f13 = f03; f14 = f04; f15 = f05; f16 = f06;
    f17 = f07; f18 = f08; f19 = f09; f110 = f010; f111 = f011; f112 = f012;
    t2 = t;
    y21 = y11; y22 = y12; y23 = y13; y24 = y14; y25 = y15; y26 = y16;
    y27 = y17; y28 = y18; y29 = y19; y210 = y110; y211 = y111; y212 = y112;
    f71 = f11; f72 = f12; f73 = f13; f74 = f14; f75 = f15; f76 = f16;
    f77 = f17; f78 = f18; f79 = f19; f710 = f110; f711 = f111; f712 = f112;
    err = 0;
    
    
    
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
            y21 = y11 + h .* (b11.*f11 );
            y22 = y12 + h .* (b11.*f12 );
            y23 = y13 + h .* (b11.*f13 );
            y24 = y14 + h .* (b11.*f14 );
            y25 = y15 + h .* (b11.*f15 );
            y26 = y16 + h .* (b11.*f16 );
            y27 = y17 + h .* (b11.*f17 );
            y28 = y18 + h .* (b11.*f18 );
            y29 = y19 + h .* (b11.*f19 );
            y210 = y110 + h .* (b11.*f110 );
            y211 = y111 + h .* (b11.*f111 );
            y212 = y112 + h .* (b11.*f112 );
            t2 = t + h .* a2;
            f21 = y24; f22 = y25; f23 = y26;                       %ode start
            f27 = y210; f28 = y211; f29 = y212;
            r = sqrt((y21-y27)^2+(y22-y28)^2+(y23-y29)^2);
            dVdr = n*Cn/r^(n+1) - m*Cm/r^(m+1);
            f24 = -1/mion*dVdr*(y21-y27)/r;
            f25 = -1/mion*dVdr*(y22-y28)/r;
            f26 = -1/mion*dVdr*(y23-y29)/r;
            f210 = -1/matom*dVdr*(y27-y21)/r;
            f211 = -1/matom*dVdr*(y28-y22)/r;
            f212 = -1/matom*dVdr*(y29-y23)/r;
            f24 = f24 - (ax+2*qx*cos(OmegaRF*t2))*OmegaRF^2/4*y21;
            f25 = f25 - (ay+2*qy*cos(OmegaRF*t2))*OmegaRF^2/4*y22;
            f26 = f26 - (az+2*qz*cos(OmegaRF*t2))*OmegaRF^2/4*y23; %ode end
            y21 = y11 + h .* (b21.*f11 + b22.*f21 );
            y22 = y12 + h .* (b21.*f12 + b22.*f22 );
            y23 = y13 + h .* (b21.*f13 + b22.*f23 );
            y24 = y14 + h .* (b21.*f14 + b22.*f24 );
            y25 = y15 + h .* (b21.*f15 + b22.*f25 );
            y26 = y16 + h .* (b21.*f16 + b22.*f26 );
            y27 = y17 + h .* (b21.*f17 + b22.*f27 );
            y28 = y18 + h .* (b21.*f18 + b22.*f28 );
            y29 = y19 + h .* (b21.*f19 + b22.*f29 );
            y210 = y110 + h .* (b21.*f110 + b22.*f210 );
            y211 = y111 + h .* (b21.*f111 + b22.*f211 );
            y212 = y112 + h .* (b21.*f112 + b22.*f212 );
            t2 = t + h .* a3;
            f31 = y24; f32 = y25; f33 = y26;                       %ode start
            f37 = y210; f38 = y211; f39 = y212;
            r = sqrt((y21-y27)^2+(y22-y28)^2+(y23-y29)^2);
            dVdr = n*Cn/r^(n+1) - m*Cm/r^(m+1);
            f34 = -1/mion*dVdr*(y21-y27)/r;
            f35 = -1/mion*dVdr*(y22-y28)/r;
            f36 = -1/mion*dVdr*(y23-y29)/r;
            f310 = -1/matom*dVdr*(y27-y21)/r;
            f311 = -1/matom*dVdr*(y28-y22)/r;
            f312 = -1/matom*dVdr*(y29-y23)/r;
            f34 = f34 - (ax+2*qx*cos(OmegaRF*t2))*OmegaRF^2/4*y21;
            f35 = f35 - (ay+2*qy*cos(OmegaRF*t2))*OmegaRF^2/4*y22;
            f36 = f36 - (az+2*qz*cos(OmegaRF*t2))*OmegaRF^2/4*y23; %ode end
            y21 = y11 + h .* (b31.*f11 + b32.*f21 + b33.*f31 );
            y22 = y12 + h .* (b31.*f12 + b32.*f22 + b33.*f32 );
            y23 = y13 + h .* (b31.*f13 + b32.*f23 + b33.*f33 );
            y24 = y14 + h .* (b31.*f14 + b32.*f24 + b33.*f34 );
            y25 = y15 + h .* (b31.*f15 + b32.*f25 + b33.*f35 );
            y26 = y16 + h .* (b31.*f16 + b32.*f26 + b33.*f36 );
            y27 = y17 + h .* (b31.*f17 + b32.*f27 + b33.*f37 );
            y28 = y18 + h .* (b31.*f18 + b32.*f28 + b33.*f38 );
            y29 = y19 + h .* (b31.*f19 + b32.*f29 + b33.*f39 );
            y210 = y110 + h .* (b31.*f110 + b32.*f210 + b33.*f310 );
            y211 = y111 + h .* (b31.*f111 + b32.*f211 + b33.*f311 );
            y212 = y112 + h .* (b31.*f112 + b32.*f212 + b33.*f312 );
            t2 = t + h .* a4;
            f41 = y24; f42 = y25; f43 = y26;                       %ode start
            f47 = y210; f48 = y211; f49 = y212;
            r = sqrt((y21-y27)^2+(y22-y28)^2+(y23-y29)^2);
            dVdr = n*Cn/r^(n+1) - m*Cm/r^(m+1);
            f44 = -1/mion*dVdr*(y21-y27)/r;
            f45 = -1/mion*dVdr*(y22-y28)/r;
            f46 = -1/mion*dVdr*(y23-y29)/r;
            f410 = -1/matom*dVdr*(y27-y21)/r;
            f411 = -1/matom*dVdr*(y28-y22)/r;
            f412 = -1/matom*dVdr*(y29-y23)/r;
            f44 = f44 - (ax+2*qx*cos(OmegaRF*t2))*OmegaRF^2/4*y21;
            f45 = f45 - (ay+2*qy*cos(OmegaRF*t2))*OmegaRF^2/4*y22;
            f46 = f46 - (az+2*qz*cos(OmegaRF*t2))*OmegaRF^2/4*y23; %ode end
            y21 = y11 + h .* (b41.*f11 + b42.*f21 + b43.*f31 + b44.*f41 );
            y22 = y12 + h .* (b41.*f12 + b42.*f22 + b43.*f32 + b44.*f42 );
            y23 = y13 + h .* (b41.*f13 + b42.*f23 + b43.*f33 + b44.*f43 );
            y24 = y14 + h .* (b41.*f14 + b42.*f24 + b43.*f34 + b44.*f44 );
            y25 = y15 + h .* (b41.*f15 + b42.*f25 + b43.*f35 + b44.*f45 );
            y26 = y16 + h .* (b41.*f16 + b42.*f26 + b43.*f36 + b44.*f46 );
            y27 = y17 + h .* (b41.*f17 + b42.*f27 + b43.*f37 + b44.*f47 );
            y28 = y18 + h .* (b41.*f18 + b42.*f28 + b43.*f38 + b44.*f48 );
            y29 = y19 + h .* (b41.*f19 + b42.*f29 + b43.*f39 + b44.*f49 );
            y210 = y110 + h .* (b41.*f110 + b42.*f210 + b43.*f310 + b44.*f410 );
            y211 = y111 + h .* (b41.*f111 + b42.*f211 + b43.*f311 + b44.*f411 );
            y212 = y112 + h .* (b41.*f112 + b42.*f212 + b43.*f312 + b44.*f412 );
            t2 = t + h .* a5;
            f51 = y24; f52 = y25; f53 = y26;                       %ode start
            f57 = y210; f58 = y211; f59 = y212;
            r = sqrt((y21-y27)^2+(y22-y28)^2+(y23-y29)^2);
            dVdr = n*Cn/r^(n+1) - m*Cm/r^(m+1);
            f54 = -1/mion*dVdr*(y21-y27)/r;
            f55 = -1/mion*dVdr*(y22-y28)/r;
            f56 = -1/mion*dVdr*(y23-y29)/r;
            f510 = -1/matom*dVdr*(y27-y21)/r;
            f511 = -1/matom*dVdr*(y28-y22)/r;
            f512 = -1/matom*dVdr*(y29-y23)/r;
            f54 = f54 - (ax+2*qx*cos(OmegaRF*t2))*OmegaRF^2/4*y21;
            f55 = f55 - (ay+2*qy*cos(OmegaRF*t2))*OmegaRF^2/4*y22;
            f56 = f56 - (az+2*qz*cos(OmegaRF*t2))*OmegaRF^2/4*y23; %ode end
            y21 = y11 + h .* (b51.*f11 + b52.*f21 + b53.*f31 + b54.*f41 + b55.*f51 );
            y22 = y12 + h .* (b51.*f12 + b52.*f22 + b53.*f32 + b54.*f42 + b55.*f52 );
            y23 = y13 + h .* (b51.*f13 + b52.*f23 + b53.*f33 + b54.*f43 + b55.*f53 );
            y24 = y14 + h .* (b51.*f14 + b52.*f24 + b53.*f34 + b54.*f44 + b55.*f54 );
            y25 = y15 + h .* (b51.*f15 + b52.*f25 + b53.*f35 + b54.*f45 + b55.*f55 );
            y26 = y16 + h .* (b51.*f16 + b52.*f26 + b53.*f36 + b54.*f46 + b55.*f56 );
            y27 = y17 + h .* (b51.*f17 + b52.*f27 + b53.*f37 + b54.*f47 + b55.*f57 );
            y28 = y18 + h .* (b51.*f18 + b52.*f28 + b53.*f38 + b54.*f48 + b55.*f58 );
            y29 = y19 + h .* (b51.*f19 + b52.*f29 + b53.*f39 + b54.*f49 + b55.*f59 );
            y210 = y110 + h .* (b51.*f110 + b52.*f210 + b53.*f310 + b54.*f410 + b55.*f510 );
            y211 = y111 + h .* (b51.*f111 + b52.*f211 + b53.*f311 + b54.*f411 + b55.*f511 );
            y212 = y112 + h .* (b51.*f112 + b52.*f212 + b53.*f312 + b54.*f412 + b55.*f512 );
            t2 = t + h;
            f61 = y24; f62 = y25; f63 = y26;                       %ode start
            f67 = y210; f68 = y211; f69 = y212;
            r = sqrt((y21-y27)^2+(y22-y28)^2+(y23-y29)^2);
            dVdr = n*Cn/r^(n+1) - m*Cm/r^(m+1);
            f64 = -1/mion*dVdr*(y21-y27)/r;
            f65 = -1/mion*dVdr*(y22-y28)/r;
            f66 = -1/mion*dVdr*(y23-y29)/r;
            f610 = -1/matom*dVdr*(y27-y21)/r;
            f611 = -1/matom*dVdr*(y28-y22)/r;
            f612 = -1/matom*dVdr*(y29-y23)/r;
            f64 = f64 - (ax+2*qx*cos(OmegaRF*t2))*OmegaRF^2/4*y21;
            f65 = f65 - (ay+2*qy*cos(OmegaRF*t2))*OmegaRF^2/4*y22;
            f66 = f66 - (az+2*qz*cos(OmegaRF*t2))*OmegaRF^2/4*y23; %ode end
            t2 = t + h;
            if done
                t2 = tf;
            end
            h = t2 - t;
            y21 = y11 + h.* ( b61.*f11 + b63.*f31 + b64.*f41 + b65.*f51 + b66.*f61 );
            y22 = y12 + h.* ( b61.*f12 + b63.*f32 + b64.*f42 + b65.*f52 + b66.*f62 );
            y23 = y13 + h.* ( b61.*f13 + b63.*f33 + b64.*f43 + b65.*f53 + b66.*f63 );
            y24 = y14 + h.* ( b61.*f14 + b63.*f34 + b64.*f44 + b65.*f54 + b66.*f64 );
            y25 = y15 + h.* ( b61.*f15 + b63.*f35 + b64.*f45 + b65.*f55 + b66.*f65 );
            y26 = y16 + h.* ( b61.*f16 + b63.*f36 + b64.*f46 + b65.*f56 + b66.*f66 );
            y27 = y17 + h.* ( b61.*f17 + b63.*f37 + b64.*f47 + b65.*f57 + b66.*f67 );
            y28 = y18 + h.* ( b61.*f18 + b63.*f38 + b64.*f48 + b65.*f58 + b66.*f68 );
            y29 = y19 + h.* ( b61.*f19 + b63.*f39 + b64.*f49 + b65.*f59 + b66.*f69 );
            y210 = y110 + h.* ( b61.*f110 + b63.*f310 + b64.*f410 + b65.*f510 + b66.*f610 );
            y211 = y111 + h.* ( b61.*f111 + b63.*f311 + b64.*f411 + b65.*f511 + b66.*f611 );
            y212 = y112 + h.* ( b61.*f112 + b63.*f312 + b64.*f412 + b65.*f512 + b66.*f612 );
            f71 = y24; f72 = y25; f73 = y26;                     %ode start
            f77 = y210; f78 = y211; f79 = y212;
            r = sqrt((y21-y27)^2+(y22-y28)^2+(y23-y29)^2);
            dVdr = n*Cn/r^(n+1) - m*Cm/r^(m+1);
            f74 = -1/mion*dVdr*(y21-y27)/r;
            f75 = -1/mion*dVdr*(y22-y28)/r;
            f76 = -1/mion*dVdr*(y23-y29)/r;
            f710 = -1/matom*dVdr*(y27-y21)/r;
            f711 = -1/matom*dVdr*(y28-y22)/r;
            f712 = -1/matom*dVdr*(y29-y23)/r;
            f74 = f74 - (ax+2*qx*cos(OmegaRF*t2))*OmegaRF^2/4*y21;
            f75 = f75 - (ay+2*qy*cos(OmegaRF*t2))*OmegaRF^2/4*y22;
            f76 = f76 - (az+2*qz*cos(OmegaRF*t2))*OmegaRF^2/4*y23; %ode end
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
            err =          abs((fE1) ./ max(max(abs(y11),abs(y21)),threshold));
            err = max(err, abs((fE2) ./ max(max(abs(y12),abs(y22)),threshold)));
            err = max(err, abs((fE3) ./ max(max(abs(y13),abs(y23)),threshold)));
            err = max(err, abs((fE4) ./ max(max(abs(y14),abs(y24)),threshold)));
            err = max(err, abs((fE5) ./ max(max(abs(y15),abs(y25)),threshold)));
            err = max(err, abs((fE6) ./ max(max(abs(y16),abs(y26)),threshold)));
            err = max(err, abs((fE7) ./ max(max(abs(y17),abs(y27)),threshold)));
            err = max(err, abs((fE8) ./ max(max(abs(y18),abs(y28)),threshold)));
            err = max(err, abs((fE9) ./ max(max(abs(y19),abs(y29)),threshold)));
            err = max(err, abs((fE10) ./ max(max(abs(y110),abs(y210)),threshold)));
            err = max(err, abs((fE11) ./ max(max(abs(y111),abs(y211)),threshold)));
            err = max(err, abs((fE12) ./ max(max(abs(y112),abs(y212)),threshold)));
            err = err * absh;
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
        if nsteps >= maxsteps
            done = true;
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
        t = t2;
        y11 = y21; y12 = y22; y13 = y23; y14 = y24; y15 = y25; y16 = y26;
        y17 = y27; y18 = y28; y19 = y29; y110 = y210; y111 = y211; y112 = y212;
        f11 = f71; f12 = f72; f13 = f73; f14 = f74; f15 = f75; f16 = f76;
        f17 = f77; f18 = f78; f19 = f79; f110 = f710; f111 = f711; f112 = f712;
        if bouncing
            if r < rmin
                tmin = t;
                rmin = r;
                pmin = sqrt(y21^2+y22^2+y23^2);
            end
            if r >= rbounce
                bouncing = false;
                bounces = bounces + 1;
                if bounces == 1
                    tfirst = tmin;
                    position = pmin;
                else
                    tlast = tmin;
                end
                rmin = 1e300;
            end
        elseif ~bouncing 
            if r < rbounce
                bouncing = true;
            end
        end
    end
    tend = t2;
    yend1 = y21; yend2 = y22; yend3 = y23; yend4 = y24; yend5 = y25; yend6 = y26;
    yend7 = y27; yend8 = y28; yend9 = y29; yend10 = y210; yend11 = y211; yend12 = y212;
    adotb = y010*yend10+y011*yend11+y012*yend12;
    norma = sqrt(y010^2+y011^2+y012^2);
    normb = sqrt(yend10^2+yend11^2+yend12^2);
    angle = real(acos(complex(adotb/(norma*normb))));
    if bounces >= 2
        lifetime = tlast - tfirst;
    end
    transfer = matom * sqrt((y010-yend10)^2+(y011-yend11)^2+(y012-yend12)^2);
    dist = r;
    KE = 1/2 * mion * (yend4^2+yend5^2+yend6^2);
end