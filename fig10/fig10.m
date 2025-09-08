clear all;

% Input start (units are a.u., unless otherwise specified):

mion          = 170.936323 * 1822.88848628276014097;
matom         = 86.9091805310 * 1822.88848628276014097;
Tion          = 10.0e-6;   % Kelvin
Tatom         = 1.0e-6;    % Kelvin
potential     = "CnDe";    % "CnCm" or "DeRe" or "CnDe"
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
ntrajectories = 10004569;  % unitless
ntheta        = 3163;      % unitless
nphi          = 3163;      % unitless
thetamin      = 0;         % radians
thetamax      = pi;        % radians
phimin        = 0;         % radians
phimax        = 2*pi;      % radians
rtol          = 1e-10;     % unitless
atol          = 1e-20;     % unitless
maxsteps      = 1000000;   % unitless
ncores        = 8;         % unitless
positions     = "uniform"; % "uniform" or "random"
collisionType = "head-on"; % "head-on" or "thermal"
processor     = "GPU";     % "GPU" or "CPU"

longestlived  = "yes";     % "yes" or "no"
custom        = "yes";     % "yes" or "no"
t0custom      = 0;
tfcustom      = tmax;
y0custom      = [0,0,0,0,0,0,4154,938,2620,...
                -6.4e-9,-1.45e-9,-4.1e-9];

saveworkspace = "no";      % "yes" or "no"
savecsvs      = "no";      % "yes" or "no"
onlyonecsv    = "angle";   % "no" or "angle" or ...

% Input end

