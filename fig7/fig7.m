%initial conditions to parallelize
ntheta = 50;
nphi = 50;
thetas = acos(1-2*linspace(0,1/2,ntheta));
phis = linspace(0,pi/2,nphi);
thetamatrix = repmat(thetas',1,nphi);
phimatrix = repmat(phis,ntheta,1);
...
x0s = gpuArray(r0*sin(thetamatrix).*cos(phimatrix));
...

%constants, tolerances
C4s = gpuArray(ones(ntheta,nphi) * C4);
...
rtols = gpuArray(ones(ntheta,nphi) * rtol);
atols = gpuArray(ones(ntheta,nphi) * atol);

%call ode45gpu to operate on gpuArrays
[tend,yend1,yend2,yend3,yend4,yend5,yend6,yend7, ...
 yend8,yend9,yend10,yend11,yend12,other] = ...
 arrayfun(@ode45gpu, C4s,C8s,mions,matoms, ...
 axs,ays,azs,qxs,qys,qzs,OmegaRFs,t0s,tmaxs, ...
 xion0s,yion0s,zion0s,vxion0s,vyion0s,vzion0s, ...
 x0s,y0s,z0s,vx0s,vy0s,vz0s,rtols,atols);

%gather outputs
yend1 = gather(yend1); yend2 = gather(yend2);
...


