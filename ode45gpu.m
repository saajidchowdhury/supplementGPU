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