function [tend,yend1,yend2,yend3,yend4,yend5,yend6, ...
    yend7,yend8,yend9,yend10,yend11,yend12,other] = ...
    ode45gpu(C4,C8,mion,matom,ax,ay,az,qx,qy,qz,OmegaRF, ...
    t0,tf,y01,y02,y03,y04,y05,y06,y07,y08,y09,y010,y011,y012,rtol,atol)
    ...
    while ~done
        ...
        while true
            ...
            f21 = y24; f22 = y25; f23 = y26;                     %ode start
            f27 = y210; f28 = y211; f29 = y212;
            r = sqrt((y21-y27)^2+(y22-y28)^2+(y23-y29)^2);
            dVdr = 4*C4/r^5 - 8*C8/r^9;
            f24 = -1/mion*dVdr*(y21-y27)/r;
            f25 = -1/mion*dVdr*(y22-y28)/r;
            f26 = -1/mion*dVdr*(y23-y29)/r;
            f210 = -1/matom*dVdr*(y27-y21)/r;
            f211 = -1/matom*dVdr*(y28-y22)/r;
            f212 = -1/matom*dVdr*(y29-y23)/r;
            f24 = f24 - (ax+2*qx*cos(OmegaRF*t2))*OmegaRF^2/4*y21;
            f25 = f25 - (ay+2*qy*cos(OmegaRF*t2))*OmegaRF^2/4*y22;
            f26 = f26 - (az+2*qz*cos(OmegaRF*t2))*OmegaRF^2/4*y23; %ode end
            ...
        end
        ...
    end
    yf1 = y21; yf2 = y22; yf3 = y23; yf4 = y24; yf5 = y25; yf6 = y26;
    yf7 = y27; yf8 = y28; yf9 = y29; yf10 = y210; yf11 = y211; yf12 = y212;
    other = nsteps;
end