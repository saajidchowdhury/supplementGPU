function g = g(C4,C8,mion,matom,ax,ay,az,qx,qy,qz,OmegaRF,t,y)
    g = zeros(12,1);
    g(1) = y(4); g(2) = y(5); g(3) = y(6);
    g(7) = y(10); g(8) = y(11); g(9) = y(12);
    r = sqrt((y(1)-y(7))^2+(y(2)-y(8))^2+(y(3)-y(9))^2);
    dVdr = 4*C4/r^5 - 8*C8/r^9;
    g(4) = -1/mion*dVdr*(y(1)-y(7))/r;
    g(5) = -1/mion*dVdr*(y(2)-y(8))/r;
    g(6) = -1/mion*dVdr*(y(3)-y(9))/r;
    g(10) = -1/matom*dVdr*(y(7)-y(1))/r;
    g(11) = -1/matom*dVdr*(y(8)-y(2))/r;
    g(12) = -1/matom*dVdr*(y(9)-y(3))/r;
    g(4) = g(4) - (ax+2*qx*cos(OmegaRF*t))*OmegaRF^2/4*y(1);
    g(5) = g(5) - (ay+2*qy*cos(OmegaRF*t))*OmegaRF^2/4*y(2);
    g(6) = g(6) - (az+2*qz*cos(OmegaRF*t))*OmegaRF^2/4*y(3);
end



