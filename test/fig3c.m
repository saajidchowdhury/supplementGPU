function h = h(C4,C8,mion,matom,ax,ay,az,qx,qy,qz,OmegaRF,t,y)
    h = zeros(12,1);
    h(1) = y(4); h(2) = y(5); h(3) = y(6);
    h(7) = y(10); h(8) = y(11); h(9) = y(12);
    r = sqrt((y(1)-y(7))^2+(y(2)-y(8))^2+(y(3)-y(9))^2);
    dVdr = 4*C4/r^5 - 8*C8/r^9;
    mi = -1/mion*dVdr;
    ma = -1/matom*dVdr;
    c = cos(OmegaRF*t);
    o = OmegaRF^2/4;
    h(4) = mi*(y(1)-y(7))/r;
    h(5) = mi*(y(2)-y(8))/r;
    h(6) = mi*(y(3)-y(9))/r;
    h(10) = ma*(y(7)-y(1))/r;
    h(11) = ma*(y(8)-y(2))/r;
    h(12) = ma*(y(9)-y(3))/r;
    h(4) = h(4) - (ax+2*qx*c)*o*y(1);
    h(5) = h(5) - (ay+2*qy*c)*o*y(2);
    h(6) = h(6) - (az+2*qz*c)*o*y(3);
end