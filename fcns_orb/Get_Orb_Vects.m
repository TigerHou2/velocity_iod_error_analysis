function [r,v] = Get_Orb_Vects(params,mu)
%GET_ORB_VECTS Finds position and velocity vectors from orbital parameters.

a = params(1);
e = params(2);
i = params(3);
omg = params(4);
w = params(5);
f = params(6);

rr = a*(1-e^2)/(1+e*cos(f));
h = sqrt(mu*a*(1-e^2));

r = rr * [ cos(omg)*cos(w+f) - sin(omg)*sin(w+f)*cos(i) ;...
           sin(omg)*cos(w+f) + cos(omg)*sin(w+f)*cos(i) ;...
           sin(w+f)*sin(i) ];
v = -mu/h * ...
    [ cos(omg)*(sin(w+f)+e*sin(w)) + sin(omg)*(cos(w+f)+e*cos(w))*cos(i) ;...
      sin(omg)*(sin(w+f)+e*sin(w)) - cos(omg)*(cos(w+f)+e*cos(w))*cos(i) ;...
    -(cos(w+f)+e*cos(w))*sin(i) ];

end

