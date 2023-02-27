function [r,v] = twoBodyOrbitProp(r0,v0,dt,mu)

% r = position vector after time dt
% v = velocity vector after time dt
% r0 = position vector of initial orbit
% v0 = velocity vector of initial orbit
% dt = time from r1 to r2 in seconds
% mu = gravitiational constant of the body we are assuming 2 body motion
% around


r0m = sqrt(dot(r0,r0));
v0m = sqrt(dot(v0,v0));

vr0 = (dot(r0,v0))/r0m;
alpha = (2/r0m)-(v0m^2/mu);

%Calls univAnomaly to find chi
chi = univAnomaly(dt,r0,vr0,mu,alpha);

%defn of z
z = alpha*(chi^2);

%find values of stumpff functions S and C 
if (z>0)
    S = (sqrt(z)-sin(sqrt(z)))/(sqrt(z)^3);
    C = (1-cos(sqrt(z)))/z;
elseif (z<0)
    S = (sinh(sqrt(-z))-sqrt(-z))/(sqrt(z)^3);
    C = (cosh(sqrt(-z))-1)/(-z);
else
    S = 1/6;
    C = 1/2;
end

% use f and g functions to solve r and v of propogated orbit
f = 1-((chi^2/r0m)*C);
g = dt - ((1/sqrt(mu))*chi^3*S);
r = f*r0 + g*v0;
rm = sqrt(dot(r,r));
fdot = (sqrt(mu)/(rm*r0m))*((alpha*(chi^3)*S)-chi);
gdot = 1 - (chi^2/rm)*C;
v= fdot*r0 + gdot*v0;
