function [h,i,rightAsc,e,argPer,theta] = orbitalElementsFromRV(r,v,mu)

% r = 1x3 position vector
% v = 1x3 velocity vector 
% mu = gravitational constant of the body we are assuming 2 body motion
% around


% This is step by step following the curtis algorithm in chapter 4

rm = sqrt(dot(r,r));
vm = sqrt(dot(v,v));

vr = dot(r,v)/rm;

hvec = cross(r,v);

h = dot(hvec,hvec);
i = acos(hvec(3)/h);

K = [0 0 1];
N = cross(K,hvec);
n = sqrt(dot(N,N));

if N(2)>=0
    rightAsc = acos(N(1)/n);
else 
    rightAsc = 2*pi - acos(N(1)/n);
end

E = (1/mu) * ((vm^2-(2*mu)*r - rm*vr*v)) ;
e = sqrt(dot(E,E));

if E(3)>=0
    argPer = acos((dot(N,E))/(n*e));
else
    argPer = 2*pi - acos((dot(N,E))/(n*e));
end

theta = acos(dot(E/e,r/rm));



