function [v1,v2] = lambertSolver(r1,r2,dt,mu)

% v1 = delta v1 of solved lambert problem
% v2 = delta v2 of solved lambert problem
% r1 = position vector of initial orbit
% r2 = position vector of final orbit
% dt = time from r1 to r2 in seconds
% mu = gravitiational constant of the body we are assuming 2 body motion
% around

r1m =sqrt(dot(r1,r1));
r2m =sqrt(dot(r2,r2));

grade = 1; %1=prograde 0=retrograde

if (grade==1)
    if (cross(r1,r2)>=0)
        dtheta = acos(dot(r1,r2)/(r1m*r2m));
    else
        dtheta = 2*pi - acos(dot(r1,r2)/(r1m*r2m));
    end
end
if (grade ==0)
    if (cross(r1,r2)>=0)
         dtheta = 2*pi - acos(dot(r1,r2)/(r1m*r2m));
    else
        dtheta = acos(dot(r1,r2)/(r1m*r2m));
    end
end

A = sin(dtheta)*sqrt((r1m*r2m)/(1-cos(dtheta)));

% Find z iteratively

ratio = 1;
tol = 10e-8;
z=0;
maxIter =1000;
count =0;
while (abs(ratio)>tol && count <= maxIter)
    if (z==0)
        S = 1/6;
        C = 1/2;
        y = r1m +r2m + A*((z*S-1)/sqrt(C));
        F = (y/C)^1.5*S + A*sqrt(y) - sqrt(mu)*dt;
        Fprime = sqrt(2)/40*y^1.5 + A/8*(sqrt(y) + A*sqrt(1/2/y));
        ratio = F/Fprime;
        z = z - ratio;
    end
    if (z > 0)
         S = (sqrt(z)-sin(sqrt(z)))/(sqrt(z)^3);
         C = (1-cos(sqrt(z)))/z;
         y = r1m +r2m + A*((z*S-1)/sqrt(C));
         F = (y/C)^1.5*S + A*sqrt(y) - sqrt(mu)*dt;
         Fprime = (y/C)^1.5*(1/2/z*(C - 3*S/2/C) ...
               + 3*S^2/4/C) + A/8*(3*S/C*sqrt(y) ...
               + A*sqrt(C/y));
         ratio = F/Fprime;
         z = z - ratio;
    end
    if (z < 0)
        S = (sinh(sqrt(-z))-sqrt(-z))/(sqrt(z)^3);
        C = (cosh(sqrt(-z))-1)/(-z);
        y = r1m +r2m + A*((z*S-1)/sqrt(C));
        F = (y/C)^1.5*S + A*sqrt(y) - sqrt(mu)*dt;
        Fprime = (y/C)^1.5*(1/2/z*(C - 3*S/2/C) ...
               + 3*S^2/4/C) + A/8*(3*S/C*sqrt(y) ...
               + A*sqrt(C/y));
        ratio = F/Fprime;
        z = z - ratio;
    end
    count = count +1;
end

y = r1m +r2m + A*((z*S-1)/sqrt(C));

f = 1 - y/r1m;
g = A*sqrt(y/mu);
%fdot = (sqrt(mu)/(r1m*r2m))*sqrt(y/C)*(z*S-1);
gdot = 1 - (y/r2m);

v1 = (1/g)*(r2-(f*r1));
v2 = (1/g)*(gdot*r2 - r1);


