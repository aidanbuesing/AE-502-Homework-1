function chi = univAnomaly(dt,r0,vr0,mu,alpha)

% chi = universal variable chi needed in twoBodyOrbitProp

% dt = time in seconds
% ro = initial position
% v0 = initial velocity
% mu = gravitiational constant of the body we are assuming 2 body motion
% around
% alpha = 1/a where a is the semimajor axis


% Initial Guess
%chi = sqrt(mu)*abs(alpha)*dt;
chi =0;
% Initialize stopping criterion
ratio = 1; 
tol = 1e-8;

%Curtis method (Newton's method) for finding chi
while (ratio>tol)
    z = alpha*chi^2;

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

    f = (r0*vr0/sqrt(mu))*chi^2*C + (1-alpha*r0)*chi^3*S + r0*chi - sqrt(mu)*dt;
    f_prime = (r0*vr0/sqrt(mu))*chi*(1-z*S)+(1-alpha*r0)*chi^2*C +r0;
    ratio = abs(f/f_prime);
    chi = chi - f/f_prime;
    
end

