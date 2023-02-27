%% AE 502 Hwk 1
% Aidan Buesing
clear all
clc
% Earth r and v
%*1.496e8 r conversion to metric if i want
%*1731.46 v conversion to metric if i want
Re =  [-1.796136509111975*10e-1,9.667949206859814*10e-1,-3.668681017942158*10e-5];
Ve = [-1.720038360888334*10e-2,-3.211186197806460*10e-3,7.927736735960840*10e-7];

% 1I/Oumouamoua r and v
R1i =  [3.515868886595499*10e-2, -3.162046390773074, 4.493983111703389];
V1i =  [-2.317577766980901*10e-3,9.843360903693031*10e-3,-1.541856855538041*10e-2];
% 2I/Borisov r and v
R2i =  [7.249472033259724, 14.61063037906177, 14.24274452216359]*1.496e8;
V2i =  [-8.241709369476881*10e-3, -1.156219024581502*10e-2, -1.317135977481448*10e-2];

% Mu of sun
%mu = 1.327e11; % this is metric units...
mu = 0.0172; % this is in our current unit system
dt =1;

%% Problem 3

%This is there the delta_v's are stored, first index is days after
%departure, second is day of arrival in days after the "start date" Jan 1
delta_v_rendezvous = zeros(365,548);
delta_v_flyby = zeros(365,548);
 
% This will store the orbit r and v for each day, making it easy to plug
% into lambert. I had decided not to do this and just do a nested for loop
% essentially where these values would end up not stored. I decided to do
% it this way becuase I thought it made it easier to debug.
Earth_pos = zeros(1278,3);
Earth_vel = zeros(1278,3);
Oum_pos = zeros(760,3);
Oum_vel =zeros(760,3);


%Propogating Earth's Orbit
[Earth_pos(1,:),Earth_vel(1,:)]=twoBodyOrbitProp(Re,Ve,dt,mu);
for dt2dep = 2:1278 %365 different departure dates
    [Earth_pos(dt2dep,:),Earth_vel(dt2dep,:)]=twoBodyOrbitProp(Earth_pos(dt2dep-1,:),Earth_vel(dt2dep-1,:),dt,mu);
end

%And now Oum's Orbit
[Oum_pos(1,:),Oum_vel(1,:)]=twoBodyOrbitProp(R1i,V1i,dt,mu);
for dt2arr = 2:760
    [Oum_pos(dt2arr,:),Oum_vel(dt2arr,:)]=twoBodyOrbitProp(Oum_pos(dt2arr-1,:),Oum_vel(dt2arr-1,:),dt,mu);
end

%Solve lamberts problem using Earth's position on the departure day (i) and
%Oum's position on the arrival day (j) 213 days is the time from Jan 1 to
%Aug 1, thus ignoring days 1-212 which we are not considering arriving on
for i = 1:365 

    for j = 213:760
        [v1,v2]=lambertSolver(Earth_pos(i,:),Oum_pos(j,:),(j-i)*dt,mu);
        delta_v_rendezvous(i,j)=(norm(v1)+norm(v2));
        delta_v_flyby(i,j)=norm(v1);
        
    end
end

%Conversion to metric for graphs
delta_v_rendezvous = delta_v_rendezvous*1731.46;
delta_v_flyby = delta_v_flyby*1731.46;

%Making the plots 
dep_dat = 1:365;
arr_dat = 213:760;
figure (1)
[c,h]=contour(arr_dat,dep_dat,delta_v_rendezvous(dep_dat,arr_dat));
xlabel('Day of Arrival in Days After August 1st 2017')
ylabel('Day of Departure in Days After January 1st 2017')
title('Delta_V for Rendez-Vous with Oumouamoua')

figure (2)
[c2,h2]=contour(arr_dat,dep_dat,delta_v_flyby(dep_dat,arr_dat));
xlabel('Day of Arrival in Days After August 1st 2017')
ylabel('Day of Departure in Days After January 1st 2017')
title('Delta_V for Fly-By of Oumouamoua')

%Will talk about these in the write up... they're clearly troubled

%% Problem 4

%This problem is the same as the last but we change to Borisov from Oum and
%we also change our dates
delta_v_rendezvous2 = zeros(1278,1855);
delta_v_flyby2 = zeros(1278,1855);
 
Bor_pos = zeros(1855,3);
Bor_vel =zeros(1855,3);

%Getting Borisov's position at all days
[Bor_pos(1,:),Bor_vel(1,:)]=twoBodyOrbitProp(R2i,V2i,dt,mu);
for dt2arr = 2:1855
    [Bor_pos(dt2arr,:),Bor_vel(dt2arr,:)]=twoBodyOrbitProp(Bor_pos(dt2arr-1,:),Bor_vel(dt2arr-1,:),dt,mu);
end

%Solving Lambert's problem again but with Borisov and the new dates
for i = 1:1278

    for j = 910:1855
        [v1,v2]=lambertSolver(Earth_pos(i,:),Bor_pos(j,:),(j-i)*dt,mu);
        delta_v_rendezvous2(i,j)=(norm(v1)+norm(v2));
        delta_v_flyby2(i,j)=norm(v1);
       
    end
end

%Conversion to metric for graphs
delta_v_rendezvous2 = delta_v_rendezvous2*1731.46;
delta_v_flyby2 = delta_v_flyby2*1731.46;

%More plots
dep_dat = 1:1278;
arr_dat = 910:1855;
figure (3)
[c3,h3]=contour(arr_dat,dep_dat,delta_v_rendezvous2(dep_dat,arr_dat));
xlabel('Day of Arrival in Days After January 1st 2017')
ylabel('Day of Departure in Days After January 1st 2017')
title('Delta_V for Rendez-Vous with Borisov')


figure (4)
[c4,h4]=contour(arr_dat,dep_dat,delta_v_flyby2(dep_dat,arr_dat));
xlabel('Day of Arrival in Days After January 1st 2017')
ylabel('Day of Departure in Days After January 1st 2017')
title('Delta_V for Fly-By of Borisov')


% Same issue with them, talking about that in the write up

%% Problem 5

%I would just make an array to dump these all in if this was a more complex
%problem, but it's not, so just variables.
[Oum_h,Oum_i,Oum_rightAsc,Oum_e,Oum_ArgPer,Oum_theta] = orbitalElementsFromRV(R1i,V1i,mu);
[Bor_h,Bor_i,Bor_rightAsc,Bor_e,Bor_ArgPer,Bor_theta] = orbitalElementsFromRV(R2i,V2i,mu);

Oum_e
Bor_e

%Both eccentricities are >1 so they are both hyperbollic and thus
%interstellar. Although they are really *really* large... too large? Maybe,
%but all of the other values make sense to me, and well I did find these to
%be interstellar so I feel as if my doubts of the eccentricity being so
%high are perhaps unfounded and it's just fine.

%Also, I answer Problem 6 in the write up.


