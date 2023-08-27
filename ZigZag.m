clear; close all; clc;
% ZigZag Manouevring in Deep water 

% Data
% Set initial time and velocity
t_initial = 0; 
u0 = 20/1.94384; 

% Set final time and number of steps
t_final = 180; 
N = 1000; 

% Calculate time step size
dt = t_final/N; 

% Define time range with a smaller time step for smoother plotting
dtt = 0.1; % Time step
Time = (0:dtt:t_final); % Time Range
L = length(Time); % Number of steps in time range
zz = zeros(L,2); % store value which used for comparison 

% Load a validation data for comparison
realdata=xlsread('ferry_zigzag.csv');
xreal=realdata(:,1);
yreal=realdata(:,2);


% Define rudder setup parameters
deltadot = -2.3; % Rate of change of rudder angle (degrees per second)
deltadot_rad = deg2rad(deltadot); % Convert rate of change to radians per second
delta_20 = deg2rad([20,20]); % Rudder angle for a 20 degree turn (in radians)


% Initialize simulation parameters
d = zeros(L,1); % Rudder angle at each time step
u = zeros(L,1); % Surge velocity at each time step
v = zeros(L,1); % Sway velocity at each time step
r = zeros(L,1); % Yaw rate at each time step
psi = zeros(L,1); % Yaw angle at each time step
ud = zeros(L,1); % Surge Accleration derivative at each time step
vd = zeros(L,1); % Sway Accleration derivative at each time step
rd = zeros(L,1); % Yaw Acceleration at each time step
xd = zeros(L,1); % xd position at each time step
yd = zeros(L,1); % yd position at each time step
x = zeros(L,1);
y = zeros(L,1); 


% Set up initial conditions before maneuvering
xi  = [u0,0,0,0,0,0,0];  % Initial state vector [u, v, r, x, y, psi, d]
u(1,1)     = xi(1);
v(1,1)     = xi(2);
r(1,1)     = xi(3);
x(1,1)     = xi(4);
y(1,1)     = xi(5);
psi(1,1)   = xi(6);
d(1,1) = xi(7);

zz = zeros(N, 2);

% Simulate maneuvering over time steps
for i = 1:L-1
    time = (i-1)*dt;
% Limit rudder angle to the maximum allowed value
    if d(i) >= delta_20(2)
       d(i) = delta_20(2);
    end 
    if d(i) <= -delta_20(2)
       d(i) = -delta_20(2);
    end
    % Execute zigzag maneuver after a certain time has passed
    if time > delta_20(2)
    if psi(i) >= -d(i) 
        % Turn rudder to the left
       d(i+1) = d(i) - dt*deltadot_rad;

    elseif psi(i) <= delta_20(2) 
        % Turn rudder to the right
        d(i+1) = d(i) + dt*deltadot_rad; 
    end 
    end

    % Compute acceleration and velocity components
    [accI]  = Motion(u(i),v(i),r(i),d(i));
    
    ud(i)    = accI(1);
    vd(i)    = accI(2);
    rd(i)    = accI(3);

    % Update velocity and position using the Euler scheme
    u(i+1)   = u(i) + dt*ud(i);
    v(i+1)   = v(i) + dt*vd(i);
    r(i+1)   = r(i) + dt*rd(i);
    psi(i+1) = psi(i) + dt*r(i);

    xd(i)  = (u(i)*cos(psi(i)) - v(i)*sin(psi(i)));
    x(i+1)   = x(i) + dt * xd(i);
   
    yd(i)  = (u(i)*sin(psi(i)) + v(i)*cos(psi(i)));
    y(i+1)   = y(i) + dt * yd(i);  
    % Store the current position in the zigzag matrix
    zz(i, :) = [x(i), y(i)];
   
end

%Convert unit 
d = rad2deg(d); 
r     = rad2deg(r);
rd    = rad2deg(rd);
psi   = rad2deg(psi);

% Data required for comparison 
save('udeep.mat', 'u', '-v7.3')
save('vdeep.mat', 'v', '-v7.3')
save('rdeep.mat', 'r', '-v7.3')
save('updeep.mat', 'ud', '-v7.3')
save('vpdeep.mat', 'vd', '-v7.3')
save('rpdeep.mat', 'rd', '-v7.3')


% All figure 
figure(1)
plot(-y,x,'k','LineWidth',1.2)
hold on
scatter(xreal,yreal, 'r', 'LineWidth', 1.2, 'SizeData',10)
xlabel('\psi [deg]'),ylabel('Time [s]')
title('20/20-ZigZag manoeuvre')
legend('Deep Water', 'Validation Data')
axis equal
pbaspect([1 4 1])
hold off

figure(2)
plot(Time, d, '--b', 'LineWidth',1)
hold on 
plot(Time,psi, 'r', 'LineWidth',1)
title('Zig-Zag Manoeuvring')
xlabel('Time [s]')
ylabel('\delta \psi [deg]')
legend('Rudder angle \psi', 'Course angle \delta')
grid on
set(gca, 'XAxisLocation', 'origin')
axis padded
hold off

figure(3)
plot(Time,ud,'k', 'LineWidth',1),xlabel('time [s]'),ylabel('u [m^2/s]'),title('Longitudinal Acceleration'),grid
figure(4)
plot(Time,vd,'b','LineWidth',1),xlabel('time [s]'),ylabel('v [m^2/s]'),title('Transverse Acceleration'),grid
figure(5)
plot(Time,rd,'r','LineWidth',1),xlabel('time [s]'),ylabel('r [rad²/s]'), title('Yaw Acceleration'),grid
figure(6)
plot(Time,d,'k','LineWidth',1),xlabel('time [s]'),ylabel('rudder angle [deg]'),title('Rudder Movement During Manoeuvring'),grid
figure(7)
plot(Time,u,'k','LineWidth',1),xlabel('time [s]'), ylabel('velocity [m/s]'),title('Longitudinal Velocity'),grid
figure(8)
plot(Time,v,'b','LineWidth',1),xlabel('time [s]'), ylabel('velocity [m/s]'), title('Transverse Velocity'),grid
figure(9)
plot(Time,r,'r','LineWidth',1),xlabel('time [s]'), ylabel('Turning rate [rad/s]'),title('Turning rate'),grid


% Function to calculate acceleration of a ship based on its current speed, heading rate and rudder angle

%accI: acceleration in ship reference frame (m/s^2)
function accI = Motion(u,v,r,d)
rudder_rate = (2.3*pi)/180;
U0 = 20/1.94384;
Lpp = 16*8.725;
U = sqrt(u^2 + v^2);
t_final = 180;
N = 1000;
dt = t_final/N;

u = u/U;
v = v/U;
r = r*Lpp/U;
du = u-U0/U;
m  = 6765*10^-6;
I = 319*10^-6;

% The following coefficients are hydrodynamic coefficients obtained from Wolf. These coefficients are used
% to calculate the forces acting on the ship based on its motion.
xg = (-116*10^-6)/m;
Xu    = -4336*10^-6;
Xuu   = -2355*10^-6;
Xuuu  = -2594*10^-6;
Xvv   = -3279*10^-6;
Xr    = -19*10^-6;
Xrr   = -571*10^-6;
Xdd   = -2879*10^-6;
Xdddd = 2185*10^-6;
Xvvu  = -2559*10^-6;
Xrru  = -734*10^-6;
Xddu  = 3425*10^-6;
Xvr   = 4627*10^-6;
Xvd   = 877*10^-6;
Xrd   = -351*10^-6;
X = Xu.*du + Xuu.*du.^2 + Xuuu.*du.^3 + Xvv.*v.^2 + Xr.*r + Xrr.*r.^2 + Xdd.*d.^2 + Xdddd.*d.^4 + Xvvu.*du.*v.^2 + Xrru.*du.*r.^2 + Xddu.*du.*d.^2 + Xvr.*v.*r + Xvd.*v.*d + Xrd.*r.*d;
Yu    = 57*10^-6;
Yv    = -12095*10^-6;
Yvvv  = -137302*10^-6;
Yvp   = -7396*10^-6;
Yr    = 1901*10^-6;
Yrrr  = -1361*10^-6;
Yrp   = -600*10^-6;
Yd    = 3587*10^-6;
Ydd   = 98*10^-6;
Yddddd= -6262*10^-6;
Yru   = -1297*10^-6;
Ydu   = -5096*10^-6;
Ydddu = 3192*10^-6;
Yvrr  = -44365*10^-6;
Yvvr  = -36490*10^-6;
Yvdd  = 2199*10^-6;
Yrdd  = -2752*10^-6;
Y = Yu.*du + Yv.*v + Yvvv.*v.^3 + Yr.*r + Yrrr.*r.^3 + Yd.*d + Ydd.*d.^2 + Yddddd.*d.^5 + Yru.*r.*du + Ydu.*du.*d + Ydddu.*du.*d.^3 + Yvrr.*v.*r.^2 + Yvvr.*r.*v.^2 + Yvdd.*v.*d.^2 + Yrdd.*r.*d.^2;
Nu    = -36*10^-6;
Nv    = -3919*10^-6;
Nvvv  = -33857*10^-6;
Nvp   = 426*10^-6;
Nvpvv = 10049*10^-6;
Nr    = -2579*10^-6;
Nrrr  = -2253*10^-6;
Nrp   = -231*10^-6;
Nd    = -1621*10^-6;
Ndd   = -73*10^-6;
Nddddd= 2886*10^-6;
Nvu   = -3666*10^-6;
Nrrru = -1322*10^-6;
Ndu   = 2259*10^-6;
Ndddu = -1382*10^-6;
Nvvr  = -60110*10^-6;
Nvdd  = 570*10^-6;
Nvvd  = -2950*10^-6;
Nrdd  = 237*10^-6;
Nrrd  = -329*10^-6;
N = Nu.*du + Nv.*v + Nvvv.*v.^3 + Nr.*r + Nrrr.*r.^3 + Nd.*d + Ndd.*d.^2 + Nddddd.*d.^5 + Nvu.*du.*v + Nrrru.*du.*r.^3 + Ndu.*d.*du + Ndddu.*du.*d.^3 + Nvvr.*r.*v.^2 + Nvdd.*v.*d.^2 + Nvvd.*d.*v.^2 + Nrdd.*r.*d.^2 + Nrrd.*d.*r.^2;
A = [m,0,0;0,m-Yvp,m*xg-Yrp;0,m*xg-Nvp-Nvpvv*v^2,I-Nrp];
b = [X+m*v*r+m*xg*r^2;Y-m*u*r;N-m*xg*u*r];
% Solve the system of equations to get acceleration in body axis
accI = A^-1*b;

accI(1)  =   accI(1)*(U^2)/Lpp;
accI(2)  =   accI(2)*(U^2)/Lpp;
accI(3)  =   accI(3)*(U^2)/(Lpp^2);
u=u*U;
v=v*U;
r=r*U/Lpp;
end
