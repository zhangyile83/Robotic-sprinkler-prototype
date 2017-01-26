%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Wind shifting Technique %%%%%%%%%%%%%%%%%%%%%%%%
% __ wind effect correction for 360 degree case
% __ Multi droplet diameter
% __ Multiple spraying circle
% __ corrections compute in advance
% __ Smooth version

clear;
tic
figure(1)
pause

%______________Configuration of nozzle_________________
alpha = pi/6;  % The angle of elevation of the nozzle 
d_range = 0.005:0.001:0.010; % The droplet diameter range
MeanD = 0.01;  % Droplet mean diameter(should be the value defined in the range)
CC = 100; % CC here is a measure of water flowing volume per unit time, while the value need be verified by experiment data

%______________Configuration of water distribution plotting_________________
RSL1 = 5; % Resolution in the water distribution figure: the unit is meter
Len = 300; % Water distribution plotting range(length): the unit is meter, note that the length and width are too big for our simulation case, therefore to show the properly show the results, tailoring is applied
Wid = 300; % Water distribution plotting range(width)
yic = 20;  % Increment at y direction in the water distribution graph to compensate the Matlab plotting discrepancy
GridsR1w = zeros(Wid/RSL1, Len/RSL1); % The matrix that collecting the water
Trejectory_width = 1; % This merely a plotting setting, the width here signify the width of the trejectory curve in the trejectory plotting

%______________Configuration of spraying_________________
gama_temp = 0:pi/32:pi;  % gama_temp is used to generated the contour of the lawn
spraying_angle = 2*pi:pi/32:4*pi-pi/32; % spraying_angle signify the spraying angle in horizontal plane, note we start from 2pi to 4pi instead of 0 to 2pi because this avoids the negative angle in the correction

% ___________________________________________________
%______________Spraying circle begin_________________
% ___________________________________________________
for spraying_C = [30, 28, 22, 16, 11]; % define the target distance in each spraying circle (Assume the lawn is rectangle)
    
    p1 = []; % p1 stores the velocity correction parameter
    p2 = []; % p2 stores the angle correction parameter
    TT = []; % TT stores the droplet flying time in the air for different spraying velocitys
     
    Trejectory_width = Trejectory_width + 0.2; % when we finish one spraying circle, we use bolder line to signify the trajectory in next spraying circle    
    Target_contour = zeros(length(gama_temp), 1);  % contour target distance

%______________Configuration of spraying target distance for each spraying circle_________________
for i = 1:find(gama_temp == pi/4)
    Target_contour(i) = spraying_C/cos(gama_temp(i));    
end
for i = find(gama_temp == pi/4) + 1:find(gama_temp == pi/2)
    Target_contour(i) = spraying_C/sin(gama_temp(i));    
end
for i = find(gama_temp == pi/2) + 1:find(gama_temp == 3*pi/4)
    Target_contour(i) = spraying_C/cos(gama_temp(i) - pi/2);    
end   
for i = find(gama_temp == 3*pi/4) + 1:find(gama_temp == pi)
    Target_contour(i) = spraying_C/cos(pi - gama_temp(i));    
end

ctemp = Target_contour(2:length(Target_contour)-1);
Target_contour_distance = [Target_contour; ctemp]; % The target distance for different spraying angles

% ________ Find the required flow speed for different spraying angles according
% to target contour distance _________
[v0S] = Characteristic(alpha, Target_contour_distance); % More details refer to Characteristic.m

%_______ Spraying loop at different angle in each spraying circle ______
for jj = 1:length(spraying_angle)
gama = spraying_angle(jj); % gama is the spraying angle
v0 = v0S(jj); % v0 is the required velocity to reach target contour distance

d = MeanD; % define the mean diameter, store it in a .mat file to pass the parameter to other functions

% __________________________________________________________________________
%______________Computing droplet dropping time without wind_________________
% __________________________________________________________________________
% _________Droplet speed in x, y and z coordinate, v0 is the initial velocity of droplet____________
vx0 = v0*cos(alpha)*cos(gama);
vy0 = v0*cos(alpha)*sin(gama);
vz0 = v0*sin(alpha);

% ________ apply ODE solver to solve the flow velocity without wind and
% corresponding time ________
[T_temp,Y_temp] = ode45(@(t, y) droplet(t, y, d),[0 10],[vx0 vy0 vz0]);

vz_temp = Y_temp(:,3);
Dz_temp = zeros(length(vz_temp),1);

for i = 2:length(T_temp)
Dz_temp(i) = trapz(T_temp(1:i),vz_temp(1:i));
end

PolyN = 4;
a = polyfit(T_temp, Dz_temp, PolyN);

t_temp = roots(a);

t = t_temp(3);  % t is the dropping time without wind

% __________________________________________________________________________
%_____________Simulation of droplet trajectory without wind_________________
% __________________________________________________________________________

vx0 = v0*cos(alpha)*cos(gama);
vy0 = v0*cos(alpha)*sin(gama);
vz0 = v0*sin(alpha);

[T,Y] = ode45(@(t,y) droplet(t,y,d),[0 t],[vx0 vy0 vz0]);

vx = Y(:,1);
vy = Y(:,2);
vz = Y(:,3);

Dx = [];
Dy = [];
Dz = [];

for i = 2:length(T)
Dx(i) = trapz(T(1:i),vx(1:i));
Dy(i) = trapz(T(1:i),vy(1:i));
Dz(i) = trapz(T(1:i),vz(1:i));
end

% ____________________________________________________________
%_____________Wind effect shifting simulation_________________
% ____________________________________________________________

%_____________Find shifting parameters_________________
% [p1(jj), p2(jj), TT(jj)] = Anti_wind_parameter_v3(v0, Dx, Dy, gama) % Anti_wind_parameter_v3.m output the shifting parameters.
end

%_____________Spraying simulation for different angles with Wind shifting parameters_________________
for jj = 1:length(spraying_angle)
gama = spraying_angle(jj);
v0 = v0S(jj);

for ii = 1:length(d_range)
    d = d_range(ii);
save('Diameter.mat','d'); % Save the different droplet diameter and pass them into the function in @droplet

vx0 = v0*cos(alpha)*cos(gama);
vy0 = v0*cos(alpha)*sin(gama);
vz0 = v0*sin(alpha);

% ____________Dropping time for droplet with different diameters_____________

[T_temp,Y_temp] = ode45(@(t,y) droplet(t,y,d),[0 10],[vx0 vy0 vz0]);

vz_temp = Y_temp(:,3);

Dz_temp = zeros(length(vz_temp),1);

for i = 2:length(T_temp)
Dz_temp(i) = trapz(T_temp(1:i),vz_temp(1:i));
end

PolyN = 4;
a = polyfit(T_temp, Dz_temp, PolyN);

t_temp = roots(a);

t3 = t_temp(3);  % t is the dropping time for droplet with certain diameter


%____Trajectory simulation under wind condition after shifting_____
[T1,Y1] = ode45(@(t,y) droplet_wind_Modi(t, y, d),[0 t3],[vx0 vy0 vz0]);

vx1 = Y1(:,1);
vy1 = Y1(:,2);
vz1 = Y1(:,3);

Dx1 = [];
Dy1 = [];
Dz1 = [];

for i = 2:length(T1)
Dx1(i) = trapz(T1(1:i),vx1(1:i));
Dy1(i) = trapz(T1(1:i),vy1(1:i));
Dz1(i) = trapz(T1(1:i),vz1(1:i));
end

%____For the plotting purpose, we only show the trajectory of the droplet with mean diameter_____
if ii == find(d_range == MeanD)

grid on

pause(0.001)
subplot(2,2,1:2);
title({['11 MPH wind without shifting']}, 'FontSize', 15);

hold
%____Plotting trajectory_____
plot3(Dx1,Dy1,Dz1,'LineWidth',Trejectory_width);

%______ plotting the black box which shows the correct range of water _____
xbr = 30*ones(61,1);  % right contour
ybr = [-30:1:30];
zbr = zeros(61,1);

xbl = -30*ones(61,1);  % left contour
ybl = [-30:1:30];
zbl = zeros(61,1);

ybb = -30*ones(61,1);  % right contour
xbb = [-30:1:30];
zbb = zeros(61,1);

ybt = 30*ones(61,1);  % right contour
xbt = [-30:1:30];
zbt = zeros(61,1);

plot3(xbr, ybr, zbr, 'k', 'linewidth', 2);
plot3(xbl, ybl, zbl, 'k', 'linewidth', 2);
plot3(xbb, ybb, zbb, 'k', 'linewidth', 2);
plot3(xbt, ybt, zbt, 'k', 'linewidth', 2);

axis([-40 40 -40 40]);
xlabel('x (m)', 'FontSize', 15);
ylabel('y (m)', 'FontSize', 15);
zlabel('z (m)', 'FontSize', 15);
set(gca,'XGrid','on')
set(gca,'YGrid','on')

hold

end

Dd1(ii) = Dx1(end);
Dxtemp1w = round(Dx1(end)/RSL1);
Dytemp1w = round(Dy1(end)/RSL1);



%______ The following setting is determined by the nozzle configurations _____
if spraying_C<=20
    CC = 60;
end

if spraying_C<= 15
    CC = 30;
end

%_______ Distribution of the number of droplet with different diameters___
Vd_temp = (pi/6)*pi*d^3*CC*distribution(d,MeanD,0.0015); %____ Normal distribution assumption is used ____
Vd_temp_1 = Vd_temp*4;

%_______ Accumulation of water at different bucket ________
GridsR1w(Dxtemp1w + Wid/2/RSL1, Dytemp1w + yic) = GridsR1w(Dxtemp1w + Wid/2/RSL1, Dytemp1w + yic) + Vd_temp_1;

end

subplot(2,2,3:4)

myplot360B(GridsR1w');
end

end