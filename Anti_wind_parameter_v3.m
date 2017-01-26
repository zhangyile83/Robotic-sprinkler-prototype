function [p1, p2, t1] = Anti_wind_parameter_v3(v0, Dx, Dy, gama, d)

%__ Multi-time corrector for both velocity and angle
%__ NON-Debug Mode

threshold = 1;
scale_p1 = 0.05;
scale_p2 = 0.002;

%_ Multi-time corrector for both velocity and angle
%__ Debug Mode


% _____________ Input parameters ________________
% v0 is the flow velocity without wind
% Dx is the displacement at x-coordinate
% Dy is the displacement at y-coordinate
% w is the wind speed
% gamaw is the horizontal angle of the wind
% MeanD is the mean diameter
% t is the dropping time of case without correction

% _____________ Output parameters ________________
% alpha1 is the velocity correction parameter
% alpha2 is the angle correction parameter
% t1 is the dropping time of the case with optimal correction

% ____________find the dropping time for correction case___________________
% find the time that the droplet touch the ground after correction
alpha = pi/6;
DxO = Dx(end);
DyO = Dy(end);
p1_optimal = 1;
p2_optimal = 1;
ErrorMin = 30;

Condition11 = 1;
% p1 = [];
% p2 = [];
% p1_index = 1;
% P2_index = 1;
NUM = 0;
%_________________ Optimization of p1 __________________________

%___________________P1 - case 1 Initialization ____________________

while(Condition11)
    Optimal_find1 = 1;
    Optimal_find2 = 1;

    NUM = NUM + 1;

n1 = 1;
n2 = 1;

v0M1 = (1+n1*scale_p1)*v0*p1_optimal;
gamaM = gama*p2_optimal;

vx0M1 = v0M1*cos(alpha)*cos(gamaM);
vy0M1 = v0M1*cos(alpha)*sin(gamaM);
vz0M1 = v0M1*sin(alpha);
[T_tempM1,Y_tempM1] = ode45(@(t, y) droplet(t, y, d),[0 10],[vx0M1 vy0M1 vz0M1]);
vz_tempM1 = Y_tempM1(:,3);
Dz_tempM1 = zeros(length(vz_tempM1),1);

for i = 2:length(T_tempM1)
Dz_tempM1(i) = trapz(T_tempM1(1:i),vz_tempM1(1:i));
end


% size(Dz_tempM1)
% size(T_tempM1)

PolyN = 4;
a = polyfit(T_tempM1, Dz_tempM1, PolyN);

t_tempM1 = roots(a);

t11 = t_tempM1(3);

[T11,Y11] = ode45(@(t, y) droplet_wind_Modi(t, y, d),[0 t11],[vx0M1 vy0M1 vz0M1]);

vx11 = Y11(:,1);
vy11 = Y11(:,2);
vz11 = Y11(:,3);

Dx11 = [];
Dy11 = [];
Dz11 = [];

for i = 2:length(T11)
Dx11(i) = trapz(T11(1:i),vx11(1:i));
Dy11(i) = trapz(T11(1:i),vy11(1:i));
Dz11(i) = trapz(T11(1:i),vz11(1:i));
end
DxC1 = Dx11(end);
DyC1 = Dy11(end);
ErrorD1 = sqrt((DxO - DxC1)^2 + (DyO - DyC1)^2);


%___________________P1 - case 2 Initialization ____________________
v0M2 = (1-n1*scale_p1)*v0*p1_optimal;
gamaM = gama*p2_optimal;

vx0M2 = v0M2*cos(alpha)*cos(gamaM);
vy0M2 = v0M2*cos(alpha)*sin(gamaM);
vz0M2 = v0M2*sin(alpha);
[T_tempM2,Y_tempM2] = ode45(@(t, y) droplet(t,y,d),[0 10],[vx0M2 vy0M2 vz0M2]);

vz_tempM2 = Y_tempM2(:,3);
Dz_tempM2 = zeros(length(vz_tempM2),1);


for i = 2:length(T_tempM2)
    
Dz_tempM2(i) = trapz(T_tempM2(1:i),vz_tempM2(1:i));
end


PolyN = 4;
a = polyfit(T_tempM2, Dz_tempM2, PolyN);

t_tempM2 = roots(a);

t12 = t_tempM2(3);

[T12,Y12] = ode45(@(t, y) droplet_wind_Modi(t,y,d),[0 t12],[vx0M2 vy0M2 vz0M2]);

vx12 = Y12(:,1);
vy12 = Y12(:,2);
vz12 = Y12(:,3);

Dx12 = [];
Dy12 = [];
Dz12 = [];

for i = 2:length(T12)
Dx12(i) = trapz(T12(1:i),vx12(1:i));
Dy12(i) = trapz(T12(1:i),vy12(1:i));
Dz12(i) = trapz(T12(1:i),vz12(1:i));
end

DxC2 = Dx12(end);
DyC2 = Dy12(end);
ErrorD2 = sqrt((DxO - DxC2)^2 + (DyO - DyC2)^2);

if ErrorD1 < ErrorD2 && ErrorD1 < ErrorMin
    aaaaa = 1;
    ErrorMin = ErrorD1;
else if ErrorD2 < ErrorD1 && ErrorD2 < ErrorMin
    aaaaa = 2;
    ErrorMin = ErrorD2;
    else
        Optimal_find1 = 0;
    end
end

ErrorD = ErrorMin;
Condition = 1;


%___________________ iteration to find optimal p1 ________________________
while Condition && Optimal_find1
        
    n1 = n1 + 1;
    ErrorMin = ErrorD;

if aaaaa == 1
v0M = (1 + n1*scale_p1)*v0*p1_optimal;
gamaM = gama*p2_optimal;

vx0M = v0M*cos(alpha)*cos(gamaM);
vy0M = v0M*cos(alpha)*sin(gamaM);
vz0M = v0M*sin(alpha);
[T_tempM,Y_tempM] = ode45(@(t, y) droplet(t, y, d),[0 10],[vx0M vy0M vz0M]);
vz_tempM = Y_tempM(:,3);
Dz_tempM = zeros(length(vz_tempM),1);

for i = 2:length(T_tempM)
Dz_tempM(i) = trapz(T_tempM(1:i),vz_tempM(1:i));
end

PolyN = 4;
a = polyfit(T_tempM, Dz_tempM, PolyN);

t_tempM = roots(a);

t11 = t_tempM(3);

[T1,Y1] = ode45(@(t, y) droplet_wind_Modi(t,y,d),[0 t11],[vx0M vy0M vz0M]);

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
DxC = Dx1(end);
DyC = Dy1(end);
ErrorD = sqrt((DxO - DxC)^2 + (DyO - DyC)^2);
end


%___________________ minus case _________________________
if aaaaa == 2
v0M = (1-n1*scale_p1)*v0*p1_optimal;
gamaM = gama*p2_optimal;

vx0M = v0M*cos(alpha)*cos(gamaM);
vy0M = v0M*cos(alpha)*sin(gamaM);
vz0M = v0M*sin(alpha);
[T_tempM,Y_tempM] = ode45(@(t,y) droplet(t,y,d),[0 10],[vx0M vy0M vz0M]);
vz_tempM = Y_tempM(:,3);
Dz_tempM = zeros(length(vz_tempM),1);

for i = 2:length(T_tempM)
Dz_tempM(i) = trapz(T_tempM(1:i),vz_tempM(1:i));
end

PolyN = 4;
a = polyfit(T_tempM, Dz_tempM, PolyN);

t_tempM = roots(a);

t12 = t_tempM(3);

[T1,Y1] = ode45(@(t,y) droplet_wind_Modi(t, y, d),[0 t12],[vx0M vy0M vz0M]);

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
DxC = Dx1(end);
DyC = Dy1(end);
ErrorD = sqrt((DxO - DxC)^2 + (DyO - DyC)^2);


end

Condition = (ErrorD < ErrorMin);
end


if aaaaa == 1
p1 = (1 + (n1-1)*scale_p1)*p1_optimal; t1 = t11;
end
if aaaaa == 2
p1 = (1 - (n1-1)*scale_p1)*p1_optimal; t1 = t12;
end

p1_optimal = p1;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%_________________ Optimization of p2 __________________________

%_______________ P2 - case 1 Initialization ____________________
v0M1 = p1*v0;
gamaMC = (1 + n2*scale_p2)*gama*p2_optimal;

vx0M1 = v0M1*cos(alpha)*cos(gamaMC);
vy0M1 = v0M1*cos(alpha)*sin(gamaMC);
vz0M1 = v0M1*sin(alpha);
[T_tempM1,Y_tempM1] = ode45(@(t, y) droplet(t, y, d),[0 10],[vx0M1 vy0M1 vz0M1]);
vz_tempM1 = Y_tempM1(:,3);
Dz_tempM1 = zeros(length(vz_tempM1),1);

for i = 2:length(T_tempM1)
Dz_tempM1(i) = trapz(T_tempM1(1:i),vz_tempM1(1:i));
end

PolyN = 4;
a = polyfit(T_tempM1, Dz_tempM1, PolyN);

t_tempM1 = roots(a);

t11 = t_tempM1(3);

[T11,Y11] = ode45(@(t, y) droplet_wind_Modi(t, y, d),[0 t11],[vx0M1 vy0M1 vz0M1]);

vx11 = Y11(:,1);
vy11 = Y11(:,2);
vz11 = Y11(:,3);

Dx11 = [];
Dy11 = [];
Dz11 = [];

for i = 2:length(T11)
Dx11(i) = trapz(T11(1:i),vx11(1:i));
Dy11(i) = trapz(T11(1:i),vy11(1:i));
Dz11(i) = trapz(T11(1:i),vz11(1:i));
end
DxC1 = Dx11(end);
DyC1 = Dy11(end);
ErrorD1 = sqrt((DxO - DxC1)^2 + (DyO - DyC1)^2);



%___________________P2 - case 2 Initialization ____________________
v0M2 = p1*v0;
gamaMC = (1 - n2*scale_p2)*gama*p2_optimal;

vx0M2 = v0M2*cos(alpha)*cos(gamaMC);
vy0M2 = v0M2*cos(alpha)*sin(gamaMC);
vz0M2 = v0M2*sin(alpha);
[T_tempM2,Y_tempM2] = ode45(@(t, y) droplet(t,y,d),[0 10],[vx0M2 vy0M2 vz0M2]);
vz_tempM2 = Y_tempM2(:,3);
Dz_tempM2 = zeros(length(vz_tempM2),1);

for i = 2:length(T_tempM1)
Dz_tempM2(i) = trapz(T_tempM2(1:i),vz_tempM2(1:i));
end

PolyN = 4;
a = polyfit(T_tempM2, Dz_tempM2, PolyN);

t_tempM2 = roots(a);

t12 = t_tempM2(3);

[T12,Y12] = ode45(@(t,y) droplet_wind_Modi(t, y, d),[0 t12],[vx0M2 vy0M2 vz0M2]);

vx12 = Y12(:,1);
vy12 = Y12(:,2);
vz12 = Y12(:,3);

Dx12 = [];
Dy12 = [];
Dz12 = [];

for i = 2:length(T12)
Dx12(i) = trapz(T12(1:i),vx12(1:i));
Dy12(i) = trapz(T12(1:i),vy12(1:i));
Dz12(i) = trapz(T12(1:i),vz12(1:i));
end

DxC2 = Dx12(end);
DyC2 = Dy12(end);
ErrorD2 = sqrt((DxO - DxC2)^2 + (DyO - DyC2)^2);


if ErrorD1 < ErrorD2 && ErrorD1 < ErrorMin
    bbbbb = 1;
    ErrorMin = ErrorD1;
else if ErrorD2 < ErrorD1 && ErrorD2 < ErrorMin
    bbbbb = 2;
    ErrorMin = ErrorD2;
    else
        Optimal_find2 = 0;
    end
end


if Optimal_find2 == 0 && NUM == 1
    p2 = 1;
    break
end

ErrorD = ErrorMin;
Condition = 1;


%___________________ iteration to find p1 and p2________________________
while Condition && Optimal_find2
        
    n2 = n2 + 1;
    ErrorMin = ErrorD;

if bbbbb == 1

    v0M = p1*v0;
gamaMC = (1 + n2*scale_p2)*gama*p2_optimal;
vx0M = v0M*cos(alpha)*cos(gamaMC);
vy0M = v0M*cos(alpha)*sin(gamaMC);
vz0M = v0M*sin(alpha);
[T_tempM,Y_tempM] = ode45(@(t, y) droplet(t, y, d),[0 10],[vx0M vy0M vz0M]);
vz_tempM = Y_tempM(:,3);
Dz_tempM = zeros(length(vz_tempM),1);

for i = 2:length(T_tempM)
Dz_tempM(i) = trapz(T_tempM(1:i),vz_tempM(1:i));
end

PolyN = 4;
a = polyfit(T_tempM, Dz_tempM, PolyN);

t_tempM = roots(a);

t11 = t_tempM(3);

[T1,Y1] = ode45(@(t, y) droplet_wind_Modi(t, y, d),[0 t11],[vx0M vy0M vz0M]);

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
DxC = Dx1(end);
DyC = Dy1(end);
ErrorD = sqrt((DxO - DxC)^2 + (DyO - DyC)^2);

end


%___________________ minus case _________________________
if bbbbb == 2
v0M = p1*v0;
gamaMC = (1 - n2*scale_p2)*gama*p2_optimal;
vx0M = v0M*cos(alpha)*cos(gamaMC);
vy0M = v0M*cos(alpha)*sin(gamaMC);
vz0M = v0M*sin(alpha);
[T_tempM,Y_tempM] = ode45(@(t, y) droplet(t, y, d),[0 10],[vx0M vy0M vz0M]);
vz_tempM = Y_tempM(:,3);
Dz_tempM = zeros(length(vz_tempM),1);

for i = 2:length(T_tempM)
Dz_tempM(i) = trapz(T_tempM(1:i),vz_tempM(1:i));
end

PolyN = 4;
a = polyfit(T_tempM, Dz_tempM, PolyN);

t_tempM = roots(a);

t12 = t_tempM(3);

[T1,Y1] = ode45(@(t, y) droplet_wind_Modi(t, y, d),[0 t12],[vx0M vy0M vz0M]);

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
DxC = Dx1(end);
DyC = Dy1(end);
ErrorD = sqrt((DxO - DxC)^2 + (DyO - DyC)^2);
end


Condition = (ErrorD < ErrorMin);
end


if bbbbb == 1
p2 = (1 + (n2-1)*scale_p2)*p2_optimal; t1 = t11;
end
if bbbbb == 2
p2 = (1 - (n2-1)*scale_p2)*p2_optimal; t1 = t12;
end

p2_optimal = p2;

Condition11 = ErrorD > threshold;


if Optimal_find1 == 0 && Optimal_find2 == 0
    break
end
% pause();
end
