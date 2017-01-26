function [y] = Characteristic(alpha, x)

% x is the sequence of distance
% y(i) is the corresponding velocity to reach the distance x(i)

v0temp = 1:200;
for v0 = v0temp
vx0 = v0*cos(alpha);
vy0 = 0;
vz0 = v0*sin(alpha);

[T_temp,Y_temp] = ode45(@drop,[0 10],[vx0 vy0 vz0]);

vx_temp = Y_temp(:,1);
vz_temp = Y_temp(:,3);

for i = 2:length(T_temp)
Dx_temp(i) = trapz(T_temp(1:i),vx_temp(1:i));
Dz_temp(i) = trapz(T_temp(1:i),vz_temp(1:i));
end
index_temp = AlterSign(Dz_temp());
V(v0) = Dx_temp(index_temp);
end

%%%%%%%%%%% using a simpler way to find the require velocity %%%%%%%
DD = V;   %DD(i) is the distance of a droplet at the speed of i, i = 100 means v = 100 m/s
y = zeros(1, length(x));

for i = 1:length(x)
    y(i) = FindSpeed(DD, x(i));
end

%%%%%%%%%%% using interpolation to find the required velocity %%%%%%%
% 
% PolyN = 4;
% a = polyfit(v0temp,V,PolyN);
% 
% for i = 1:length(x)
%     atemp = a;
%     atemp(PolyN + 1) = atemp(PolyN + 1) - x(i);
%     
%     ytemp = roots(atemp);
%     y(i) = ytemp(end);
% end