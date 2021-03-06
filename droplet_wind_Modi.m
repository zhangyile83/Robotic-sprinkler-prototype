function dy = droplet_wind_Modi(t, y, d)
% droplet function with different droplet diameter


g=9.8; % gravitation constant
fai=0.18; % fai is a constant related to the Reynolds number

pa=1.29; % pa is the atmospheric density
pw=1000; % pw is the water density
m=pi/6*pw*d^3; % mass of the droplet
k=fai*pa*d^2; % friction constant


w = 4.92;  % The unit is m/s
beta = pi/2 + pi/6; % try pi and pi/6
wx = w * cos(beta);
wy = w * sin(beta);


dy = zeros(3,1);    % a column vector
dy(1) = -(k/m) * sqrt((y(1) - wx)^2 + (y(2) - wy)^2 + y(3)^2) * (y(1)-wx);

dy(2) = -(k/m) * sqrt((y(1) - wx)^2 + (y(2) - wy)^2 + y(3)^2) * (y(2)-wy);

dy(3) = -(k/m) * sqrt((y(1) - wx)^2 + (y(2) - wy)^2 + y(3)^2) * y(3) - g;