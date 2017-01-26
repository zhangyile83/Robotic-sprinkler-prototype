function dy = droplet(t, y, d)
% d = 0.02;
g=9.8; % gravitation constant
fai=0.18; % fai is a constant related to the Reynolds number

pa=1.29; % pa is the atmospheric density
pw=1000; % pw is the water density
m=pi/6*pw*d^3; % mass of the droplet
k=fai*pa*d^2; % friction constant

dy = zeros(3,1);    % a column vector
dy(1) = -(k/m) * sqrt(y(1)^2 + y(2)^2 + y(3)^2) * y(1);

dy(2) = -(k/m) * sqrt(y(1)^2 + y(2)^2 + y(3)^2) * y(2);

dy(3) = -(k/m) * sqrt(y(1)^2 + y(2)^2 + y(3)^2) * y(3) - g;