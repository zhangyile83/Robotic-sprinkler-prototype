function [y] = Characteristic(alpha, x, d)
% x is the sequence of distance
% y(i) is the corresponding velocity to reach the distance x(i)
scale=0.2;
v=1:scale:150;
persistent DD;
if isempty(DD)
    DD= zeros(1,size(v,2));   %DD(i) is the distance of a droplet at the speed of i, i = 100 means v = 100 m/s    
    for i=1:size(v,2)
        v0 = v(i);
        [DD(i),~]=droppoint( v0, 360, alpha, 0, d, 0);        
    end
end
temp = DD - x;
s = min(abs(temp)) == abs(temp);
y=v(s);

