function s = FindSpeed(V, x)
% this function is used to find the required speed for a given
% speed data in form of:
%DD(i) is the distance of a droplet at the speed of i, i = 100 means v = 100 m/s
% 
v0temp = 1:200;
a=polyfit(v0temp,V,4);
DD =polyval(a,v0temp);


temp = DD - x;
s = find(min(abs(temp)) == abs(temp));
% s = s/10;