function [ x,y, t] = droppoint( v0, gama, alpha, beta, d, w)
    vx = v0*cosd(alpha)*cosd(gama);
    vy = v0*cosd(alpha)*sind(gama);
    vz = v0*sind(alpha);
    [T_tempM1,Y_tempM1] = ode45(@(t, y) drop(t, y, d),[0:0.2:10],[vx vy vz]);
    vz_tempM1 = Y_tempM1(:,3);
    Dz_tempM1 = zeros(length(vz_tempM1),1);    
    for i = 2:length(T_tempM1)
        Dz_tempM1(i) = trapz(T_tempM1(1:i),vz_tempM1(1:i));
    end        
    PolyN = 4;
    a = polyfit(T_tempM1, Dz_tempM1, PolyN);    
    t_tempM1 = roots(a);        
    t = max(t_tempM1(3),t_tempM1(4));
    [T11,Y11] = ode45(@(t, y) droplet_wind_Modi(t, y, d, w, beta),[0 t],[vx vy vz]);    
    vx11 = Y11(:,1);
    vy11 = Y11(:,2);   
    x = trapz(T11(1:length(T11)),vx11(1:length(T11)));
    y = trapz(T11(1:length(T11)),vy11(1:length(T11)));        
end