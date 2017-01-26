%%
% %% Wind shifting Technique
% __ wind effect correction for 360 degree case
% __ Multi droplet diameter
% __ Multiple spraying circle
% __ corrections compute in advance
% __ Smooth version
clear;
tic
figure(1);
%%
%Configuration of water distribution plotting
RSL1 = 5; % Resolution in the water distribution figure: the unit is meter
Len = 300; % Water distribution plotting range(length): the unit is meter, note that the length and width are too big for our simulation case, therefore to show the properly show the results, tailoring is applied
Wid = 300; % Water distribution plotting range(width)
yic = 20;  % Increment at y direction in the water distribution graph to compensate the Matlab plotting discrepancy
GridsR1w = zeros(Wid/RSL1, Len/RSL1); % The matrix that collecting the water
Trejectory_width = 1; % This merely a plotting setting, the width here signify the width of the trejectory curve in the trejectory plotting

%%
%input
CC = 100; % CC here is a measure of water flowing volume per unit time, while the value need be verified by experiment data
alpha = 30;  % The angle of elevation of the nozzle
d_range = 0.005:0.001:0.010; % The droplet diameter range
MeanD = 0.01;  % Droplet mean diameter(should be the value defined in the range)
vertex_x=[-23;23;23;-23];
vertex_y=[-23;-23;23;23];
wind_speed = 0;
beta = 0;
distance_proportion_start = 1;
distance_proportion_step_length = 0.1;
distance_proportion_end = 1;
angle_step_length=20;

%%
%dubug information
count=1;
%%
%如果需要测试C code的结果，就把这个循环注释掉，把输入存入control_series
% %calculation of control series
[contour_distance, angle_start, angle_end]=calc_contour(vertex_x, vertex_y,angle_step_length);
for distance_proportion=distance_proportion_start:distance_proportion_step_length:distance_proportion_end
    index_gama = 1; % 181 in total
    %_______ Spraying loop at different angle in each spraying circle ______
    for gama=angle_start:angle_step_length:angle_end
        ticMark_total=tic;%debugging information
        distance=contour_distance(index_gama)*distance_proportion;
        v0=Characteristic(alpha, distance, MeanD); %Find the required flow speed for different spraying angles according
        ticMark_para=tic;
        [p1, p2] = Anti_wind_parameter_v4(v0, distance*cosd(gama), distance*sind(gama), gama, MeanD, wind_speed, beta);%_____________Find shifting parameters_________________
        fprintf('%d: %f %f\n',count, p1*v0, gama*p2-360);
        control_series(count,1)=p1*v0;
        control_series(count,2)=gama*p2-360;
        control_series_without_shifting(count,1)=v0;
        control_series_without_shifting(count,2)=gama;
        index_gama=index_gama+1;
        count=count+1;
    end    
end

%%
%plotting contour
for i=[1,3,4,6]
subplot(2,3,i);
plot([vertex_x;vertex_x(1)],[vertex_y;vertex_y(1)],'-.ok', 'LineWidth',2,'MarkerFaceColor','k');
grid on;
hold on;
end
% subplot(2,3,2);
% hold on;
% subplot(2,3,5);
% hold on;

%%
%Spraying simulation without Wind shifting parameters
for count=1:size(control_series_without_shifting,1)
    v0 = control_series_without_shifting(count,1);
    gama = control_series_without_shifting(count,2);
    for ii = 1:length(d_range)
        d = d_range(ii);
        vx0 = v0*cosd(alpha)*cosd(gama);
        vy0 = v0*cosd(alpha)*sind(gama);
        vz0 = v0*sind(alpha);
        %%
        %Dropping time for droplet
        [T_temp,Y_temp] = ode45(@(t, y) drop(t, y, d),[0 10],[vx0 vy0 vz0]);
        vz_temp = Y_temp(:,3);
        Dz_temp = zeros(length(vz_temp),1);
        for i = 2:length(T_temp)
            Dz_temp(i) = trapz(T_temp(1:i),vz_temp(1:i));
        end
        PolyN = 4;
        a = polyfit(T_temp, Dz_temp, PolyN);
        t_temp = roots(a);
        t3 = max(t_temp(3),t_temp(4));
        %%
        %Trajectory simulation under wind condition without shifting
        [T1,Y1] = ode45(@(t, y) droplet_wind_Modi(t, y, d, wind_speed, beta),[0 t3],[vx0 vy0 vz0]);
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
        drop_x=Dx1(end);
        drop_y=Dy1(end);
        drop_z=Dz1(end);
        %%
        %For the plotting purpose, we only show the trajectory of the droplet with mean diameter
        if ii == find(d_range == MeanD)            
            subplot(2,3,3);
            pause(0.001);
            axis([-40 40 -40 40]);
            title({['Dropping points without shifting']}, 'FontSize', 15);
            plot(drop_x,drop_y,'d');
            subplot(2,3,1);
            pause(0.001);
            axis([-40 40 -40 40]);
            title({['Trajectory without shifting']}, 'FontSize', 15);
            %%
            %Plotting trajectory
            plot3(Dx1,Dy1,Dz1,'LineWidth',Trejectory_width);
            xlabel('x (m)', 'FontSize', 15);
            ylabel('y (m)', 'FontSize', 15);
            zlabel('z (m)', 'FontSize', 15);
            set(gca,'XGrid','on')
            set(gca,'YGrid','on')
        end
        Dd1(ii) = Dx1(end);
        Dxtemp1w = round(Dx1(end)/RSL1);
        Dytemp1w = round(Dy1(end)/RSL1);
        %%
        %             %The following setting is determined by the nozzle configurations
        %             if spraying_C<=20
        %                 CC = 60;
        %             end
        %             if spraying_C<= 15
        %                 CC = 30;
        %             end
        %%
        %Distribution of the number of droplet with different diameters
        Vd_temp = (pi/6)*pi*d^3*CC*distribution(d,MeanD,0.0015); %Normal distribution assumption is used
        Vd_temp_1 = Vd_temp*4;
        %Accumulation of water at different bucket
        GridsR1w(Dxtemp1w + Wid/2/RSL1, Dytemp1w + yic) = GridsR1w(Dxtemp1w + Wid/2/RSL1, Dytemp1w + yic) + Vd_temp_1;
    end
    subplot(2,3,2)
    myplot360B(GridsR1w');
end

%%
%Spraying simulationwith with shifting parameters
for count=1:size(control_series,1)
    v0 = control_series(count,1);
    gama = control_series(count,2);
    for ii = 1:length(d_range)
        d = d_range(ii);
        vx0 = v0*cosd(alpha)*cosd(gama);
        vy0 = v0*cosd(alpha)*sind(gama);
        vz0 = v0*sind(alpha);
        %%
        %Dropping time for droplet with different diameters
        [T_temp,Y_temp] = ode45(@(t, y) drop(t, y, d),[0 10],[vx0 vy0 vz0]);
        vz_temp = Y_temp(:,3);
        Dz_temp = zeros(length(vz_temp),1);
        for i = 2:length(T_temp)
            Dz_temp(i) = trapz(T_temp(1:i),vz_temp(1:i));
        end
        PolyN = 4;
        a = polyfit(T_temp, Dz_temp, PolyN);
        t_temp = roots(a);
        t3 = max(t_temp(3),t_temp(4));
        %%
        %Trajectory simulation under wind condition after shifting
        [T1,Y1] = ode45(@(t, y) droplet_wind_Modi(t, y, d, wind_speed, beta),[0 t3],[vx0 vy0 vz0]);
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
        drop_x=Dx1(end);
        drop_y=Dy1(end);
        drop_z=Dz1(end);
        %%
        %For the plotting purpose, we only show the trajectory of the droplet with mean diameter
        if ii == find(d_range == MeanD)            
            subplot(2,3,6);
            pause(0.001);
            axis([-40 40 -40 40]);
            title({['Dropping points after shifting']}, 'FontSize', 15);
            plot(drop_x,drop_y,'d');
            subplot(2,3,4);
            pause(0.001);
            axis([-40 40 -40 40]);
            title({['Trajectory after shifting']}, 'FontSize', 15);
            %%
            %Plotting trajectory
            plot3(Dx1,Dy1,Dz1,'LineWidth',Trejectory_width);            
            xlabel('x (m)', 'FontSize', 15);
            ylabel('y (m)', 'FontSize', 15);
            zlabel('z (m)', 'FontSize', 15);
            set(gca,'XGrid','on')
            set(gca,'YGrid','on')
        end
        Dd1(ii) = Dx1(end);
        Dxtemp1w = round(Dx1(end)/RSL1);
        Dytemp1w = round(Dy1(end)/RSL1);
        %%
        %             %The following setting is determined by the nozzle configurations _____
        %             if spraying_C<=20
        %                 CC = 60;
        %             end
        %             if spraying_C<= 15
        %                 CC = 30;
        %             end
        %%
        %Distribution of the number of droplet with different diameter
        Vd_temp = (pi/6)*pi*d^3*CC*distribution(d,MeanD,0.0015); %Normal distribution assumption is used
        Vd_temp_1 = Vd_temp*4;
        %Accumulation of water at different bucket
        GridsR1w(Dxtemp1w + Wid/2/RSL1, Dytemp1w + yic) = GridsR1w(Dxtemp1w + Wid/2/RSL1, Dytemp1w + yic) + Vd_temp_1;
    end
    subplot(2,3,5)
    myplot360B(GridsR1w');
end
toc