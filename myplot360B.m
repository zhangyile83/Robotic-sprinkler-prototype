%

function myplot360B(A)
% [m,n]=size(A);
% if m==n
% else
%     A=reshape(A,sqrt(m),sqrt(m));
% end
 
% A = A(1:15, 16:45);
grid on;
A = A(11:30, 21:40);
colormap(flipud(bone));

h = pcolor(A);shading flat;
% set(h, 'EdgeColor', 'b');

% colorbar('eastoutside');
% hcb=colorbar;

% ylabel(hcb,'Precipitation level (m)',  'FontSize', 15);
set(gca,'xticklabel',[-25:25:50]);
% set(gca,'yticklabel',[-45:10:50]);
set(gca,'yticklabel',[-25:25:50]);

% set(gca,'yticklabel',[0:8:54 - 22 + 8 + 8  ]);

caxis([0 2*10^-3]);
title('Water Distribution',  'FontSize', 15);
% title({['Water Distribution with changing wind'];[' from w = 2m/s to 1 m/s']},  'FontSize', 20);
xlabel('x (m)', 'FontSize', 15);
ylabel('y (m)', 'FontSize', 15);