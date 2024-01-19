% t_model,T_E_3R(range_model)
% L = length(range_model);
% Hov = zeros(101,L);
% Hov(26,:) = T_C_3R(range_model)';
% Hov(76,:) = T_E_3R(range_model)';
% dx = (T_E_3R(range_model)' - T_C_3R(range_model)')/50;
% for i = 1:101
%     Hov(i,:) = Hov(26,:) + (i-26)*dx;
% end
% figure
% [xx,yy] = meshgrid(t_model,1:101);
% contourf(yy,xx,Hov,30,'linestyle','none')
% caxis([-2,2])
% colorbar

L = length(range_model);
Hov = zeros(150,L);
for i = 1:50
    Hov(i,:) = h_W_3R(range_model)'/30;
end
for i = 51:100
    Hov(i,:) = T_C_3R(range_model)';
end
for i = 101:150
    Hov(i,:) = T_E_3R(range_model)';
end
for j = 1:L
    Hov(:,j) = smooth(Hov(:,j),30);
end
figure
colormap(jet)
[xx,yy] = meshgrid(t_model,1:150);
contourf(yy,xx,Hov,30,'linestyle','none')
caxis([-2,2])
colorbar
