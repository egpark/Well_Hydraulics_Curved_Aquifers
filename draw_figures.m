%% Point data plots (heights and orientation)
figure('color','w','position',[200 200 800 700])
% imagesc(imcomplement(Io))
hold on
gnorm=sqrt(grad_pnt(:,3).^2+grad_pnt(:,4).^2);
quiver(dat_pnt(:,1),dat_pnt(:,2), grad_pnt(:,4)./gnorm, -grad_pnt(:,3)./gnorm,0.125, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0);
quiver(dat_pnt(:,1),dat_pnt(:,2), -grad_pnt(:,4)./gnorm, grad_pnt(:,3)./gnorm,0.125, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0);
set(gca,'fontsize',24,'linewidth',2,'fontname','times new roman')
xlabel('\it\bfu\rm (m)','fontweight','bold','fontsize',32)
ylabel('\it\bfv\rm (m)','fontweight','bold','fontsize',32)
zlabel('\it\bfw\rm (m)','fontweight','bold','fontsize',32)
set(gca,'xtick',[1 50:50:1000],'xticklabel',[0 50:50:1000])
set(gca,'ytick',[1 50:50:1000],'yticklabel',[0 50:50:1000])
% set(gca,'ytick',[imH-ny:50:imH],'yticklabel',fliplr([0 50:50:ny]))
grid on
axis equal
axis tight
box on
set(gca,'ydir','normal')

plot3(dat_pnt(:,1),dat_pnt(:,2),dat_pnt(:,3),'ko','markersize',8,'markerfacecolor','w','linewidth',2)
plot(dat_pnt(:,1),dat_pnt(:,2),'ks','markersize',4,'markerfacecolor','k','linewidth',2)
view(3)
axis([1 308 1 218 0 max(dat_pnt(:,3))+eps])
for ii=1:17
    plot3([dat_pnt(ii,1) dat_pnt(ii,1)],[dat_pnt(ii,2) dat_pnt(ii,2)],[0 dat_pnt(ii,3)],'k--','linewidth',0.25)
end

%% Aquifer topography 2D
figure('position',[250 250 800 700],'color','w')
x = 1:1:nx;y = 1:1:ny;
contourf(x,y,T_est,25,'LineColor','none')
hold on
plot(dat_pnt(1:16,1),dat_pnt(1:16,2),'ko','markersize',8,'markerfacecolor','w','linewidth',2)
p=plot(dat_pnt(17,1),dat_pnt(17,2),'k','markersize',20,'markerfacecolor','w','linewidth',2);
p.Marker='pentagram';
hc=colorbar;
axis equal
axis tight
set(gca,'fontsize',24,'linewidth',2,'fontname','times new roman')
xlabel('\it\bfu \rm(m)','fontweight','bold','fontsize',32)
ylabel('\it\bfv \rm(m)','fontweight','bold','fontsize',32)
set(gca,'xtick',[1 50:50:1000],'xticklabel',[0 50:50:1000])
set(gca,'ytick',[1 50:50:1000],'yticklabel',[0 50:50:1000])
set(hc,'linewidth',2)
demcmap([min(T_est(:)) max(T_est(:))])

%% Analytical drawdowns (of this study)
figure('position',[250 250 800 700],'color','w')
hold on
% ddn1=expint(reshape((2*d_w).^2*5e-4/(4*1*1),ny,nx));
% contour(x,y,ddn1,[0.05 0.5:0.5:5],'--','LineColor',[0.75 0.75 0.75],'linewidth',1)

[~,h]=contourf(x,y,ddn,[0.1 0.5:0.5:10],'LineColor','w');
set(h,'facealpha',1)
box on
hc=colorbar;
axis equal
axis tight
set(gca,'fontsize',24,'linewidth',2,'fontname','times new roman')
xlabel('\it\bfu \rm(m)','fontweight','bold','fontsize',32)
ylabel('\it\bfv \rm(m)','fontweight','bold','fontsize',32)
set(gca,'xtick',[1 50:50:1000],'xticklabel',[0 50:50:1000])
set(gca,'ytick',[1 50:50:1000],'yticklabel',[0 50:50:1000])
set(hc,'linewidth',2)

%% Geodesic distances on aquifer topography
figure('position',[250 250 800 700],'color','w')
surf(x, y, T_est,reshape(2*dg_w,ny,nx),'FaceAlpha',1)
shading interp
hold on
plot3(dat_pnt(1:16,1),dat_pnt(1:16,2),dat_pnt(1:16,3),'ko','markersize',8,'markerfacecolor','w','linewidth',2)
p=plot3(dat_pnt(17,1),dat_pnt(17,2),dat_pnt(17,3),'k','markersize',20,'markerfacecolor','w','linewidth',2);
p.Marker='pentagram';
hc=colorbar;
set(hc,'location','northoutside')
axis equal
axis tight
set(gca,'fontsize',24,'linewidth',2,'fontname','times new roman')
xlabel('\it\bfu \rm(m)','fontweight','bold','fontsize',32)
ylabel('\it\bfv \rm(m)','fontweight','bold','fontsize',32)
zlabel('\it\bfw\rm (m)','fontweight','bold','fontsize',32)
set(gca,'xtick',[1 50:50:1000],'xticklabel',[0 50:50:1000])
set(gca,'ytick',[1 50:50:1000],'yticklabel',[0 50:50:1000])
set(hc,'linewidth',2)
camlight 
colormap(lines(7))
set(gca,"CameraPosition",[-423.0656 -1.0150e+03 1.6560e+03])
box on

%% Benchmark to the corresponding COMSOL result
% figure('position',[250 250 800 700],'color','w')
% data = readmatrix('Untitled.txt');
% xd = data(:, 1);
% yd = data(:, 2);
% zd = data(:, 3);
% values = data(:, 4);
% 
% [X,Y]=meshgrid(1:nx,1:ny);
% ddn_comsol = griddata(xd, yd, values, X, Y, 'cubic'); % You can change 'cubic' to another method if needed
% 
% [~,h]=contourf(x,y,ddn_comsol,[0.1 0.5:0.5:10],'LineColor','w');
% set(h,'facealpha',1)
% 
% box on
% % hold on
% % plot(dat_trn(:,1),dat_trn(:,2),'ko','markersize',8,'markerfacecolor','w','linewidth',2)
% hc=colorbar;
% axis equal
% axis tight
% set(gca,'fontsize',24,'linewidth',2,'fontname','times new roman')
% xlabel('\it\bfu \rm(m)','fontweight','bold','fontsize',32)
% ylabel('\it\bfv \rm(m)','fontweight','bold','fontsize',32)
% set(gca,'xtick',[1 50:50:1000],'xticklabel',[0 50:50:1000])
% set(gca,'ytick',[1 50:50:1000],'yticklabel',[0 50:50:1000])
% set(hc,'linewidth',2)


