%% 画不含潮流和地转流的误差图，按照上下堆叠的方式

figure(1)
set(gcf,'position',[77.4,57.8,560,583.6])% 建立一个绘图图窗
%% 第1个坐标轴
ax = axes;% 生成一个坐标轴，命名为ax
surf(ax,sea_LON',sea_LAT',0*ones(size(sea_LON')),error_U');

% surf(ax,x',y',-40*ones(size(x')),Depth_core');% 在z=-1的位置绘制图层
colormap(ax,slanCL(194))% 给予ax坐标轴特定的colormap
axis([117.5 118.2 23 23.6 0 40]);% 调整坐标轴范围
shading interp
box off % 隐藏坐标轴轮廓
grid off % 隐藏坐标轴轮廓
hidden off % 隐藏坐标轴轮廓
set(ax,'zcolor','w','zticklabel','') % 隐藏z轴（将z轴设置为白色）
view(30,15) % 修改图片显示视角
xlabel('Longitude','Rotation',-8,'fontsize',12) % 调整xlabel方向  与x轴平行
ylabel('Latitude','Rotation',33,'fontsize',12) % 调整ylabel方向  与y轴平行
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %调整刻度和 

% title(date_save(1))
c1 = colorbar('horiz','AxisLocation','in','Limits',[-20 20]);% 纵向向colorbar，横向：horiz
set(c1,'Position',[0.7 0.30 0.2 0.01]);
set(get(c1,'title'),'string','\fontname{Times New Roman}Error (cm/s)','fontsize',10); % 这里展示了添加colorbar标题的另一个方法

%% 第2个坐标轴,没有地转流的u
ax2 = axes;
surf(ax2,sea_LON',sea_LAT',20*ones(size(sea_LON')),(U_sw-u)');
hold on 
% plot3(xx,yy,-ones(size(xx)),'k.','markersize',10)
colormap(ax2,slanCL(194))
axis([117.5 118.2 23 23.6 0 40]);% 调整坐标轴范围
shading interp
box off
axis off
hidden off
view(30,15)
% c2 = colorbar('horiz','Limits',[-40 40]);% 纵向colorbar，横向：horiz
c2 = colorbar('horiz');% 纵向colorbar，横向：horiz
set(c2,'Position',[0.7 0.6 0.2 0.01]);
set(get(c2,'title'),'string','\fontname{Times New Roman}Error (cm/s)','fontsize',10);
% title('\fontname{Times New Roman}U component error at 10:29 ','fontsize',14,'fontweight','bold','color','k')
% set(get(ax2,'title'),'string','U component error at 10:29','fontsize',10);
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %调整刻度和 
set(gca,'Clim',[-20 20]);
%% 第3个坐标轴,没有潮流的u
ax3 = axes;
surf(ax3,sea_LON',sea_LAT',40*ones(size(sea_LON')),(U_sw-uu)');
hold on 
% plot3(xx,yy,-ones(size(xx)),'k.','markersize',10)
colormap(ax3,slanCL(194))
% c3 = colorbar('horiz','Limits',[-40 40]);% 纵向colorbar，横向：horiz
c3 = colorbar('horiz');% 纵向colorbar，横向：horiz
set(c3,'Position',[0.7 0.9 0.2 0.01]);
set(get(c3,'title'),'string','\fontname{Times New Roman}Error (cm/s)','fontsize',10);
axis([117.5 118.2 23 23.6 0 40]);% 调整坐标轴范围
shading interp
box off
axis off
hidden off
view(30,15)
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %调整刻度和 
set(gca,'Clim',[-40 40]);
% set(get(ax3,'title'),'string','V component error at 10:29','fontsize',10);

% c3 = colorbar('horiz','Limits',[-20 20]);% 横向colorbar
% set(get(c3,'title'),'string','error (cm/s)','fontsize',10);
% set(c3,'Position',[0.76 0.88 0.18 0.02]);
% set(c3,'Limits',[-15 15]);  % 调整colorbar的范围
% print(figure(1),'-dpng','D:\博士研究\小论文\矢量流合成\投稿\图片\error_0932U','-r800')




%% 画不含潮流和地转流的误差图，按照上下堆叠的方式

figure(2)
set(gcf,'position',[77.4,57.8,560,583.6])% 建立一个绘图图窗
%% 第1个坐标轴
ax = axes;% 生成一个坐标轴，命名为ax
surf(ax,sea_LON',sea_LAT',0*ones(size(sea_LON')),error_V');

% surf(ax,x',y',-40*ones(size(x')),Depth_core');% 在z=-1的位置绘制图层
colormap(ax,slanCL(194))% 给予ax坐标轴特定的colormap
axis([117.5 118.2 23 23.6 0 40]);% 调整坐标轴范围
shading interp
box off % 隐藏坐标轴轮廓
grid off % 隐藏坐标轴轮廓
hidden off % 隐藏坐标轴轮廓
set(ax,'zcolor','w','zticklabel','') % 隐藏z轴（将z轴设置为白色）
view(30,15) % 修改图片显示视角
xlabel('Longitude','Rotation',-8,'fontsize',12) % 调整xlabel方向  与x轴平行
ylabel('Latitude','Rotation',33,'fontsize',12) % 调整ylabel方向  与y轴平行
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %调整刻度和 

% title(date_save(1))
c1 = colorbar('horiz','AxisLocation','in','Limits',[-20 20]);% 纵向向colorbar，横向：horiz
set(c1,'Position',[0.7 0.30 0.2 0.01]);
set(get(c1,'title'),'string','\fontname{Times New Roman}Error (cm/s)','fontsize',10); % 这里展示了添加colorbar标题的另一个方法

%% 第2个坐标轴,没有地转流的u
ax2 = axes;
surf(ax2,sea_LON',sea_LAT',20*ones(size(sea_LON')),(V_sw-v)');
hold on 
% plot3(xx,yy,-ones(size(xx)),'k.','markersize',10)
colormap(ax2,slanCL(194))
axis([117.5 118.2 23 23.6 0 40]);% 调整坐标轴范围
shading interp
box off
axis off
hidden off
view(30,15)
% c2 = colorbar('horiz','Limits',[-40 40]);% 纵向colorbar，横向：horiz
c2 = colorbar('horiz');% 纵向colorbar，横向：horiz
set(c2,'Position',[0.7 0.6 0.2 0.01]);
set(get(c2,'title'),'string','\fontname{Times New Roman}Error (cm/s)','fontsize',10);
% title('\fontname{Times New Roman}U component error at 10:29 ','fontsize',14,'fontweight','bold','color','k')
% set(get(ax2,'title'),'string','U component error at 10:29','fontsize',10);
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %调整刻度和 
set(gca,'Clim',[-20 20]);
%% 第3个坐标轴,没有潮流的u
ax3 = axes;
surf(ax3,sea_LON',sea_LAT',40*ones(size(sea_LON')),(V_sw-vv)');
hold on 
% plot3(xx,yy,-ones(size(xx)),'k.','markersize',10)
colormap(ax3,slanCL(194))
% c3 = colorbar('horiz','Limits',[-40 40]);% 纵向colorbar，横向：horiz
c3 = colorbar('horiz');% 纵向colorbar，横向：horiz
set(c3,'Position',[0.7 0.9 0.2 0.01]);
set(get(c3,'title'),'string','\fontname{Times New Roman}Error (cm/s)','fontsize',10);
axis([117.5 118.2 23 23.6 0 40]);% 调整坐标轴范围
shading interp
box off
axis off
hidden off
view(30,15)
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %调整刻度和 
set(gca,'Clim',[-40 40]);
% set(get(ax3,'title'),'string','V component error at 10:29','fontsize',10);

% c3 = colorbar('horiz','Limits',[-20 20]);% 横向colorbar
% set(get(c3,'title'),'string','error (cm/s)','fontsize',10);
% set(c3,'Position',[0.76 0.88 0.18 0.02]);
% set(c3,'Limits',[-15 15]);  % 调整colorbar的范围
% print(figure(2),'-dpng','D:\博士研究\小论文\矢量流合成\投稿\图片\error_0932V','-r800')














