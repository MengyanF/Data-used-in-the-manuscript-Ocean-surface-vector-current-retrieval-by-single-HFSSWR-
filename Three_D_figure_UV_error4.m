%% �����������͵�ת�������ͼ���������¶ѵ��ķ�ʽ

figure(1)
set(gcf,'position',[77.4,57.8,560,583.6])% ����һ����ͼͼ��
%% ��1��������
ax = axes;% ����һ�������ᣬ����Ϊax
surf(ax,sea_LON',sea_LAT',0*ones(size(sea_LON')),error_U');

% surf(ax,x',y',-40*ones(size(x')),Depth_core');% ��z=-1��λ�û���ͼ��
colormap(ax,slanCL(194))% ����ax�������ض���colormap
axis([117.5 118.2 23 23.6 0 40]);% ���������᷶Χ
shading interp
box off % ��������������
grid off % ��������������
hidden off % ��������������
set(ax,'zcolor','w','zticklabel','') % ����z�ᣨ��z������Ϊ��ɫ��
view(30,15) % �޸�ͼƬ��ʾ�ӽ�
xlabel('Longitude','Rotation',-8,'fontsize',12) % ����xlabel����  ��x��ƽ��
ylabel('Latitude','Rotation',33,'fontsize',12) % ����ylabel����  ��y��ƽ��
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %�����̶Ⱥ� 

% title(date_save(1))
c1 = colorbar('horiz','AxisLocation','in','Limits',[-20 20]);% ������colorbar������horiz
set(c1,'Position',[0.7 0.30 0.2 0.01]);
set(get(c1,'title'),'string','\fontname{Times New Roman}Error (cm/s)','fontsize',10); % ����չʾ�����colorbar�������һ������

%% ��2��������,û�е�ת����u
ax2 = axes;
surf(ax2,sea_LON',sea_LAT',20*ones(size(sea_LON')),(U_sw-u)');
hold on 
% plot3(xx,yy,-ones(size(xx)),'k.','markersize',10)
colormap(ax2,slanCL(194))
axis([117.5 118.2 23 23.6 0 40]);% ���������᷶Χ
shading interp
box off
axis off
hidden off
view(30,15)
% c2 = colorbar('horiz','Limits',[-40 40]);% ����colorbar������horiz
c2 = colorbar('horiz');% ����colorbar������horiz
set(c2,'Position',[0.7 0.6 0.2 0.01]);
set(get(c2,'title'),'string','\fontname{Times New Roman}Error (cm/s)','fontsize',10);
% title('\fontname{Times New Roman}U component error at 10:29 ','fontsize',14,'fontweight','bold','color','k')
% set(get(ax2,'title'),'string','U component error at 10:29','fontsize',10);
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %�����̶Ⱥ� 
set(gca,'Clim',[-20 20]);
%% ��3��������,û�г�����u
ax3 = axes;
surf(ax3,sea_LON',sea_LAT',40*ones(size(sea_LON')),(U_sw-uu)');
hold on 
% plot3(xx,yy,-ones(size(xx)),'k.','markersize',10)
colormap(ax3,slanCL(194))
% c3 = colorbar('horiz','Limits',[-40 40]);% ����colorbar������horiz
c3 = colorbar('horiz');% ����colorbar������horiz
set(c3,'Position',[0.7 0.9 0.2 0.01]);
set(get(c3,'title'),'string','\fontname{Times New Roman}Error (cm/s)','fontsize',10);
axis([117.5 118.2 23 23.6 0 40]);% ���������᷶Χ
shading interp
box off
axis off
hidden off
view(30,15)
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %�����̶Ⱥ� 
set(gca,'Clim',[-40 40]);
% set(get(ax3,'title'),'string','V component error at 10:29','fontsize',10);

% c3 = colorbar('horiz','Limits',[-20 20]);% ����colorbar
% set(get(c3,'title'),'string','error (cm/s)','fontsize',10);
% set(c3,'Position',[0.76 0.88 0.18 0.02]);
% set(c3,'Limits',[-15 15]);  % ����colorbar�ķ�Χ
% print(figure(1),'-dpng','D:\��ʿ�о�\С����\ʸ�����ϳ�\Ͷ��\ͼƬ\error_0932U','-r800')




%% �����������͵�ת�������ͼ���������¶ѵ��ķ�ʽ

figure(2)
set(gcf,'position',[77.4,57.8,560,583.6])% ����һ����ͼͼ��
%% ��1��������
ax = axes;% ����һ�������ᣬ����Ϊax
surf(ax,sea_LON',sea_LAT',0*ones(size(sea_LON')),error_V');

% surf(ax,x',y',-40*ones(size(x')),Depth_core');% ��z=-1��λ�û���ͼ��
colormap(ax,slanCL(194))% ����ax�������ض���colormap
axis([117.5 118.2 23 23.6 0 40]);% ���������᷶Χ
shading interp
box off % ��������������
grid off % ��������������
hidden off % ��������������
set(ax,'zcolor','w','zticklabel','') % ����z�ᣨ��z������Ϊ��ɫ��
view(30,15) % �޸�ͼƬ��ʾ�ӽ�
xlabel('Longitude','Rotation',-8,'fontsize',12) % ����xlabel����  ��x��ƽ��
ylabel('Latitude','Rotation',33,'fontsize',12) % ����ylabel����  ��y��ƽ��
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %�����̶Ⱥ� 

% title(date_save(1))
c1 = colorbar('horiz','AxisLocation','in','Limits',[-20 20]);% ������colorbar������horiz
set(c1,'Position',[0.7 0.30 0.2 0.01]);
set(get(c1,'title'),'string','\fontname{Times New Roman}Error (cm/s)','fontsize',10); % ����չʾ�����colorbar�������һ������

%% ��2��������,û�е�ת����u
ax2 = axes;
surf(ax2,sea_LON',sea_LAT',20*ones(size(sea_LON')),(V_sw-v)');
hold on 
% plot3(xx,yy,-ones(size(xx)),'k.','markersize',10)
colormap(ax2,slanCL(194))
axis([117.5 118.2 23 23.6 0 40]);% ���������᷶Χ
shading interp
box off
axis off
hidden off
view(30,15)
% c2 = colorbar('horiz','Limits',[-40 40]);% ����colorbar������horiz
c2 = colorbar('horiz');% ����colorbar������horiz
set(c2,'Position',[0.7 0.6 0.2 0.01]);
set(get(c2,'title'),'string','\fontname{Times New Roman}Error (cm/s)','fontsize',10);
% title('\fontname{Times New Roman}U component error at 10:29 ','fontsize',14,'fontweight','bold','color','k')
% set(get(ax2,'title'),'string','U component error at 10:29','fontsize',10);
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %�����̶Ⱥ� 
set(gca,'Clim',[-20 20]);
%% ��3��������,û�г�����u
ax3 = axes;
surf(ax3,sea_LON',sea_LAT',40*ones(size(sea_LON')),(V_sw-vv)');
hold on 
% plot3(xx,yy,-ones(size(xx)),'k.','markersize',10)
colormap(ax3,slanCL(194))
% c3 = colorbar('horiz','Limits',[-40 40]);% ����colorbar������horiz
c3 = colorbar('horiz');% ����colorbar������horiz
set(c3,'Position',[0.7 0.9 0.2 0.01]);
set(get(c3,'title'),'string','\fontname{Times New Roman}Error (cm/s)','fontsize',10);
axis([117.5 118.2 23 23.6 0 40]);% ���������᷶Χ
shading interp
box off
axis off
hidden off
view(30,15)
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %�����̶Ⱥ� 
set(gca,'Clim',[-40 40]);
% set(get(ax3,'title'),'string','V component error at 10:29','fontsize',10);

% c3 = colorbar('horiz','Limits',[-20 20]);% ����colorbar
% set(get(c3,'title'),'string','error (cm/s)','fontsize',10);
% set(c3,'Position',[0.76 0.88 0.18 0.02]);
% set(c3,'Limits',[-15 15]);  % ����colorbar�ķ�Χ
% print(figure(2),'-dpng','D:\��ʿ�о�\С����\ʸ�����ϳ�\Ͷ��\ͼƬ\error_0932V','-r800')














