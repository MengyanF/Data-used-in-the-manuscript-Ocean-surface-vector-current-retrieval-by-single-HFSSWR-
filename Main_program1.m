%% 进行单站天地波雷达矢量流合成算法的主程序
%% 作者：冯梦延
%% 时间：2023.3.15
clear;
clc;
%% 加载风场数据和海面高度场数据
% 风场数据注意事项
%   风场数据选择了从2017年3月1日到4月31日每小时的ECMWF数据
%   由于不需要考虑风场的时延特征，且考虑到 LT = UT + 7.5h
%   4月4日9：32和9：44的风场数据用 V10(:,:,819) 和 U10(:,:,819)
%   10：19和10：29的用 V10(:,:,820) 和 U(:,:,820)
% 海面高度注意事项
%   海面高度的时间分辨率只有每日，由于台湾海峡常年存在一个东北
%   方向10cm/s的地转流，因此该时间分辨率的海面高度数据足够

% load wind_data.mat
load wind_data_4_4.mat
load SSH.mat
% 海水深度
lat_h = ncread('D:\博士研究\小论文\矢量流合成\program\海面高度\gebco_2022_n26.0_s20.0_w115.0_e122.0.nc','lat');
lon_h = ncread('D:\博士研究\小论文\矢量流合成\program\海面高度\gebco_2022_n26.0_s20.0_w115.0_e122.0.nc','lon');
h = ncread('D:\博士研究\小论文\矢量流合成\program\海面高度\gebco_2022_n26.0_s20.0_w115.0_e122.0.nc','elevation');

%% 依据海面高度计算地转流

SSH = (SSH(:,:,4) + SSH(:,:,3))./2;
[x,y] = size(SSH);
% 南北分量
for j =1 : y
    for i = 1 : x-1
        f = double(7.272e-5 * 2 * sind(latitude_SSH(j)));
        SSH(i+1,j)- SSH(i,j)
        double(((longitude_SSH(i+1)-longitude_SSH(i))*100000))
        V_f(i,j) = (SSH(i+1,j)- SSH(i,j))./ double(((longitude_SSH(i+1)-longitude_SSH(i))*100000)) * 9.8 ./ f;
        V_f(i,j)
    end
        V_f(x,j) = 0;
end
% 东西分量
for i =1 : x
    for j = 1 : y-1
        f = double(7.272e-5 * 2 * sind(latitude_SSH(j)));
        U_f(i,j) = -1 * (SSH(i,j+1)- SSH(i,j))./ double(((latitude_SSH(j+1)-latitude_SSH(j))*100000)) * 9.8 ./ f;
    end
        U_f(i,y) = 0;
end
% 画图
figure
S = load('EastCoastline.dat');
plot(S(:,1),S(:,2),'k-','linewidth',0.5);
hold on
for i = 1 : length(V_f(:,1))
    for j =  1 : length(V_f(1,:))
        if isnan(V_f(i,j)) || isnan(U_f(i,j))
            continue
        else
            quiver(longitude_SSH(i),latitude_SSH(j),0.8*U_f(i,j),0.8*V_f(i,j),'r','AutoScale','off','MaxHeadSize',2);hold on
        end    
    end
end
axis ([115 120 21 25])
% axis ([117 118.8 22.4 24.2])
set(gcf,'Position',[100 100 800 400]);

% 地转流速 
cd_speed = sqrt(U_f.^2 + V_f.^2)*100;
figure
pcolor(longitude_SSH,latitude_SSH,cd_speed');hold on
shading interp
hc=colorbar;
set(gca,'Clim',[0 60]);  % 调整colorbar的范围
S = load('EastCoastline.dat');
plot(S(:,1),S(:,2),'k-','linewidth',0.5);
hold on
% axis ([115 120 21 25])
axis ([117.5 118.2 23 23.6]) % 核心区域
set(gcf,'Position',[100 100 800 400]);
ylabel(hc,'\fontname{Times New Roman}Current speed (cm/s)','fontweight','bold','fontsize',14);
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %调整刻度和 
xlabel('\fontname{Times New Roman}Longitude (°E)','fontsize',14,'fontweight','bold','color','k');
ylabel('\fontname{Times New Roman}Latitude (°N)','fontsize',14,'fontweight','bold','color','k')
title('\fontname{Times New Roman}2017/04/04','fontsize',14,'fontweight','bold','color','k')




%  海面高度画图
figure
% 东山经纬度
Sta1Lon=117.4863;Sta1Lat=23.6575; 
% 赤湖经纬度
SiteLatChS = 24.03702;SiteLonChS = 117.9025;
S = load('EastCoastline.dat');
plot(S(:,1),S(:,2),'k-','linewidth',0.5);
hold on
pcolor(longitude_SSH,latitude_SSH,SSH');hold on
shading interp
hc=colorbar;
axis ([115 120 21 25]) % 核心区域
colormap(jet)
set(gca,'Clim',[0.3 0.9]);  % 调整colorbar的范围
plot(Sta1Lon,Sta1Lat,'p','linewidth',2);
text(Sta1Lon,Sta1Lat,'\fontname{Times New Roman}Dongshan','fontsize',12,'Fontname','Times New Man','FontWeight','Bold');
plot(SiteLonChS,SiteLatChS,'p','linewidth',2);
text(SiteLonChS,SiteLatChS,'\fontname{Times New Roman}Chihu','fontsize',12,'Fontname','Times New Man','FontWeight','Bold');
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %调整刻度和 
xlabel('\fontname{Times New Roman}Longitude (°E)','fontsize',14,'fontweight','bold','color','k');
ylabel('\fontname{Times New Roman}Latitude (°N)','fontsize',14,'fontweight','bold','color','k')
title('\fontname{Times New Roman}2017/04/04 \fontname{宋体}海面高度','fontsize',14,'fontweight','bold','color','k')
ylabel(hc,'\fontname{Times New Roman}SSH (m)','fontweight','bold','fontsize',14);



%% 风生流
% U_wind = U10(:,:,819);
% V_wind = V10(:,:,819);
%% Ekman流
pa = 1.29; % 空气密度
pw = 1.05*1e3; % 海水密度
% 水深
h = 20;
% 涡流粘度,单位 m^2/s
Az = 0.01;
z=0;
clear j
% latitude_wind = double(latitude_wind);
% longitude_wind = double(longitude_wind);
% V_wind = double(V_wind);
for i = 1 : length(latitude_wind)
    for k = 1 : length(longitude_wind)
        VE(k,i) = Func_Ekman(U_wind(k,i), V_wind(k,i), latitude_wind(i), Az, h, z);
    end
end
Ue = real(VE);
Ve = imag(VE);

% Ekman流速
figure
pcolor(longitude_wind,latitude_wind,abs(VE)');hold on
shading interp
hc=colorbar;
% set(gca,'Clim',[0 1]);  % 调整colorbar的范围
S = load('EastCoastline.dat');
plot(S(:,1),S(:,2),'k-','linewidth',0.5);
hold on
axis ([115 120 21 25]);
set(gcf,'Position',[100 100 800 400]);
ylabel(hc,'\fontname{Times New Roman}Ekman speed (m)','fontweight','bold','fontsize',14);
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %调整刻度和 
xlabel('\fontname{Times New Roman}Longitude (°E)','fontsize',14,'fontweight','bold','color','k');
ylabel('\fontname{Times New Roman}Latitude (°N)','fontsize',14,'fontweight','bold','color','k');
title('\fontname{Times New Roman}2017/04/04','fontsize',14,'fontweight','bold','color','k');

% Ekman 矢量
figure
S = load('EastCoastline.dat');
plot(S(:,1),S(:,2),'k-','linewidth',0.5);
hold on
for i = 1 : length(Ve(:,1))
    for j =  1 : length(Ve(1,:))
        if isnan(Ve(i,j)) || isnan(Ue(i,j))
            continue
        else
            quiver(longitude_wind(i),latitude_wind(j),Ue(i,j),Ve(i,j),'r','AutoScale','off','MaxHeadSize',2);hold on
            quiver(longitude_wind(i),latitude_wind(j),0.005*U_wind(i,j),0.005*V_wind(i,j),'b','AutoScale','off','MaxHeadSize',2);hold on            
        end    
    end
end
axis ([115 120 21 25])
% axis ([117 118.8 22.4 24.2])
set(gcf,'Position',[100 100 800 400]);

%% Stokes流
% stokes流计划采用参数化方案
% JONSWAP谱
F =100e3;% 风区大小
z = -0.5;
wind_speed = sqrt(U_wind.^2 +  V_wind.^2);
for i = 1 : length(wind_speed(:,1))
    for j =  1 : length(wind_speed(1,:))
        [VS0(i,j), VSZ(i,j)] = Func_Stokes(wind_speed(i,j), F, z);
    end
end
for i = 1 : length(wind_speed(:,1))
    for j =  1 : length(wind_speed(1,:))
        if U_wind(i,j) > 0 & V_wind(i,j) > 0
            aoe(i,j) = atand(V_wind(i,j)/U_wind(i,j));
        elseif U_wind(i,j) > 0 & V_wind(i,j) < 0
            aoe(i,j) = atand(V_wind(i,j)/U_wind(i,j));
        elseif U_wind(i,j) < 0 & V_wind(i,j) < 0
            aoe(i,j) = atand(V_wind(i,j)/U_wind(i,j))-180;
        elseif U_wind(i,j) < 0 & V_wind(i,j) > 0
            aoe(i,j) = atand(V_wind(i,j)/U_wind(i,j))+180;
        end
    end
end
U_VS0 = VS0 .* cosd(aoe);
V_VS0 = VS0 .* sind(aoe);
U_VSZ = VSZ .* cosd(aoe);
V_VSZ = VSZ .* sind(aoe);

% Stokes 矢量
figure
S = load('EastCoastline.dat');
plot(S(:,1),S(:,2),'k-','linewidth',0.5);
hold on
for i = 1 : length(V_VS0(:,1))
    for j =  1 : length(V_VS0(1,:))
        if isnan(V_VS0(i,j)) || isnan(U_VS0(i,j))
            continue
        else
            quiver(longitude_wind(i),latitude_wind(j),U_VS0(i,j),V_VS0(i,j),'r','AutoScale','off','MaxHeadSize',2);hold on
            quiver(longitude_wind(i),latitude_wind(j),0.005*U_wind(i,j),0.005*V_wind(i,j),'b','AutoScale','off','MaxHeadSize',2);hold on            
        end    
    end
end
axis ([115 120 21 25])
% axis ([117 118.8 22.4 24.2])
set(gcf,'Position',[100 100 800 400]);


%% 天地波雷达径向流
% load('D:\博士研究\小论文\矢量流合成\program\径向流\scatter_point_information_9_32_new_eta_0.5853_GP_GIA_all_information.mat')
% load('D:\博士研究\小论文\矢量流合成\program\径向流\scatter_point_information_9_44_new_eta_0.5609_GP_GIA_all_information.mat')
% load('D:\博士研究\小论文\矢量流合成\program\径向流\scatter_point_information_10_19_new_eta_0.4456_GP_GIA_all_information.mat')
load('D:\博士研究\小论文\矢量流合成\program\径向流\scatter_point_information_10_29_new_eta_0.3995_GP_GIA_all_information.mat')

% 东山经纬度
Sta1Lon=117.4863;Sta1Lat=23.6575; 
% 赤湖经纬度
SiteLatChS = 24.03702;SiteLonChS = 117.9025;
% 武大经纬度
SiteLatWhU = 30.541093;SiteLonWhU = 114.360734;

fr = fieldnames(w_20170404093225);
scale=2.3*1e-3;
figure
for j = 1 : length(fr)
    data = w_20170404093225_smooth.(fr{j});
    vb_IRI = data(8,:);
    dir_IRI = data(9,:);
    v_sur_Dir = data(11,:);
    Sea_scatter_lon = data(3,:);
    Sea_scatter_lat = data(2,:);
    v_sur_speed = data(12,:);
    for i = 1 : length(vb_IRI)
        if vb_IRI(i) < 0 && vb_IRI(i) > -50 && v_sur_speed(i)<50
            dir = dir_IRI(i) + 180;
            % 红色是反演值，蓝色是观测值
            quiver(Sea_scatter_lon(i),Sea_scatter_lat(i),scale*abs(vb_IRI(i)).*sind(dir),scale*abs(vb_IRI(i)).*cosd(dir),'r','AutoScale','off','MaxHeadSize',2);hold on
            % 地波矢量流在 dir_IRI 上的投影，dir_IRI与v_sur_Dir之间是锐角，v_sur大于0 
            v_sur = v_sur_speed(i) * cosd(dir_IRI(i) - v_sur_Dir(i));
            if v_sur > 0 % 说明地波的径向流方向为 dir_IRI(i)
                quiver(Sea_scatter_lon(i),Sea_scatter_lat(i),scale*abs(v_sur).*sind(dir_IRI(i)),scale*abs(v_sur).*cosd(dir_IRI(i)),'b','AutoScale','off','MaxHeadSize',2);
                hold on
            elseif v_sur < 0 % 说明地波的径向流方向为 dir
                quiver(Sea_scatter_lon(i),Sea_scatter_lat(i),scale*abs(v_sur).*sind(dir),scale*abs(v_sur).*cosd(dir),'b','AutoScale','off','MaxHeadSize',2);
                hold on
            end
        elseif vb_IRI(i) > 0 && vb_IRI(i) < 50 && v_sur_speed(i)<50
            dir = dir_IRI(i);
            quiver(Sea_scatter_lon(i),Sea_scatter_lat(i),scale*abs(vb_IRI(i)).*sind(dir),scale*abs(vb_IRI(i)).*cosd(dir),'r','AutoScale','off','MaxHeadSize',2);hold on
            v_sur = v_sur_speed(i) * cosd(dir_IRI(i) - v_sur_Dir(i));
            if v_sur > 0  % 说明地波的径向流方向为 dir
                quiver(Sea_scatter_lon(i),Sea_scatter_lat(i),scale*abs(v_sur).*sind(dir),scale*abs(v_sur).*cosd(dir),'b','AutoScale','off','MaxHeadSize',2);
                hold on
            elseif v_sur < 0  % 说明地波的径向流方向为 dir+180
                quiver(Sea_scatter_lon(i),Sea_scatter_lat(i),scale*abs(v_sur).*sind(dir+180),scale*abs(v_sur).*cosd(dir+180),'b','AutoScale','off','MaxHeadSize',2);
                hold on
            end
        else
            continue
        end
        if i==1 && j==1
            legend('\fontname{Times New Roman}V_E of sky-surface wave','\fontname{Times New Roman}V_E of surface wave','fontsize',10) 
        end
    end
     clear data vb_IRI dir_IRI v_sur_Dir
end
S = load('EastCoastline.dat');
plot(S(:,1),S(:,2),'k-','linewidth',0.5);
hold on
% plot(SiteLon_LoH,SiteLat_LoH,'r*','linewidth',3);
% text(SiteLon_LoH,SiteLat_LoH,'Longhai','fontsize',12,'Fontname','Times New Man','FontWeight','Bold')
plot(SiteLonChS,SiteLatChS,'p','linewidth',2);
text(SiteLonChS,SiteLatChS,'\fontname{Times New Roman}Chihu','fontsize',12,'Fontname','Times New Man','FontWeight','Bold');
plot(Sta1Lon,Sta1Lat,'p','linewidth',2);
text(Sta1Lon,Sta1Lat,'\fontname{Times New Roman}Dongshan','fontsize',12,'Fontname','Times New Man','FontWeight','Bold');

axis ([116.5 119 22.5 24.5])
% set(gcf,'Position',[100 100 800 400]);
set(gca,'fontweight','Bold','fontsize',12,'linewi',1,'FontName','Times New Roman'); %调整刻度和 
xlabel('\fontname{Times New Roman}Longitude (°E)','fontsize',14,'fontweight','bold','color','k');
ylabel('\fontname{Times New Roman}Latitude (°N)','fontsize',14,'fontweight','bold','color','k')
title('\fontname{Times New Roman}Elliptical velocity at 10:29','fontsize',14,'fontweight','bold','color','k');
legend('\fontname{Times New Roman}V_E of sky-surface wave','\fontname{Times New Roman}V_E of surface wave')
% print(figure(1),'-dpng','D:\博士研究\小论文\电离层污染分析与建模\IRI_散射点定位\原始图\vb1029','-r800')
% 显示核心区域
axis ([117.5 118.2 23 23.6])
% print(figure(1),'-dpng','D:\博士研究\小论文\电离层污染分析与建模\IRI_散射点定位\原始图\vb1029_核心区域','-r800')

% 将不同距离元的径向流写入一个数据中

fr = fieldnames(w_20170404093225);
lat_v = []; % 纬度 
lon_v = []; % 经度
v_SSW = []; % 天地波流速
dir_SSW = []; % 天地波流向
dir_SW = []; % 地波流向（矢量）
v_SW = []; % 地波流速（矢量）
vr_SW = []; % 地波流速（径向）
for j = 1 : length(fr)
    data = w_20170404093225_smooth.(fr{j});
    lat_v = [lat_v,data(2,:)]; % 径向流纬度
    lon_v = [lon_v,data(3,:)]; % 径向流经度
    v_SSW = [v_SSW,data(8,:)]; % 天地波流速
    dir_SSW = [dir_SSW,data(9,:)]; % 天地波流向
    dir_SW = [dir_SW,data(11,:)]; % 地波流向
    v_SW = [v_SW,data(12,:)]; % 地波流速
    vr_SW = [vr_SW,data(13,:)];% 地波流速在径向方向上的投影
    clear data
end
% save Current_data_0944 lat_v lon_v v_SSW dir_SSW dir_SW v_SW vr_SW

%% 统一网格（极坐标）
EARTHRIDIUS = 6371; % 地球半径
len = 5:5:200; % 距离间隔
fai(1)=Sta1Lon/180*pi;       % 起点球坐标系中的方位角（弧度） 
sita(1)=pi/2-Sta1Lat/180*pi; % 终点纬度,冯：转换为球坐标系中的sita（弧度）
seaDOA_with_north = (72:0.5:237)/180*pi;
for D = 1:length(seaDOA_with_north)
% 遍历整个len和seaDOA_with_north，计算网格经纬度
    for count = 1:length(len)
        TOA(D,count) = len(count) / EARTHRIDIUS; % 终点和起点与球心之间的夹角（弧度）
        fai(2)=fai(1)+TOA(D,count)*sin(seaDOA_with_north(D))/sin(sita(1));   % 在三角锥中运用的角度近似转换，不严谨
        sita(2)=sita(1)-TOA(D,count)*cos(seaDOA_with_north(D)); % 在三角锥中运用的角度近似转换，不严谨，但误差不大
        % 沿回波方向前进len距离后海面上的经纬度
        sea_LAT(D,count) = (pi/2 - sita(2))/pi * 180;
        sea_LON(D,count) = fai(2)/pi * 180;
    end
end
sea_LAT = 23:0.002:23.6;
sea_LON = 117.5:0.002:118.2;
[sea_LON,sea_LAT] = meshgrid(sea_LON,sea_LAT);
% 采用平均滤波的方法来同一网格点
[x2,y2]=meshgrid(double(latitude_wind),double(longitude_wind));
[x3,y3]=meshgrid(double(latitude_SSH),double(longitude_SSH));
[x4,y4]=meshgrid(double(lat_h),double(lon_h));

for D = 1:length(sea_LAT(:,1))
% 遍历整个len和seaDOA_with_north，计算网格经纬度
    for count = 1:length(sea_LAT(1,:))
        % 逐点计算距离大小
        Sea_scatter_lon = sea_LON(D,count);
        Sea_scatter_lat = sea_LAT(D,count);
        [dis1,~] = distance(Sea_scatter_lon,Sea_scatter_lat,y3,x3,6371);
        % 地转
        V_f_comgrid1(D,count) = nanmean(V_f(find(dis1<5)));
        U_f_comgrid1(D,count) = nanmean(U_f(find(dis1<5)));
        [dis2,~] = distance(Sea_scatter_lon,Sea_scatter_lat,y2,x2,6371);
        % stokes
        V_VS0_comgrid1(D,count) = nanmean(V_VS0(find(dis2<26)));
        U_VS0_comgrid1(D,count) = nanmean(U_VS0(find(dis2<26)));
        V_VSZ_comgrid1(D,count) = nanmean(V_VSZ(find(dis2<26)));
        U_VSZ_comgrid1(D,count) = nanmean(U_VSZ(find(dis2<26)));
        % ekman
        Ve_comgrid1(D,count) = nanmean(Ve(find(dis2<26)));
        Ue_comgrid1(D,count) = nanmean(Ue(find(dis2<26)));
        % 风速
        U_wind_comgrid(D,count) = nanmean(U_wind(find(dis2<26)));
        V_wind_comgrid(D,count) = nanmean(V_wind(find(dis2<26)));
        % 海水深度
%         [dis3,~] = distance(Sea_scatter_lon,Sea_scatter_lat,y4,x4,6371);
%         h_comgrid(D,count) = nanmean(h(find(dis3<3)));
    end
end
% 画滤波后地转流
figure
S = load('EastCoastline.dat');
plot(S(:,1),S(:,2),'k-','linewidth',0.5);
hold on
for i = 1 :10: length(V_f_comgrid1(:,1))
    for j =  1 :10: length(V_f_comgrid1(1,:))
        if isnan(V_f_comgrid1(i,j)) || isnan(U_f_comgrid1(i,j))
            continue
        else
            quiver(sea_LON(i,j),sea_LAT(i,j),0.8*U_f_comgrid1(i,j),0.8*V_f_comgrid1(i,j),'r','AutoScale','off','MaxHeadSize',2);hold on
        end    
    end
end
axis ([115 120 21 25])
% axis ([117 118.8 22.4 24.2])
set(gcf,'Position',[100 100 800 400]);

% 画风生流ekman
figure
S = load('EastCoastline.dat');
plot(S(:,1),S(:,2),'k-','linewidth',0.5);
hold on
for i = 1 :10: length(Ve_comgrid1(:,1))
    for j =  1 :10: length(Ve_comgrid1(1,:))
        if isnan(Ve_comgrid1(i,j)) || isnan(Ue_comgrid1(i,j))
            continue
        else
            quiver(sea_LON(i,j),sea_LAT(i,j),Ue_comgrid1(i,j),Ve_comgrid1(i,j),'r','AutoScale','off','MaxHeadSize',2);hold on
        end    
    end
end
axis ([115 120 21 25])
set(gcf,'Position',[100 100 800 400]);

% 画风生流stokes
figure
S = load('EastCoastline.dat');
plot(S(:,1),S(:,2),'k-','linewidth',0.5);
hold on
for i = 1 :10: length(V_VS0_comgrid1(:,1))
    for j =  1 :10: length(V_VS0_comgrid1(1,:))
        if isnan(V_VS0_comgrid1(i,j)) || isnan(U_VS0_comgrid1(i,j))
            continue
        else
            quiver(sea_LON(i,j),sea_LAT(i,j),U_VS0_comgrid1(i,j),V_VS0_comgrid1(i,j),'r','AutoScale','off','MaxHeadSize',2);hold on
        end    
    end
end
axis ([115 120 21 25])
set(gcf,'Position',[100 100 800 400]);

% 天地波径向流 
load Current_data_0932.mat
for D = 1:length(sea_LAT(:,1))
% 遍历整个len和seaDOA_with_north，计算网格经纬度
    for count = 1:length(sea_LAT(1,:))
        Sea_scatter_lon = sea_LON(D,count);
        Sea_scatter_lat = sea_LAT(D,count);
        [dis3,~] = distance(Sea_scatter_lon,Sea_scatter_lat,lon_v,lat_v,6371);
        % 天地波径向速度
        v_SSW_comgrid(D,count) = nanmean(v_SSW(find(dis3<10)));
        % 天地波径向方向
        dir_SSW_comgrid(D,count) = nanmean(dir_SSW(find(dis3<10)));
        % 地波矢量速度
        v_SW_comgrid(D,count) = nanmean(v_SW(find(dis3<10)));
        % 地波矢量方向
        dir_SW_comgrid(D,count) = nanmean(dir_SW(find(dis3<10)));
        vr_SW_comgrid(D,count) = nanmean(vr_SW(find(dis3<10)));
    end
end
% save Comgrid_current_inform0932 wv_SSW_comgrid dir_SSW_comgrid v_SW_comgrid v_SW_comgrid V_f_comgrid1 U_f_comgrid1...
%     V_VS0_comgrid1 U_VS0_comgrid1 V_VSZ_comgrid1 U_VSZ_comgrid1 Ve_comgrid1 Ue_comgrid1 sea_LAT sea_LON
% 
% save Comgrid_current_inform1029_new v_SSW_comgrid dir_SSW_comgrid v_SW_comgrid dir_SW_comgrid V_f_comgrid1 U_f_comgrid1...
%     U_wind_comgrid V_wind_comgrid sea_LAT sea_LON h_comgrid




