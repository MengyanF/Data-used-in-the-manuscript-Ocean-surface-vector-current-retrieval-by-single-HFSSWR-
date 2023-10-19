%% ���е�վ��ز��״�ʸ�����ϳ��㷨��������
%% ���ߣ�������
%% ʱ�䣺2023.3.15
clear;
clc;
%% ���ط糡���ݺͺ���߶ȳ�����
% �糡����ע������
%   �糡����ѡ���˴�2017��3��1�յ�4��31��ÿСʱ��ECMWF����
%   ���ڲ���Ҫ���Ƿ糡��ʱ���������ҿ��ǵ� LT = UT + 7.5h
%   4��4��9��32��9��44�ķ糡������ V10(:,:,819) �� U10(:,:,819)
%   10��19��10��29���� V10(:,:,820) �� U(:,:,820)
% ����߶�ע������
%   ����߶ȵ�ʱ��ֱ���ֻ��ÿ�գ�����̨�庣Ͽ�������һ������
%   ����10cm/s�ĵ�ת������˸�ʱ��ֱ��ʵĺ���߶������㹻

% load wind_data.mat
load wind_data_4_4.mat
load SSH.mat
% ��ˮ���
lat_h = ncread('D:\��ʿ�о�\С����\ʸ�����ϳ�\program\����߶�\gebco_2022_n26.0_s20.0_w115.0_e122.0.nc','lat');
lon_h = ncread('D:\��ʿ�о�\С����\ʸ�����ϳ�\program\����߶�\gebco_2022_n26.0_s20.0_w115.0_e122.0.nc','lon');
h = ncread('D:\��ʿ�о�\С����\ʸ�����ϳ�\program\����߶�\gebco_2022_n26.0_s20.0_w115.0_e122.0.nc','elevation');

%% ���ݺ���߶ȼ����ת��

SSH = (SSH(:,:,4) + SSH(:,:,3))./2;
[x,y] = size(SSH);
% �ϱ�����
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
% ��������
for i =1 : x
    for j = 1 : y-1
        f = double(7.272e-5 * 2 * sind(latitude_SSH(j)));
        U_f(i,j) = -1 * (SSH(i,j+1)- SSH(i,j))./ double(((latitude_SSH(j+1)-latitude_SSH(j))*100000)) * 9.8 ./ f;
    end
        U_f(i,y) = 0;
end
% ��ͼ
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

% ��ת���� 
cd_speed = sqrt(U_f.^2 + V_f.^2)*100;
figure
pcolor(longitude_SSH,latitude_SSH,cd_speed');hold on
shading interp
hc=colorbar;
set(gca,'Clim',[0 60]);  % ����colorbar�ķ�Χ
S = load('EastCoastline.dat');
plot(S(:,1),S(:,2),'k-','linewidth',0.5);
hold on
% axis ([115 120 21 25])
axis ([117.5 118.2 23 23.6]) % ��������
set(gcf,'Position',[100 100 800 400]);
ylabel(hc,'\fontname{Times New Roman}Current speed (cm/s)','fontweight','bold','fontsize',14);
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %�����̶Ⱥ� 
xlabel('\fontname{Times New Roman}Longitude (��E)','fontsize',14,'fontweight','bold','color','k');
ylabel('\fontname{Times New Roman}Latitude (��N)','fontsize',14,'fontweight','bold','color','k')
title('\fontname{Times New Roman}2017/04/04','fontsize',14,'fontweight','bold','color','k')




%  ����߶Ȼ�ͼ
figure
% ��ɽ��γ��
Sta1Lon=117.4863;Sta1Lat=23.6575; 
% �����γ��
SiteLatChS = 24.03702;SiteLonChS = 117.9025;
S = load('EastCoastline.dat');
plot(S(:,1),S(:,2),'k-','linewidth',0.5);
hold on
pcolor(longitude_SSH,latitude_SSH,SSH');hold on
shading interp
hc=colorbar;
axis ([115 120 21 25]) % ��������
colormap(jet)
set(gca,'Clim',[0.3 0.9]);  % ����colorbar�ķ�Χ
plot(Sta1Lon,Sta1Lat,'p','linewidth',2);
text(Sta1Lon,Sta1Lat,'\fontname{Times New Roman}Dongshan','fontsize',12,'Fontname','Times New Man','FontWeight','Bold');
plot(SiteLonChS,SiteLatChS,'p','linewidth',2);
text(SiteLonChS,SiteLatChS,'\fontname{Times New Roman}Chihu','fontsize',12,'Fontname','Times New Man','FontWeight','Bold');
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %�����̶Ⱥ� 
xlabel('\fontname{Times New Roman}Longitude (��E)','fontsize',14,'fontweight','bold','color','k');
ylabel('\fontname{Times New Roman}Latitude (��N)','fontsize',14,'fontweight','bold','color','k')
title('\fontname{Times New Roman}2017/04/04 \fontname{����}����߶�','fontsize',14,'fontweight','bold','color','k')
ylabel(hc,'\fontname{Times New Roman}SSH (m)','fontweight','bold','fontsize',14);



%% ������
% U_wind = U10(:,:,819);
% V_wind = V10(:,:,819);
%% Ekman��
pa = 1.29; % �����ܶ�
pw = 1.05*1e3; % ��ˮ�ܶ�
% ˮ��
h = 20;
% ����ճ��,��λ m^2/s
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

% Ekman����
figure
pcolor(longitude_wind,latitude_wind,abs(VE)');hold on
shading interp
hc=colorbar;
% set(gca,'Clim',[0 1]);  % ����colorbar�ķ�Χ
S = load('EastCoastline.dat');
plot(S(:,1),S(:,2),'k-','linewidth',0.5);
hold on
axis ([115 120 21 25]);
set(gcf,'Position',[100 100 800 400]);
ylabel(hc,'\fontname{Times New Roman}Ekman speed (m)','fontweight','bold','fontsize',14);
set(gca,'fontweight','Bold','fontsize',10,'linewi',1); %�����̶Ⱥ� 
xlabel('\fontname{Times New Roman}Longitude (��E)','fontsize',14,'fontweight','bold','color','k');
ylabel('\fontname{Times New Roman}Latitude (��N)','fontsize',14,'fontweight','bold','color','k');
title('\fontname{Times New Roman}2017/04/04','fontsize',14,'fontweight','bold','color','k');

% Ekman ʸ��
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

%% Stokes��
% stokes���ƻ����ò���������
% JONSWAP��
F =100e3;% ������С
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

% Stokes ʸ��
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


%% ��ز��״ﾶ����
% load('D:\��ʿ�о�\С����\ʸ�����ϳ�\program\������\scatter_point_information_9_32_new_eta_0.5853_GP_GIA_all_information.mat')
% load('D:\��ʿ�о�\С����\ʸ�����ϳ�\program\������\scatter_point_information_9_44_new_eta_0.5609_GP_GIA_all_information.mat')
% load('D:\��ʿ�о�\С����\ʸ�����ϳ�\program\������\scatter_point_information_10_19_new_eta_0.4456_GP_GIA_all_information.mat')
load('D:\��ʿ�о�\С����\ʸ�����ϳ�\program\������\scatter_point_information_10_29_new_eta_0.3995_GP_GIA_all_information.mat')

% ��ɽ��γ��
Sta1Lon=117.4863;Sta1Lat=23.6575; 
% �����γ��
SiteLatChS = 24.03702;SiteLonChS = 117.9025;
% ���γ��
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
            % ��ɫ�Ƿ���ֵ����ɫ�ǹ۲�ֵ
            quiver(Sea_scatter_lon(i),Sea_scatter_lat(i),scale*abs(vb_IRI(i)).*sind(dir),scale*abs(vb_IRI(i)).*cosd(dir),'r','AutoScale','off','MaxHeadSize',2);hold on
            % �ز�ʸ������ dir_IRI �ϵ�ͶӰ��dir_IRI��v_sur_Dir֮������ǣ�v_sur����0 
            v_sur = v_sur_speed(i) * cosd(dir_IRI(i) - v_sur_Dir(i));
            if v_sur > 0 % ˵���ز��ľ���������Ϊ dir_IRI(i)
                quiver(Sea_scatter_lon(i),Sea_scatter_lat(i),scale*abs(v_sur).*sind(dir_IRI(i)),scale*abs(v_sur).*cosd(dir_IRI(i)),'b','AutoScale','off','MaxHeadSize',2);
                hold on
            elseif v_sur < 0 % ˵���ز��ľ���������Ϊ dir
                quiver(Sea_scatter_lon(i),Sea_scatter_lat(i),scale*abs(v_sur).*sind(dir),scale*abs(v_sur).*cosd(dir),'b','AutoScale','off','MaxHeadSize',2);
                hold on
            end
        elseif vb_IRI(i) > 0 && vb_IRI(i) < 50 && v_sur_speed(i)<50
            dir = dir_IRI(i);
            quiver(Sea_scatter_lon(i),Sea_scatter_lat(i),scale*abs(vb_IRI(i)).*sind(dir),scale*abs(vb_IRI(i)).*cosd(dir),'r','AutoScale','off','MaxHeadSize',2);hold on
            v_sur = v_sur_speed(i) * cosd(dir_IRI(i) - v_sur_Dir(i));
            if v_sur > 0  % ˵���ز��ľ���������Ϊ dir
                quiver(Sea_scatter_lon(i),Sea_scatter_lat(i),scale*abs(v_sur).*sind(dir),scale*abs(v_sur).*cosd(dir),'b','AutoScale','off','MaxHeadSize',2);
                hold on
            elseif v_sur < 0  % ˵���ز��ľ���������Ϊ dir+180
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
set(gca,'fontweight','Bold','fontsize',12,'linewi',1,'FontName','Times New Roman'); %�����̶Ⱥ� 
xlabel('\fontname{Times New Roman}Longitude (��E)','fontsize',14,'fontweight','bold','color','k');
ylabel('\fontname{Times New Roman}Latitude (��N)','fontsize',14,'fontweight','bold','color','k')
title('\fontname{Times New Roman}Elliptical velocity at 10:29','fontsize',14,'fontweight','bold','color','k');
legend('\fontname{Times New Roman}V_E of sky-surface wave','\fontname{Times New Roman}V_E of surface wave')
% print(figure(1),'-dpng','D:\��ʿ�о�\С����\�������Ⱦ�����뽨ģ\IRI_ɢ��㶨λ\ԭʼͼ\vb1029','-r800')
% ��ʾ��������
axis ([117.5 118.2 23 23.6])
% print(figure(1),'-dpng','D:\��ʿ�о�\С����\�������Ⱦ�����뽨ģ\IRI_ɢ��㶨λ\ԭʼͼ\vb1029_��������','-r800')

% ����ͬ����Ԫ�ľ�����д��һ��������

fr = fieldnames(w_20170404093225);
lat_v = []; % γ�� 
lon_v = []; % ����
v_SSW = []; % ��ز�����
dir_SSW = []; % ��ز�����
dir_SW = []; % �ز�����ʸ����
v_SW = []; % �ز����٣�ʸ����
vr_SW = []; % �ز����٣�����
for j = 1 : length(fr)
    data = w_20170404093225_smooth.(fr{j});
    lat_v = [lat_v,data(2,:)]; % ������γ��
    lon_v = [lon_v,data(3,:)]; % ����������
    v_SSW = [v_SSW,data(8,:)]; % ��ز�����
    dir_SSW = [dir_SSW,data(9,:)]; % ��ز�����
    dir_SW = [dir_SW,data(11,:)]; % �ز�����
    v_SW = [v_SW,data(12,:)]; % �ز�����
    vr_SW = [vr_SW,data(13,:)];% �ز������ھ������ϵ�ͶӰ
    clear data
end
% save Current_data_0944 lat_v lon_v v_SSW dir_SSW dir_SW v_SW vr_SW

%% ͳһ���񣨼����꣩
EARTHRIDIUS = 6371; % ����뾶
len = 5:5:200; % ������
fai(1)=Sta1Lon/180*pi;       % ���������ϵ�еķ�λ�ǣ����ȣ� 
sita(1)=pi/2-Sta1Lat/180*pi; % �յ�γ��,�룺ת��Ϊ������ϵ�е�sita�����ȣ�
seaDOA_with_north = (72:0.5:237)/180*pi;
for D = 1:length(seaDOA_with_north)
% ��������len��seaDOA_with_north����������γ��
    for count = 1:length(len)
        TOA(D,count) = len(count) / EARTHRIDIUS; % �յ�����������֮��ļнǣ����ȣ�
        fai(2)=fai(1)+TOA(D,count)*sin(seaDOA_with_north(D))/sin(sita(1));   % ������׶�����õĽǶȽ���ת�������Ͻ�
        sita(2)=sita(1)-TOA(D,count)*cos(seaDOA_with_north(D)); % ������׶�����õĽǶȽ���ת�������Ͻ���������
        % �ػز�����ǰ��len��������ϵľ�γ��
        sea_LAT(D,count) = (pi/2 - sita(2))/pi * 180;
        sea_LON(D,count) = fai(2)/pi * 180;
    end
end
sea_LAT = 23:0.002:23.6;
sea_LON = 117.5:0.002:118.2;
[sea_LON,sea_LAT] = meshgrid(sea_LON,sea_LAT);
% ����ƽ���˲��ķ�����ͬһ�����
[x2,y2]=meshgrid(double(latitude_wind),double(longitude_wind));
[x3,y3]=meshgrid(double(latitude_SSH),double(longitude_SSH));
[x4,y4]=meshgrid(double(lat_h),double(lon_h));

for D = 1:length(sea_LAT(:,1))
% ��������len��seaDOA_with_north����������γ��
    for count = 1:length(sea_LAT(1,:))
        % ����������С
        Sea_scatter_lon = sea_LON(D,count);
        Sea_scatter_lat = sea_LAT(D,count);
        [dis1,~] = distance(Sea_scatter_lon,Sea_scatter_lat,y3,x3,6371);
        % ��ת
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
        % ����
        U_wind_comgrid(D,count) = nanmean(U_wind(find(dis2<26)));
        V_wind_comgrid(D,count) = nanmean(V_wind(find(dis2<26)));
        % ��ˮ���
%         [dis3,~] = distance(Sea_scatter_lon,Sea_scatter_lat,y4,x4,6371);
%         h_comgrid(D,count) = nanmean(h(find(dis3<3)));
    end
end
% ���˲����ת��
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

% ��������ekman
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

% ��������stokes
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

% ��ز������� 
load Current_data_0932.mat
for D = 1:length(sea_LAT(:,1))
% ��������len��seaDOA_with_north����������γ��
    for count = 1:length(sea_LAT(1,:))
        Sea_scatter_lon = sea_LON(D,count);
        Sea_scatter_lat = sea_LAT(D,count);
        [dis3,~] = distance(Sea_scatter_lon,Sea_scatter_lat,lon_v,lat_v,6371);
        % ��ز������ٶ�
        v_SSW_comgrid(D,count) = nanmean(v_SSW(find(dis3<10)));
        % ��ز�������
        dir_SSW_comgrid(D,count) = nanmean(dir_SSW(find(dis3<10)));
        % �ز�ʸ���ٶ�
        v_SW_comgrid(D,count) = nanmean(v_SW(find(dis3<10)));
        % �ز�ʸ������
        dir_SW_comgrid(D,count) = nanmean(dir_SW(find(dis3<10)));
        vr_SW_comgrid(D,count) = nanmean(vr_SW(find(dis3<10)));
    end
end
% save Comgrid_current_inform0932 wv_SSW_comgrid dir_SSW_comgrid v_SW_comgrid v_SW_comgrid V_f_comgrid1 U_f_comgrid1...
%     V_VS0_comgrid1 U_VS0_comgrid1 V_VSZ_comgrid1 U_VSZ_comgrid1 Ve_comgrid1 Ue_comgrid1 sea_LAT sea_LON
% 
% save Comgrid_current_inform1029_new v_SSW_comgrid dir_SSW_comgrid v_SW_comgrid dir_SW_comgrid V_f_comgrid1 U_f_comgrid1...
%     U_wind_comgrid V_wind_comgrid sea_LAT sea_LON h_comgrid




