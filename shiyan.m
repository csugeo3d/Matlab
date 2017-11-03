clc,clear,close all;
%��˸���������
point=[
  3.7, 1.7, 2.1; 
  4.1, 3.8, 3.2;
  4.7, 2.9, 1.9;
  5.2, 2.8, 5.8;
  6.0, 4.0, 7.7;
  6.3, 3.6, 3.9;
  9.7, 6.3, 5.1;
  10.0, 4.9, 4.6;
  11.0, 3.6, 9.2;
  12.5, 6.4, 8.8
];

%ԭʼ�����
x = point(:,1);
y = point(:,2);
z = point(:,3);


%�������ĵ�
minX = min(x);
maxX = max(x);
minY = min(y);
maxY = max(y);
minZ = min(z);
maxZ = max(z);

centerX = (minX+maxX)./2;
centerY = (minY+maxY)./2;
centerZ = (minZ+maxZ)./2;

%���Ʒ�Χ
dx = (maxX-minX).*1/4;
dy = (maxY-minY).*1/4;


%����Э������
CONV_xyz=cov(point);
%eig����������ֵ����������,����u������������������;
%u�Ѿ��ǵ�λ��������,����������u(:,1)'*u(:,2)=0,u�����µ������᷽��
%v������ֵ�������Խ�Ԫ�ؾ�������ֵ�������������͵��ж�Ӧ
[u,v] = eig(CONV_xyz);
% dd = max(point(:,1))-min(point(:,1));
% xx = linspace(centerX-dd,centerX+dd,50);  %������ʾ��Χ
% yy = linspace(centerX-dd,centerX+dd,50); 

t = [-10 10];

%��ת�����
%x����ת������ᣬ�ߵĵ㷨ʽ����
x_x=u(1).*t+centerX;
x_y=u(2).*t+centerY;
x_z=u(3).*t+centerZ;
%y����ת������ᣬ�ߵĵ㷨ʽ����
y_x=u(4).*t+centerX;
y_y=u(5).*t+centerY;
y_z=u(6).*t+centerZ;
%z����ת������ᣬ�ߵĵ㷨ʽ����
z_x=u(7).*t+centerX;
z_y=u(8).*t+centerY;
z_z=u(9).*t+centerZ;


[len,wid] = size(point);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %ѭ�����ƴ�ֱ�����ƽ�棬���ֱ����ͶӰ��
% %�õ�u(1),u(2),u(3)˵��Ϊx���ϵ�ͶӰ��
for i=1:len
    A(i,:) = point(i,:);  %����ȡ��ÿһ�У�x,y,z��
    xx(i) = A(i,1);    %��i�еĵ�1�У���xx�д�����е�x
    yy(i) = A(i,2);    %��i�еĵ�2�У���yy�д�����е�y
    zz(i) = A(i,3);    %��i�еĵ�3�У���zz�д�����е�z
    
    %�������е��ƽ�棨Ϊ�˹۲췽�㣩
%     xArea = (minX-dx):0.25:(maxX+dx);
%     yArea = (minY-dy):0.25:(maxY+dy);
%     [X,Y] = meshgrid(xArea,yArea);
%     Z = -((u(1).*(X-xx(i))+u(2).*(Y-yy(i)))./u(3))+zz(i);
%     mesh(X,Y,Z);
%     hold on
    
    %����tt,���ͶӰ������
    tx(i) = -(u(1).*(centerX-xx(i))+u(2).*(centerY-yy(i))+u(3).*(centerZ-zz(i)))/(u(1).^2+u(2).^2+u(3).^2);
    xpro_x(i)=u(1).*tx(i)+centerX;
    xpro_y(i)=u(2).*tx(i)+centerY;
    xpro_z(i)=u(3).*tx(i)+centerZ;
    
    %����ͶӰ������
%     plot3(xpro_x(i),xpro_y(i),xpro_z(i),'xg');
%     hold on
%     
%     axis equal
%     grid on

end

%�ڻ�ȡ��ͶӰ����һ����ת�ᣨ��x��ת����ᣩ�ϵ������֮��
%ȡ��ͶӰ������x����Сx�������㣬�����Χ���ص�һ����İ���
[ma_prox,I_prox]=max(xpro_x);  %��������ֵ�����ֵ��λ��
pro_x_maxPoint = [ma_prox,xpro_y(I_prox),xpro_z(I_prox)];
[mi_prox,J_prox]=min(xpro_x); 
pro_x_minPoint = [mi_prox,xpro_y(J_prox),xpro_z(J_prox)];
dis_pro_x = (norm(pro_x_maxPoint-pro_x_minPoint))./2;  %������ת���ϵİ볤��

%��ʾ�����ĵ�����ת���x�ᴹֱ��ƽ��
xArea = (minX-dx):0.25:(maxX+dx);
yArea = (minY-dy):0.25:(maxY+dy);
[X,Y] = meshgrid(xArea,yArea);
Z = -((u(1).*(X-centerX)+u(2).*(Y-centerY))./u(3))+centerZ;
% mesh(X,Y,Z);
% hold on

%����ƽ��ֱ��ط�������������ƽ�ư볤������1/4���ܹ�����������1/4
x_xMAX = X + 5/4.*dis_pro_x .* u(1);
y_xMAX = Y + 5/4.*dis_pro_x .* u(2);
z_xMAX = Z + 5/4.*dis_pro_x .* u(3);

x_xMIN = X - 5/4.*dis_pro_x .* u(1);
y_xMIN = Y - 5/4.*dis_pro_x .* u(2);
z_xMIN = Z - 5/4.*dis_pro_x .* u(3);








% mesh(x_xMAX,y_xMAX,z_xMAX);
% alpha(.2)
% % hidden off
% hold on
% mesh(x_xMIN,y_xMIN,z_xMIN);
% alpha(.2)
% % hidden off
% hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ѭ�����ƴ�ֱ�����ƽ�棬���ֱ����ͶӰ��
%�õ�u(4),u(5),u(6)˵��Ϊx���ϵ�ͶӰ��
for i=1:len
    A(i,:) = point(i,:);  %����ȡ��ÿһ�У�x,y,z��
    xx(i) = A(i,1);    %��i�еĵ�1�У���xx�д�����е�x
    yy(i) = A(i,2);    %��i�еĵ�2�У���yy�д�����е�y
    zz(i) = A(i,3);    %��i�еĵ�3�У���zz�д�����е�z
    
    %�������е��ƽ�棨Ϊ�˹۲췽�㣩
%     xArea = minX:0.25:maxX;
%     yArea = minY:0.25:maxY;
%     [X,Y] = meshgrid(xArea,yArea);
%     Z = -((u(1).*(X-xx(i))+u(2).*(Y-yy(i)))./u(3))+zz(i);
%     mesh(X,Y,Z);
%     hold on
    
    %����tt,���ͶӰ������
    ty(i) = -(u(4).*(centerX-xx(i))+u(5).*(centerY-yy(i))+u(6).*(centerZ-zz(i)))/(u(4).^2+u(5).^2+u(6).^2);
    ypro_x(i)=u(4).*ty(i)+centerX;
    ypro_y(i)=u(5).*ty(i)+centerY;
    ypro_z(i)=u(6).*ty(i)+centerZ;
    
    %����ͶӰ������
%     plot3(ypro_x(i),ypro_y(i),ypro_z(i),'*g');
%     hold on
%     
%     axis equal
%     grid on

end

%�ڻ�ȡ��ͶӰ����һ����ת�ᣨ��x��ת����ᣩ�ϵ������֮��
%ȡ��ͶӰ������x����Сx�������㣬�����Χ���ص�һ����İ���
[ma_proy,I_proy]=max(ypro_y);  %��������ֵ�����ֵ��λ��
pro_y_maxPoint = [ypro_x(I_proy),ma_proy,ypro_z(I_proy)];
[mi_proy,J_proy]=min(ypro_y); 
pro_y_minPoint = [ypro_x(J_proy),mi_proy,ypro_z(J_proy)];
dis_pro_y = (norm(pro_y_maxPoint-pro_y_minPoint))./2;  %������ת���ϵİ볤��

%��ʾ�����ĵ�����ת���x�ᴹֱ��ƽ��
xArea = (minX-dx):0.25:(maxX+dx);
yArea = (minY-dy):0.25:(maxY+dy);
[X,Y] = meshgrid(xArea,yArea);
Z = -((u(4).*(X-centerX)+u(5).*(Y-centerY))./u(6))+centerZ;
% mesh(X,Y,Z);
% hold on

%����ƽ��ֱ��ط�������������ƽ�ư볤�����
x_yMAX = X + 5/4.*dis_pro_y .* u(4);
y_yMAX = Y + 5/4.*dis_pro_y .* u(5);
z_yMAX = Z + 5/4.*dis_pro_y .* u(6);

x_yMIN = X - 5/4.*dis_pro_y .* u(4);
y_yMIN = Y - 5/4.*dis_pro_y .* u(5);
z_yMIN = Z - 5/4.*dis_pro_y .* u(6);

% mesh(x_yMAX,y_yMAX,z_yMAX);
% alpha(.2)
% % hidden off
% hold on
% alpha(.2)
% % hidden off
% mesh(x_yMIN,y_yMIN,z_yMIN);
% hold on


%��x���������y�����Ľ���


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ѭ�����ƴ�ֱ�����ƽ�棬���ֱ����ͶӰ��
%�õ�u(7),u(8),u(9)˵��Ϊx���ϵ�ͶӰ��
for i=1:len
    A(i,:) = point(i,:);  %����ȡ��ÿһ�У�x,y,z��
    xx(i) = A(i,1);    %��i�еĵ�1�У���xx�д�����е�x
    yy(i) = A(i,2);    %��i�еĵ�2�У���yy�д�����е�y
    zz(i) = A(i,3);    %��i�еĵ�3�У���zz�д�����е�z
    
    %�������е��ƽ�棨Ϊ�˹۲췽�㣩
%     xArea = minX:0.25:maxX;
%     yArea = minY:0.25:maxY;
%     [X,Y] = meshgrid(xArea,yArea);
%     Z = -((u(1).*(X-xx(i))+u(2).*(Y-yy(i)))./u(3))+zz(i);
%     mesh(X,Y,Z);
%     hold on
    
    %����tt,���ͶӰ������
    tz(i) = -(u(7).*(centerX-xx(i))+u(8).*(centerY-yy(i))+u(9).*(centerZ-zz(i)))/(u(7).^2+u(8).^2+u(9).^2);
    zpro_x(i)=u(7).*tz(i)+centerX;
    zpro_y(i)=u(8).*tz(i)+centerY;
    zpro_z(i)=u(9).*tz(i)+centerZ;
    
    %����ͶӰ������
%     plot3(zpro_x(i),zpro_y(i),zpro_z(i),'xg');
%     hold on
%     
%     axis equal
%     grid on

end

%�ڻ�ȡ��ͶӰ����һ����ת�ᣨ��x��ת����ᣩ�ϵ������֮��
%ȡ��ͶӰ������x����Сx�������㣬�����Χ���ص�һ����İ���
[ma_proz,I_proz]=max(zpro_z);  %��������ֵ�����ֵ��λ��
pro_z_maxPoint = [zpro_x(I_proz),zpro_y(I_proz),ma_proz];
[mi_proz,J_proz]=min(zpro_z);
pro_z_minPoint = [zpro_x(J_proz),zpro_y(J_proz),mi_proz];
dis_pro_z = (norm(pro_z_maxPoint-pro_z_minPoint))./2;  %������ת���ϵİ볤��

%��ʾ�����ĵ�����ת���x�ᴹֱ��ƽ��
xArea = (minX-dx):0.25:(maxX+dx);
yArea = (minY-dy):0.25:(maxY+dy);
[X,Y] = meshgrid(xArea,yArea);
Z = -((u(7).*(X-centerX)+u(8).*(Y-centerY))./u(9))+centerZ;
% mesh(X,Y,Z);
% hold on

%����ƽ��ֱ��ط�������������ƽ�ư볤�����
x_zMAX = X + 5/4.*dis_pro_z .* u(7);
y_zMAX = Y + 5/4.*dis_pro_z .* u(8);
z_zMAX = Z + 5/4.*dis_pro_z .* u(9);

x_zMIN = X - 5/4.*dis_pro_z .* u(7);
y_zMIN = Y - 5/4.*dis_pro_z .* u(8);
z_zMIN = Z - 5/4.*dis_pro_z .* u(9);

% mesh(x_zMAX,y_zMAX,z_zMAX);
% alpha(.2)
% % hidden off
% hold on
% alpha(.2)
% % hidden off
% mesh(x_zMIN,y_zMIN,z_zMIN);
% hold on

%��ƽ��Ľ���x�����ϵ���С�棬�ֱ���y�������С��z�������С�涼�ཻ

%ƽ���ᣬ�����ת���x�ᣬ�ֱ�����ת���y,z��ֱ�ƽ�ƶ�Ӧ�İ��᳤

%x����ת������ᣬ�ߵĵ㷨ʽ����
%��ת���x����y�᷽������ƽ�ƺ�����������ϵı�

% xy_xMAX = x_x + 5/4.*dis_pro_y .* u(4);
% xy_yMAX = x_y + 5/4.*dis_pro_y .* u(5);
% xy_zMAX = x_z + 5/4.*dis_pro_y .* u(6);
% 
% xy_xMIN = x_x - 5/4.*dis_pro_y .* u(4);
% xy_yMIN = x_y - 5/4.*dis_pro_y .* u(5);
% xy_zMIN = x_z - 5/4.*dis_pro_y .* u(6);
% 
% %����ת���z�᷽���Ϸֱ�����ƽ��
% xyz_xmaxT = xy_xMAX + 5/4.*dis_pro_z .* u(7);
% xyz_ymaxT = xy_yMAX + 5/4.*dis_pro_z .* u(8);
% xyz_zmaxT = xy_zMAX + 5/4.*dis_pro_z .* u(9);
% 
% xyz_xmaxL = xy_xMAX - 5/4.*dis_pro_z .* u(7);
% xyz_ymaxL = xy_yMAX - 5/4.*dis_pro_z .* u(8);
% xyz_zmaxL = xy_zMAX - 5/4.*dis_pro_z .* u(9);
% 
% xyz_xminT = xy_xMIN + 5/4.*dis_pro_z .* u(7);
% xyz_yminT = xy_yMIN + 5/4.*dis_pro_z .* u(8);
% xyz_zminT = xy_zMIN + 5/4.*dis_pro_z .* u(9);
% 
% xyz_xminL = xy_xMIN - 5/4.*dis_pro_z .* u(7);
% xyz_yminL = xy_yMIN - 5/4.*dis_pro_z .* u(8);
% xyz_zminL = xy_zMIN - 5/4.*dis_pro_z .* u(9);

%�����ϼ������̣��ֱ����һ��ֻ���в���t1�Ĳ���������ʽ

xyz_xmaxT = u(1).*t+centerX + 5/4.*dis_pro_y .* u(4) + 5/4.*dis_pro_z .* u(7);
xyz_ymaxT = u(2).*t+centerY + 5/4.*dis_pro_y .* u(5) + 5/4.*dis_pro_z .* u(8);
xyz_zmaxT = u(3).*t+centerZ + 5/4.*dis_pro_y .* u(6) + 5/4.*dis_pro_z .* u(9);

xyz_xmaxL = u(1).*t+centerX + 5/4.*dis_pro_y .* u(4) - 5/4.*dis_pro_z .* u(7);
xyz_ymaxL = u(2).*t+centerY + 5/4.*dis_pro_y .* u(5) - 5/4.*dis_pro_z .* u(8);
xyz_zmaxL = u(3).*t+centerZ + 5/4.*dis_pro_y .* u(6) - 5/4.*dis_pro_z .* u(9);
% 
% xyz_xminT = xy_xMIN + 5/4.*dis_pro_z .* u(7);
% xyz_yminT = xy_yMIN + 5/4.*dis_pro_z .* u(8);
% xyz_zminT = xy_zMIN + 5/4.*dis_pro_z .* u(9);
% 
% xyz_xminL = xy_xMIN - 5/4.*dis_pro_z .* u(7);
% xyz_yminL = xy_yMIN - 5/4.*dis_pro_z .* u(8);
% xyz_zminL = xy_zMIN - 5/4.*dis_pro_z .* u(9);



%��ת���y����x����������ƽ�ƺ�ı�
yx_xMAX = y_x + 5/4.*dis_pro_x .* u(1);
yx_yMAX = y_y + 5/4.*dis_pro_x .* u(2);
yx_zMAX = y_z + 5/4.*dis_pro_x .* u(3);

yx_xMIN = y_x - 5/4.*dis_pro_x .* u(1);
yx_yMIN = y_y - 5/4.*dis_pro_x .* u(2);
yx_zMIN = y_z - 5/4.*dis_pro_x .* u(3);

%����ת���z�᷽���Ϸֱ�����ƽ��
% yxz_xmaxT = yx_xMAX + 5/4.*dis_pro_z .* u(7);
% yxz_ymaxT = yx_yMAX + 5/4.*dis_pro_z .* u(8);
% yxz_zmaxT = yx_zMAX + 5/4.*dis_pro_z .* u(9);
% 
% yxz_xmaxL = yx_xMAX - 5/4.*dis_pro_z .* u(7);
% yxz_ymaxL = yx_yMAX - 5/4.*dis_pro_z .* u(8);
% yxz_zmaxL = yx_zMAX - 5/4.*dis_pro_z .* u(9);
% 
% yxz_xminT = yx_xMIN + 5/4.*dis_pro_z .* u(7);
% yxz_yminT = yx_yMIN + 5/4.*dis_pro_z .* u(8);
% yxz_zminT = yx_zMIN + 5/4.*dis_pro_z .* u(9);
% 
% yxz_xminL = yx_xMIN - 5/4.*dis_pro_z .* u(7);
% yxz_yminL = yx_yMIN - 5/4.*dis_pro_z .* u(8);
% yxz_zminL = yx_zMIN - 5/4.*dis_pro_z .* u(9);

%������ת���y,�����ϼ������̣��ֱ����һ��ֻ���в���t1�Ĳ���������ʽ
yxz_xmaxT = u(4).*t +centerX + 5/4.*dis_pro_x .* u(1) + 5/4.*dis_pro_z .* u(7);
yxz_ymaxT = u(5).*t+centerY + 5/4.*dis_pro_x .* u(2) + 5/4.*dis_pro_z .* u(8);
yxz_zmaxT = u(6).*t+centerZ + 5/4.*dis_pro_x .* u(3) + 5/4.*dis_pro_z .* u(9);

yxz_xmaxL = u(4).*t+centerX + 5/4.*dis_pro_x .* u(1) - 5/4.*dis_pro_z .* u(7);
yxz_ymaxL = u(5).*t+centerY + 5/4.*dis_pro_x .* u(2) - 5/4.*dis_pro_z .* u(8);
yxz_zmaxL = u(6).*t+centerZ + 5/4.*dis_pro_x .* u(3) - 5/4.*dis_pro_z .* u(9);


%��x+T��y+T�Ľ�������(x0,y0,z0)
t = 5/4.*(dis_pro_x.*u(1)-dis_pro_y.*u(4))./(u(1)-u(4));
x0 = u(1).*t+centerX + 5/4.*dis_pro_y .* u(4) + 5/4.*dis_pro_z .* u(7);
% x1 = u(1).*t+centerX + 5/4.*dis_pro_y .* u(4) - 5/4.*dis_pro_z .* u(7);
t = 5/4.*(dis_pro_x.*u(2)-dis_pro_y.*u(5))./(u(2)-u(5));
y0 = u(2).*t+centerY + 5/4.*dis_pro_y .* u(5) + 5/4.*dis_pro_z .* u(8);
% y1 = u(2).*t+centerY + 5/4.*dis_pro_y .* u(5) - 5/4.*dis_pro_z .* u(8);
t = 5/4.*(dis_pro_x.*u(3)-dis_pro_y.*u(6))./(u(3)-u(6));
z0 = u(3).*t+centerZ + 5/4.*dis_pro_y .* u(6) + 5/4.*dis_pro_z .* u(9);
% z1 = u(6).*t+centerZ + 5/4.*dis_pro_x .* u(3) - 5/4.*dis_pro_z .* u(9);
%��x+L��y+L�Ľ�������(x1,y1,z1)
% t = 5/4.*(dis_pro_x.*u(1)-dis_pro_y.*u(4))./(u(1)-u(4));
% x1 = u(1).*t+centerX + 5/4.*dis_pro_y .* u(4) - 5/4.*dis_pro_z .* u(7);



%��ת���z��ƽ�ƺ�����������ϵı�
% zx_xMAX = z_x + 5/4.*dis_pro_x .* u(1);
% zx_yMAX = z_y + 5/4.*dis_pro_x .* u(2);
% zx_zMAX = z_z + 5/4.*dis_pro_x .* u(3);
% 
% zx_xMIN = z_x - 5/4.*dis_pro_x .* u(1);
% zx_yMIN = z_y - 5/4.*dis_pro_x .* u(2);
% zx_zMIN = z_z - 5/4.*dis_pro_x .* u(3);
%����ת���y�᷽���Ϸֱ�����ƽ��
% zxy_xmaxT = zx_xMAX + 5/4.*dis_pro_y .* u(4);
% zxy_ymaxT = zx_yMAX + 5/4.*dis_pro_y .* u(5);
% zxy_zmaxT = zx_zMAX + 5/4.*dis_pro_y .* u(6);
% 
% zxy_xmaxL = zx_xMAX - 5/4.*dis_pro_y .* u(4);
% zxy_ymaxL = zx_yMAX - 5/4.*dis_pro_y .* u(5);
% zxy_zmaxL = zx_zMAX - 5/4.*dis_pro_y .* u(6);
% 
% zxy_xminT = zx_xMIN + 5/4.*dis_pro_y .* u(4);
% zxy_yminT = zx_yMIN + 5/4.*dis_pro_y .* u(5);
% zxy_zminT = zx_zMIN + 5/4.*dis_pro_y .* u(6);
% 
% zxy_xminL = zx_xMIN - 5/4.*dis_pro_y .* u(4);
% zxy_yminL = zx_yMIN - 5/4.*dis_pro_y .* u(5);
% zxy_zminL = zx_zMIN - 5/4.*dis_pro_y .* u(6);

%���彻������
% pointBox = [
%     x0,y0,z0;
%     x1,y1,z1          
% ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����չʾ����

%���Ƽ�������
plot3(x0,y0,z0,'*black');
hold on


%������ת��xƽ�ƺ�ı�
plot3(xyz_xmaxT,xyz_ymaxT,xyz_zmaxT,'g');
hold on
% plot3(xyz_xmaxL,xyz_ymaxL,xyz_zmaxL,'g');
% hold on
% plot3(xyz_xminT,xyz_yminT,xyz_zminT,'g');
% hold on
% plot3(xyz_xminL,xyz_yminL,xyz_zminL,'g');
% hold on

%������ת��yƽ�ƺ�ı�
plot3(yxz_xmaxT,yxz_ymaxT,yxz_zmaxT,'r');
hold on
% plot3(yxz_xmaxL,yxz_ymaxL,yxz_zmaxL,'r');
% hold on
% plot3(yxz_xminT,yxz_yminT,yxz_zminT,'r');
% hold on
% plot3(yxz_xminL,yxz_yminL,yxz_zminL,'r');
% hold on

%������ת��zƽ�ƺ�ı�
% plot3(zxy_xmaxT,zxy_ymaxT,zxy_zmaxT,'b');
% hold on
% 
% plot3(zxy_xmaxL,zxy_ymaxL,zxy_zmaxL,'b');
% hold on
% plot3(zxy_xminT,zxy_yminT,zxy_zminT,'b');
% hold on
% 
% plot3(zxy_xminL,zxy_yminL,zxy_zminL,'b');
% hold on

%����������
plot3(x_x,x_y,x_z,'b');
hold on
plot3(y_x,y_y,y_z,'b');
hold on
plot3(z_x,z_y,z_z,'b');
hold on


%���ĵ�����
plot3(centerX,centerY,centerZ,'*b');
hold on
%ԭʼ������
plot3(x,y,z,'or');
hold on

xlabel('x');
ylabel('y');
zlabel('z');
axis equal
grid on

