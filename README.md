# Matlab
repository
function main20170326_1
clc;
clear;
close all;
warning off;
%Step 1: Load image 
[filename, pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';...
          '*.*','All Files' });
RGB = imread([ pathname,filename]);
imshow(RGB);
title('原始图像');

I=imnoise(RGB,'salt',0.01);
figure
imshow(I);
title('灰尘噪声干扰');
% Store (X,Y) offsets for later use; subtract 1 so that each offset will
%Step 3: Threshold the image
I = rgb2gray(I);
threshold = graythresh(I);
BW = im2bw(I,threshold);
BW = ~BW;  
figure
imshow(BW);
title('二值化图像');
BW = medfilt2(BW,[7 7]);  %对二值图像进行中值滤波，去掉细小的干扰点
figure
imshow(BW)
title('去干扰后的图像')
se=strel('rectangle',[5,5]);
BW=imclose(BW,se);
figure
imshow(BW);
title('闭运算后的图像'); %闭运算通常用来填充目标内细小空洞，连接断开的邻近目标，平滑其边界的同时不明显改变其面积

se=strel('disk',5);%这里是创建一个半径为5的平坦型圆盘结构元素,参数可调
BW=imopen(BW,se);%直接开运算
figure
imshow(BW);
title('开运算后的图像');%开运算属于形态学图像处理，是先腐蚀后膨胀，作用是：可以使边界平滑，消除细小的尖刺，断开窄小的连接,保持面积大小不变等。

%找到左边界和右边边界
[L,N]=size(BW);
Array=[];
for i=1:L
    guodu=BW(i,:);
    index=find(guodu==0);
    if length(index>0)
        Array=[Array;[index(1),index(end),i]];
    end
end
%%找出上边界和下边界
Brray=[];
for j=1:N
    guodu=BW(:,j);
    index=find(guodu==0);
    if length(index>0)
        Brray=[Brray;[index(1),index(end),j]];
    end
end
%%做出拟合曲线
%%左边界
%%消除差别过大的点
data=Array(:,1);
xscale=Array(:,3);
[tc,xcenter]=hist(data);
figure
hist(data);
Pro=cumsum(tc)/sum(tc);
index=find(Pro>0.8);
Dc=xcenter(index(1));%
Deta_Vatlue=abs(data-Dc);
De=(xcenter(2)-xcenter(1));
index=find(Deta_Vatlue<=De);
Value_Left=data(index);
Left_Scale=xscale(index);
figure
plot(Left_Scale,Value_Left,'ro');
hold on;
plot(xscale,data);
%%拟合
p = polyfit(Left_Scale,Value_Left,1);
low=min(xscale);
up=max(xscale);
N=1000;
x=linspace(low,up,N);
y = polyval(p,x);
plot(x,y,'-k');
title('左边界');
legend('剔除后的边界点','原始边界点','拟合直线');
Crray{1}=[x;y];
Drray{1}=p;
%%
%%右边界
%%消除差别过大的点
data=Array(:,2);
xscale=Array(:,3);
[tc,xcenter]=hist(data);
figure
hist(data);
Pro=cumsum(tc)/sum(tc);
index=find(Pro>0.8);
Dc=xcenter(index(1));%
Deta_Vatlue=abs(data-Dc);
De=(xcenter(2)-xcenter(1));
index=find(Deta_Vatlue<=De);
Value_Left=data(index);
Left_Scale=xscale(index);
figure
plot(Left_Scale,Value_Left,'ro');
hold on;
plot(xscale,data);
%%拟合
p = polyfit(Left_Scale,Value_Left,1);
low=min(xscale);
up=max(xscale);
N=1000;
x=linspace(low,up,N);
y = polyval(p,x);
plot(x,y,'-k');
title('右边界');
legend('剔除后的边界点','原始边界点','拟合直线');
Crray{2}=[x;y];
Drray{2}=p;
%%上边界
%%消除差别过大的点
data=Brray(:,1);
xscale=Brray(:,3);
[tc,xcenter]=hist(data);
figure
hist(data);
Pro=cumsum(tc)/sum(tc);
index=find(Pro>0.8);
Dc=xcenter(index(1));%
Deta_Vatlue=abs(data-Dc);
De=(xcenter(2)-xcenter(1));
index=find(Deta_Vatlue<=De);
Value_Left=data(index);
Left_Scale=xscale(index);
figure
plot(Left_Scale,Value_Left,'ro');
hold on;
plot(xscale,data);
%%拟合
p = polyfit(Left_Scale,Value_Left,1); %除去偏离的点
low=min(xscale);
up=max(xscale);
N=1000;
x=linspace(low,up,N);
y = polyval(p,x);
plot(x,y,'-k');
title('上边界');
legend('剔除后的边界点','原始边界点','拟合直线');
Crray{3}=[x;y];
Drray{3}=p;
%%
%%下边界
%%消除差别过大的点
data=Brray(:,2);
xscale=Brray(:,3);
[tc,xcenter]=hist(data);
figure
hist(data);
Pro=cumsum(tc)/sum(tc);
index=find(Pro>0.8);
Dc=xcenter(index(1));%
Deta_Vatlue=abs(data-Dc);
De=(xcenter(2)-xcenter(1));
index=find(Deta_Vatlue<=De);
Value_Left=data(index);
Left_Scale=xscale(index);
figure
plot(Left_Scale,Value_Left,'ro');
hold on;
plot(xscale,data);
%%拟合
p = polyfit(Left_Scale,Value_Left,1);
low=min(xscale);
up=max(xscale);
N=1000;
x=linspace(low,up,N);
y = polyval(p,x);
plot(x,y,'-k');
title('下边界');
legend('剔除后的边界点','原始边界点','拟合直线');
Crray{4}=[x;y];
Drray{4}=p;
figure
imshow(I);
hold on;
guodu=Crray{1};
plot(guodu(2,:),guodu(1,:),'-r','linewidth',2);
guodu=Crray{2};
plot(guodu(2,:),guodu(1,:),'-r','linewidth',2);
guodu=Crray{3};
plot(guodu(1,:),guodu(2,:),'-r','linewidth',2);
guodu=Crray{4};
plot(guodu(1,:),guodu(2,:),'-r','linewidth',2);
parameter=Drray{1};
disp(['左边界直线方程为：y=',num2str(parameter(1)),'*x+',num2str(parameter(2))]);
parameter=Drray{2};
disp(['右边界直线方程为：y=',num2str(parameter(1)),'*x+',num2str(parameter(2))]);
parameter=Drray{3};
disp(['上边界直线方程为：y=',num2str(parameter(1)),'*x+',num2str(parameter(2))]);
parameter=Drray{4};
disp(['下边界直线方程为：y=',num2str(parameter(1)),'*x+',num2str(parameter(2))]);
%% 给出准确的纸币中心位置坐标模型以及长宽模型
%%左上角
guodu1=Crray{1};
guodu3=Crray{3};
%%取均值
X=mean([guodu1(2,1),guodu3(1,1)]);
Y=mean([guodu1(1,1),guodu3(2,1)]);
Left_Up=[X,Y];
plot(Left_Up(1),Left_Up(2),'g^');
%%左下角
guodu1=Crray{1};
guodu4=Crray{4};
%%取均值
X=mean([guodu1(2,end),guodu4(1,1)]);
Y=mean([guodu1(1,end),guodu4(2,1)]);
Left_Down=[X,Y];
plot(Left_Down(1),Left_Down(2),'g^');
%%右上角
guodu2=Crray{2};
guodu3=Crray{3};
%%取均值
X=mean([guodu2(2,1),guodu3(1,end)]);
Y=mean([guodu2(1,1),guodu3(2,end)]);
Right_Up=[X,Y];
plot(Right_Up(1),Right_Up(2),'g^');
%%右下角
guodu2=Crray{2};
guodu4=Crray{4};
%%取均值
X=mean([guodu2(2,end),guodu4(1,end)]);
Y=mean([guodu2(1,end),guodu4(2,end)]);
Right_Down=[X,Y];
plot(Right_Down(1),Right_Down(2),'g^');
%%
Center=mean([Left_Up;Left_Down;Right_Up;Right_Down],1);
disp('坐标中心');
disp(Center);
L=mean([norm(Right_Up-Left_Up),norm(Right_Down-Left_Up)]);
N=mean([norm(Left_Up-Left_Down),norm(Right_Up-Right_Down)]);
disp(['纸币长度为 ',num2str(max(L,N))]);
disp(['纸币宽度为 ',num2str(min(L,N))]);
