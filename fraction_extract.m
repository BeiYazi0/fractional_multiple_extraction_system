close all;clear all;clc;
%% 输入信号和相关参数
fs=15e6;                         % 采样率
n=512;
f=@(t) sin(fs*pi/2*t)./t/pi/fs;  % 输入信号
x=[0.5,f((1:n-1)/fs)];           % 采样信号
L=3;M=5;                         % 0.6倍抽取
wp=0.15;                         % 通带截止频率(关于pi归一化）
wst=0.25;                        % 阻带截止频率(关于pi归一化）

%% L倍内插
x1=zeros(1,n*L);
x1(1:L:end)=x;                  

%% 低通滤波
h=myfilter(wp,wst,2);
hplot(h)
x2=filter(h,1,x1);

%% M倍抽取
x2=L*x2;
y=x2(1:M:end);

%% 功能函数

% 方案选择
function h=myfilter(wp,wst,index)
if index==1
    h=hammingwin(wp,wst);
elseif index==2
    h=freqsample(wp,wst);
else
    h=freqsample_plus(wp,wst);
end
end

% 窗函数设计低通FIR滤波器
function h=hammingwin(wp,wst)
% wp:通带截止频率(关于pi归一化）
% wst:阻带截止频率(关于pi归一化）
% h:返回设计好的FIR滤波器系数
wc=(wp+wst)/2;                % 理想低通滤波器截止频率(关于pi归一化）
wm=wst-wp;                    % 过渡带带宽(关于pi归一化）
N=8/wm*0.825;                 % 计算汉明窗的阶数
N=fix(N)+1;
n=0:N;
w=0.54-0.46*cos(2*pi*n/N);    %窗函数
alpha=N/2;
p=pi*(n-alpha);
hd=sin(wc*p)./p;              % 理想的频率响应(截取段)
if fix(alpha)==alpha
    hd(alpha+1)=wc;
end
h=hd.*w;
end

% 频率采样法设计低通FIR滤波器
function h=freqsample(wp,wst)
% wp:通带截止频率(关于pi归一化）
% wst:阻带截止频率(关于pi归一化）
% h:返回设计好的FIR滤波器系数
wc=(wp+wst)/2;                % 理想低通滤波器截止频率(关于pi归一化）
N=182;                         % 阶数+1
kc=fix(wc*N/2)+1;             
alpha=(N-1)/2;
Hrs=ones(1,N);             
Hrs(kc:(N+1-kc))=0;           % Hg(k)
k1=0:floor((N-1)/2);
k2=floor((N-1)/2)+1:N-1; 
angH=[-alpha*(2*pi)/N*k1,alpha*(2*pi)/N*(N-k2)]; % θ(k)
H=Hrs.*exp(1i*angH);          % H(k)
h=real(ifft(H,N)); 
end

% 频率采样法设计低通FIR滤波器(增加过渡采样点)
function h=freqsample_plus(wp,wst)
% wp:通带截止频率(关于pi归一化）
% wst:阻带截止频率(关于pi归一化）
% h:返回设计好的FIR滤波器系数
wc=(wp+wst)/2;                % 理想低通滤波器截止频率(关于pi归一化）
N=48;                         % 阶数+1
kc=fix(wc*N/2)+1;             
alpha=(N-1)/2;
Hrs=ones(1,N);             
Hrs(kc+1:(N+1-kc))=0;           % Hg(k)
Hrs(kc+1)=0.5925;Hrs(N+1-kc)=0.5925; % 增加过渡采样点;
Hrs(kc+2)=0.1099;Hrs(N-kc)=0.1099;
k1=0:floor((N-1)/2);
k2=floor((N-1)/2)+1:N-1; 
angH=[-alpha*(2*pi)/N*k1,alpha*(2*pi)/N*(N-k2)]; % θ(k)
H=Hrs.*exp(1i*angH);          % H(k)
h=real(ifft(H,N)); 
end

% 绘制滤波器的脉冲响应和频率响应
function hplot(h)
n=size(h,2)-1;
m=0:n;
figure(1);grid on;
subplot(211);
plot(m,h);	                         % 绘制脉冲响应
title('脉冲响应');
axis([0 n 1.1*min(h) 1.1*max(h)]);
ylabel('h(n)');xlabel('n');
m=0:511; f=m/511;
H=abs(fft(h,1024));                  % 获取幅频响应
subplot(212);
plot(f,H(1:512));                    % 绘制幅频响应                 
title('幅频响应');
axis([0 1 1.1*min(H) 1.1*max(H)]);
ylabel('幅值');xlabel('频率(关于pi归一化）');
figure(2);grid on;
freqz(h,1);                          % 绘制幅频特性和相频特性
subplot(211);
title('幅频特性');
subplot(212);
title('相频特性');
end

