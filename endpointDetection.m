close all;
clear;
clc;

[data,Fs]=audioread('3.mp3');
M=size(data,1);
y = data;
x = zscore(y);
figure;
subplot(211);plot(y);
subplot(212);plot(x);
%% 计算短时能量
N=length(x);
%%短时窗长, 窗长选用10, 即帧长
L=10;
%帧移
k=4;
%帧数
num=ceil(N/k);
%求窗函数
for i=1:L
   %w(i)=0.54-0.46*cos(2*pi*(i-1)/(L-1));%汉明窗
   %w(i)=0.5*(1-cos(2*pi*(i-1)/(L-1)))
   w(i)=1; % 矩形窗
end
% 计算每个窗的短时能量，并记录在E中
for i=1:num 
    sm=0;
    k1=(i-1)*k;
    for j=1:L
        if(k1+j)>N 
            break;
        end
        sm = ( x(k1+j)*w(j) )^2 + sm;
    end
    E(i)=sm;
end

%% 计算过零率
Z = zeros(1,num);
for i = 1:num 
    sm = 0;
    k1 = (i-1)*k;
    for j = 1:L
        if (k1+j) > N 
            break;
        end
        if (k1+j-1)==0 
            continue;
        end
        sm = abs( sign(x(k1+j)) - sign(x(k1+j-1)) ) * w(j) + sm;
    end
    Z(i) = sm;
end

%% 设置门限
key = max(E)/10;
ZeroLevel = 0.2;  % 设置短时过零率门限,0.3是自已通过试验设置的
% 计算短时能量和短时过零率的门限
A = max(E);
B = min(E);
EL1 = ( min(E) * ( 1 + 2*log10(max(E)/min(E)) ) )*10;
js = 0;
zh = 0;
for i=1:num 
    if (E(i) > EL1)
        js = js + 1;
        zh = zh + E(i);
    end
end
SL = zh / js;
EL2 = EL1 / 20 + 0.25*(SL - EL1);
ZL = max(Z) * ZeroLevel;

%% 根据短时能量判断起止点
% 根据较高的门限初步判断
voicedIndex = find( E >= EL2 );
sound=[];
kk=1;
if( length(voicedIndex) ~= 0 )
    sound(kk).begin = voicedIndex(1);
    for i = 2: length(voicedIndex)-1
        if voicedIndex(i+1) - voicedIndex(i) > 1
            sound(kk).end = voicedIndex(i);
            sound(kk+1).begin = voicedIndex(i+1);
            kk=kk+1;
        end
    end
    sound(kk).end=voicedIndex(end);
    index=[];
% 忽略一些小的细节
    for i=1:length(sound)
        if sound(i).end - sound(i).begin < 3
            index = [index,i];
        end
    end
    sound(index) = [];
end

% 根据较低的门限进一步判断（将根据较高门限判断出的波段前后的小波幅震动包含进来）
for i=1:length(sound)
    head = sound(i).begin;
    while (head-1) >= 1 && E(head-1) >= EL1
        head = head-1;
    end
    sound(i).begin = head;
    tail=sound(i).end;
    while (tail+1) <= length(E) && E(tail+1) > EL1
        tail = tail+1;
    end
    sound(i).end = tail;
end

%% 根据过零率门限判断起止点
for i=1:length(sound)
    head = sound(i).begin;
    while (head-1) >= 1 && Z(head-1) >= ZL
        head = head-1;
    end
    sound(i).begin = head;
    tail = sound(i).end;
    while (tail+1) <= length(Z) && Z(tail+1) > ZL;
        tail = tail+1;
    end
    sound(i).end = tail;
end

% 去掉重复的语音帧
if (length(sound) ~= 0)
    index=[];
    for i=1:length(sound)-1
        if sound(i).begin == sound(i+1).begin && sound(i).end == sound(i+1).end
            index=[index,i];
        end
    end
    sound(index)=[];
end

%% 将序列变换成整段语音的样点序列
if length(sound) ~= 0
    for i = 1:length(sound)
        out(i).begin = (sound(i).begin - 1) * (L-k) + 1;
        out(i).end = (sound(i).end) * (L-k) + k;
    end
else
    out=[];
end

%端点检测图
%时域波形
figure(1);
subplot(311);
plot(x,'LineWidth',2);
title('幅度图 红线:始点 绿线:止点');
xlabel('时间 t');
ylabel('幅度 A');
y=[-1 1];
for i=1:length(sound);
    line(sound(i).begin * k * [1,1], y, 'color', 'r');
    line(sound(i).end * k * [1,1] + 1, y, 'color', 'g');
end

%基于短时能量的端点检测图
subplot(312);
%能量图
[t1, t2] = size(E);
kk = linspace(1, t2, t2);
subplot(312);
plot(kk, E);
title('能量图');
xlabel('帧序列 n');
ylabel('能量 E');
y=[min(E) max(E)];
for i=1:length(sound);
    line(sound(i).begin*[1,1], y, 'color', 'r');
    line(sound(i).end*[1,1], y, 'color', 'g');
end

%基于短时过零率的端点检测图
subplot(313);
%过零率图
[t1,t2] = size(Z);
km = linspace(1, t2, t2);
plot(km, Z);
title('过零率图');
xlabel('帧序列 n');
ylabel('过零率 Z');
y=[0 max(Z)];
for i=1:length(sound);
    line(sound(i).begin*[1,1], y, 'color', 'r');
    line(sound(i).end*[1,1], y, 'color', 'g');
end