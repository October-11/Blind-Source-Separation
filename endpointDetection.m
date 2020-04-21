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
%% �����ʱ����
N=length(x);
%%��ʱ����, ����ѡ��10, ��֡��
L=10;
%֡��
k=4;
%֡��
num=ceil(N/k);
%�󴰺���
for i=1:L
   %w(i)=0.54-0.46*cos(2*pi*(i-1)/(L-1));%������
   %w(i)=0.5*(1-cos(2*pi*(i-1)/(L-1)))
   w(i)=1; % ���δ�
end
% ����ÿ�����Ķ�ʱ����������¼��E��
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

%% ���������
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

%% ��������
key = max(E)/10;
ZeroLevel = 0.2;  % ���ö�ʱ����������,0.3������ͨ���������õ�
% �����ʱ�����Ͷ�ʱ�����ʵ�����
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

%% ���ݶ�ʱ�����ж���ֹ��
% ���ݽϸߵ����޳����ж�
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
% ����һЩС��ϸ��
    for i=1:length(sound)
        if sound(i).end - sound(i).begin < 3
            index = [index,i];
        end
    end
    sound(index) = [];
end

% ���ݽϵ͵����޽�һ���жϣ������ݽϸ������жϳ��Ĳ���ǰ���С�����𶯰���������
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

%% ���ݹ����������ж���ֹ��
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

% ȥ���ظ�������֡
if (length(sound) ~= 0)
    index=[];
    for i=1:length(sound)-1
        if sound(i).begin == sound(i+1).begin && sound(i).end == sound(i+1).end
            index=[index,i];
        end
    end
    sound(index)=[];
end

%% �����б任��������������������
if length(sound) ~= 0
    for i = 1:length(sound)
        out(i).begin = (sound(i).begin - 1) * (L-k) + 1;
        out(i).end = (sound(i).end) * (L-k) + k;
    end
else
    out=[];
end

%�˵���ͼ
%ʱ����
figure(1);
subplot(311);
plot(x,'LineWidth',2);
title('����ͼ ����:ʼ�� ����:ֹ��');
xlabel('ʱ�� t');
ylabel('���� A');
y=[-1 1];
for i=1:length(sound);
    line(sound(i).begin * k * [1,1], y, 'color', 'r');
    line(sound(i).end * k * [1,1] + 1, y, 'color', 'g');
end

%���ڶ�ʱ�����Ķ˵���ͼ
subplot(312);
%����ͼ
[t1, t2] = size(E);
kk = linspace(1, t2, t2);
subplot(312);
plot(kk, E);
title('����ͼ');
xlabel('֡���� n');
ylabel('���� E');
y=[min(E) max(E)];
for i=1:length(sound);
    line(sound(i).begin*[1,1], y, 'color', 'r');
    line(sound(i).end*[1,1], y, 'color', 'g');
end

%���ڶ�ʱ�����ʵĶ˵���ͼ
subplot(313);
%������ͼ
[t1,t2] = size(Z);
km = linspace(1, t2, t2);
plot(km, Z);
title('������ͼ');
xlabel('֡���� n');
ylabel('������ Z');
y=[0 max(Z)];
for i=1:length(sound);
    line(sound(i).begin*[1,1], y, 'color', 'r');
    line(sound(i).end*[1,1], y, 'color', 'g');
end