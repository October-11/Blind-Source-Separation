clc;
[y,Fs]=audioread('wash10.wav');
%Y2=y(((Fs*0+1):Fs*10),:);
%Y2=Y2(:,1);
y=y(:,1);
                             
dt=1/Fs                         %采样间隔时间1/fs
n=length(y)                     %得到序列的长度
t=[0:n-1]*dt

a=fft(y)
b=[0:n-1]/(n*dt);               %频率序列
c=abs(a);                       %傅里叶变换后的振幅


figure(1)
subplot(3,1,1)
plot(t,y)
xlabel('Time')
ylabel('Audio Signal')
title('原始序列')
subplot(3,1,2)
plot(b,c*2/n)                   %绘制信号的振幅谱
xlabel('频率/HZ'),title('幅频图')
ylabel('振幅')
subplot(3,1,3)
plot(b(1:n/2),c(1:n/2));        %绘制nyquist频率之前随频率变化的振幅
xlabel('频率/HZ');
ylabel('振幅');
title('Nyquist频率')
