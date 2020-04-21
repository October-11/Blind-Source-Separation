clc;
[y,Fs]=audioread('wash5.wav');
%Y2=y(((Fs*0+1):Fs*10),:);
%Y2=Y2(:,1);


dt=1/Fs                         %采样间隔时间1/fs
n=length(y)                     %得到序列的长度
t=[0:n-1]*dt

figure(1)
plot(y)
xlabel('Time')
ylabel('Audio Signal')
title('原始序列')
   