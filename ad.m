clc;
[y,Fs]=audioread('wash10.wav');
%Y2=y(((Fs*0+1):Fs*10),:);
%Y2=Y2(:,1);
y=y(:,1);
                             
dt=1/Fs                         %�������ʱ��1/fs
n=length(y)                     %�õ����еĳ���
t=[0:n-1]*dt

a=fft(y)
b=[0:n-1]/(n*dt);               %Ƶ������
c=abs(a);                       %����Ҷ�任������


figure(1)
subplot(3,1,1)
plot(t,y)
xlabel('Time')
ylabel('Audio Signal')
title('ԭʼ����')
subplot(3,1,2)
plot(b,c*2/n)                   %�����źŵ������
xlabel('Ƶ��/HZ'),title('��Ƶͼ')
ylabel('���')
subplot(3,1,3)
plot(b(1:n/2),c(1:n/2));        %����nyquistƵ��֮ǰ��Ƶ�ʱ仯�����
xlabel('Ƶ��/HZ');
ylabel('���');
title('NyquistƵ��')
