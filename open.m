function [ S ] = open( ~ )
%open()�����������������Ҫ���ڶ�ȡ��Ƶ�ļ�
%�������룬�Լ��������Ƶ�ļ�����һ������S��
%SΪ�ú��������

[y1,fs1]= audioread('wash5.wav');   %��ȡ��Ƶ�ļ���yΪ��Ϣ����fsΪ����Ƶ��
[y2,fs2]= audioread('dragen5.wav');
[y3,fs3]= audioread('music5.wav');

S1=y1(:,1); % ��ȡ�� 1 ����
S2=y2(:,1); % ��ȡ�� 1 ����
S3=y3(:,1); % ��ȡ�� 1 ����
% audiowrite('1234.WAV',S1,fs1); % ʵ�� 1 ��������
% audiowrite('2234.WAV',S2,fs2); % ʵ�� 1 ��������
% audiowrite('3234.WAV',S3,fs3); % ʵ�� 1 ��������

I1=S1';                          %��һ������Ϣ��Ϊ���������Ա�����ϳ�S����
I2=S2';
I3=S3';

[~,aIcol]=size(I1);              %S����ÿһ�д���һ����Ƶ����Ϣ
[~,bIcol]=size(I2);
[~,cIcol]=size(I3);
dI=[aIcol,bIcol,cIcol];
dI=max(dI);                      %�����Ƶ��Ϣ��������

S=zeros(3,dI);                   %��������Ƶ�źŰ�������S�У���Ϣ����ʱ����
S(1,1:aIcol)=I1;
S(2,1:bIcol)=I2;
S(3,1:cIcol)=I3;

Si1=S(1,1:aIcol);                
Si2=S(2,1:bIcol);
Si3=S(3,1:cIcol);

subplot(331);plot(Si1),title('�����ź�1');
subplot(332);plot(Si2),title('�����ź�2');%axis([0,1000,-2,2]);
subplot(333);plot(Si3),title('�����ź�3');%axis([0,1000,-2,2]);

% figure;plot(Si1),title('�����ź�1');axis([0,50000,-0.4,0.4]);
% figure;plot(Si2),title('�����ź�2');%axis([0,1000,-2,2]);
% figure;plot(Si3),title('�����ź�3');%axis([0,1000,-2,2]);



end

