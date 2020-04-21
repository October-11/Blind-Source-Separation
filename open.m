function [ S ] = open( ~ )
%open()函数无输入参数，主要用于读取音频文件
%声道分离，以及将多个音频文件放入一个矩阵S中
%S为该函数的输出

[y1,fs1]= audioread('wash5.wav');   %读取音频文件，y为信息矩阵，fs为抽样频率
[y2,fs2]= audioread('dragen5.wav');
[y3,fs3]= audioread('music5.wav');

S1=y1(:,1); % 抽取第 1 声道
S2=y2(:,1); % 抽取第 1 声道
S3=y3(:,1); % 抽取第 1 声道
% audiowrite('1234.WAV',S1,fs1); % 实现 1 声道分离
% audiowrite('2234.WAV',S2,fs2); % 实现 1 声道分离
% audiowrite('3234.WAV',S3,fs3); % 实现 1 声道分离

I1=S1';                          %将一声道信息变为行向量，以便于组合成S矩阵
I2=S2';
I3=S3';

[~,aIcol]=size(I1);              %S矩阵每一行代表一个音频的信息
[~,bIcol]=size(I2);
[~,cIcol]=size(I3);
dI=[aIcol,bIcol,cIcol];
dI=max(dI);                      %求得音频信息最大的列数

S=zeros(3,dI);                   %将三个音频信号按行收入S中，信息不够时补零
S(1,1:aIcol)=I1;
S(2,1:bIcol)=I2;
S(3,1:cIcol)=I3;

Si1=S(1,1:aIcol);                
Si2=S(2,1:bIcol);
Si3=S(3,1:cIcol);

subplot(331);plot(Si1),title('输入信号1');
subplot(332);plot(Si2),title('输入信号2');%axis([0,1000,-2,2]);
subplot(333);plot(Si3),title('输入信号3');%axis([0,1000,-2,2]);

% figure;plot(Si1),title('输入信号1');axis([0,50000,-0.4,0.4]);
% figure;plot(Si2),title('输入信号2');%axis([0,1000,-2,2]);
% figure;plot(Si3),title('输入信号3');%axis([0,1000,-2,2]);



end

