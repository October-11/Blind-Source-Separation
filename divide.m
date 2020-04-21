
%%%%%%%%%%%%%%%%%%%%%%%%%%  初始化  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
clc;clear all;close all;

 
%%%%%%%%%%%%%%  读入原始音频文件，混合，并输出混合音频  %%%%%%%%%%%%%%%%%% 

[S]=open();

d=size(S);
Sweight=randn(size(S,1));             % 取一随机矩阵，作为信号混合的权矩阵 
MixedS=Sweight*S;                      % 得到三个混合信号矩阵 
 
% 将混合矩阵重新排列并输出 
subplot(3,3,4),plot(MixedS(1,:)),title('混合信号1'), axis([0,50000,-0.5,0.5]);
subplot(3,3,5),plot(MixedS(2,:)),title('混合信号2'), axis([0,50000,-0.5,0.5]);
subplot(3,3,6),plot(MixedS(3,:)),title('混合信号3'), axis([0,50000,-0.5,0.5]);

% figure,plot(MixedS(1,:)),title('混合信号1'),% axis([0,1000,-4,4]);
% figure,plot(MixedS(2,:)),title('混合信号2'), %axis([0,1000,-4,4]);
% figure,plot(MixedS(3,:)),title('混合信号3'), %axis([0,1000,-4,4]);

fs=44100;
audiowrite('123.WAV',MixedS(1,:),fs);  
audiowrite('223.WAV',MixedS(2,:),fs); 
audiowrite('323.WAV',MixedS(3,:),fs);

 
MixedS_bak=MixedS;                      % 将混合后的数据备份，以便在恢复时直接调用 


%%%%%%%%%%%%%%%%%%%%%%%%%%  标准化  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
MixedS_mean=zeros(3,1); 
for i=1:3 
    MixedS_mean(i)=mean(MixedS(i,:)); 
end                                        % 计算MixedS的均值 
 
for i=1:3 
    for j=1:size(MixedS,2) 
        MixedS(i,j)=MixedS(i,j)-MixedS_mean(i);    %去均值
    end 
end                                        
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%  白化  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
MixedS_cov=cov(MixedS');                    % cov为求协方差的函数 
[E,D]=eig(MixedS_cov);                      % 对图片矩阵的协方差函数进行特征值分解 E特征向量，D特征值
Q=inv(sqrt(D))*(E)';                        % Q为白化矩阵 
MixedS_white=Q*MixedS;                      % MixedS_white为白化后的音频信息矩阵 
IsI=cov(MixedS_white');                     % IsI应为单位阵  
% figure;
% subplot(131);plot(MixedS_white(1,:)),title('白化输出信号1'),  %axis([0,1000,-2,2]);
% subplot(132);plot(MixedS_white(2,:)),title('白化输出信号2'),   axis([0,50000,-5,5]);
% subplot(133);plot(MixedS_white(3,:)),title('白化输出信号3'),  %axis([0,1000,-2,2]);
 
%%%%%%%%%%%%%%%%%%%%%%%%　FASTICA算法  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
X=MixedS_white;                            % 以下算法将对X进行操作 
[VariableNum,SampleNum]=size(X); 
numofIC=VariableNum;                       % 在此应用中，独立元个数等于变量个数 
B=zeros(numofIC,VariableNum) ;            % 初始化列向量w的寄存矩阵,B=[b1  b2  ...   bd] 
for r=1:numofIC                            % 迭代求取每一个独立元 
    i=1;
    maxIterationsNum=100;                  % 设置最大迭代次数（即对于每个独立分量而言迭代均不超过此次数） 
   % IterationsNum=0; 
    b=rand(numofIC,1);                     % 随机设置b初值 
   
    b=b/norm(b);                           % 对b标准化 
    while i<=maxIterationsNum+1 
        if i == maxIterationsNum           % 循环结束处理 
            fprintf('\n第%d分量在%d次迭代内并不收敛。', r,maxIterationsNum); 
            break; 
        end 
        bOld=b;                        
        t=X'*b; 
        g=(exp(2.*t)-1)./(exp(2.*t)+1); 
        dg=4*exp(2.*t)./(exp(2.*t)+1).^2;   %导数
        b=(X*g)/SampleNum-mean(dg)*b; 
                                            % 核心公式
        b=b-B*B'*b;                         % 对b正交化 
        b=b/norm(b); 
        
        if abs(abs(b'*bOld)-1)<1e-9         % 如果收敛，则保存b 
               B(:,r)=b; 
               break; 
        end 
        i=i+1;         
     end 
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%  数据复原  %%%%%%%%%%%%%%%%%%%%%%%% 
 
ICAedS=B\X;                    % x(t)=B*s(t);s(t)=B\x(t); (左除)
 
% 将混合矩阵重新排列并输出 
subplot(3,3,7),plot(ICAedS(1,:)),title('ICA输出信号1'),  %axis([0,1000,-2,2]);
subplot(3,3,8),plot(ICAedS(2,:)),title('ICA输出信号2'),   axis([0,50000,-5,5]);
subplot(3,3,9),plot(ICAedS(3,:)),title('ICA输出信号3'),  %axis([0,1000,-2,2]);

% figure,plot(ICAedS(1,:)),title('ICA输出信号1'),  %axis([0,1000,-2,2]);
% figure,plot(ICAedS(2,:)),title('ICA输出信号2'),   axis([0,50000,-5,5]);
% figure,plot(ICAedS(3,:)),title('ICA输出信号3'),  %axis([0,1000,-2,2]);

%fs=44100;
audiowrite('12345.WAV',ICAedS(1,:),fs);  
audiowrite('22345.WAV',ICAedS(2,:),fs); 
audiowrite('32345.WAV',ICAedS(3,:),fs);
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%  PCA计算并构图  %%%%%%%%%%%%%%%%%%%%%%%% 