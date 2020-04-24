
%%%%%%%%%%%%%%%%%%%%%%%%%%  初始化  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all;close all;

%%%%%%%%%%%%%%  读入原始图像，混合，并输出混合图像  %%%%%%%%%%%%%%%%%%

% 读入混合前的原始图片并显示
I1=audioread ('washmix3.wav')';
figure
plot(I1),title('输入声音1'),

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  进行EMD分解  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[imf] = emd(I1,'fix',10)
F=1:length(I1);
[G H]=size(imf);

figure
for i=1:G
    subplot(G,1,i);
    plot(F,imf(i,:));
    ylabel (['IMF ' num2str(i)]);
    set(gca,'xtick',[])
    xlim([1 length(I1)])
end;

%figure(2)
%subplot(G,1,G)
%plot(F,imf(G,:))
%xlabel('Sampling Site') ;
%ylabel(['IMF ' num2str(G)]);
%xlim([1 length(I1)])

% 将EMD分解后的分量组成矩阵
S=[imf(1,:);imf(2,:);imf(3,:);imf(4,:);imf(5,:);imf(6,:);imf(7,:);imf(8,:);imf(9,:);imf(10,:);imf(11,:);imf(12,:);imf(13,:);imf(14,:)]; 
                                       % 因此S_all是一个变量个数＊采样个数的矩阵
Sweight=rand(size(S,1));               % 取一随机矩阵，作为信号混合的权矩阵
MixedS=Sweight*S;                      % 得到三个图像的混合信号矩阵
MixedS_bak=MixedS;                     % 将混合后的数据备份，以便在恢复时直接调用
%%%%%%%%%%%%%%%%%%%%%%%%%%  标准化  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MixedS_mean=zeros(14,1);
for i=1:14
    MixedS_mean(i)=mean(MixedS(i,:));
end                                        % 计算MixedS的均值

for i=1:14
    for j=1:size(MixedS,2)
        MixedS(i,j)=MixedS(i,j)-MixedS_mean(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%  白化  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MixedS_cov=cov(MixedS');                    % cov为求协方差的函数
[E,D]=eig(MixedS_cov);                      % 对图片矩阵的协方差函数进行特征值分解
Q=inv(sqrt(D))*(E)';                        % Q为白化矩阵
MixedS_white=Q*MixedS;                      % MixedS_white为白化后的图片矩阵
IsI=cov(MixedS_white');                     % IsI应为单位阵            

%%%%%%%%%%%%%%%%%%%%%%%%　FASTICA算法  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=MixedS_white;                            % 以下算法将对X进行操作
[VariableNum,SampleNum]=size(X);
numofIC=VariableNum;                       % 在此应用中，独立元个数等于变量个数
B=zeros(numofIC,VariableNum);              % 初始化列向量w的寄存矩阵,B=[b1  b2  ...   bd]
for r=1:numofIC
    i=1;maxIterationsNum=100;               % 设置最大迭代次数（即对于每个独立分量而言迭代均不超过此次数）
    IterationsNum=0;
    b=rand(numofIC,1)-.5;                  % 随机设置b初值
    b=b/norm(b);                           % 对b标准化 norm(b):向量元素平方和开根号
    while i<=maxIterationsNum+1
        if i == maxIterationsNum           % 循环结束处理
            fprintf('\n第%d分量在%d次迭代内并不收敛。', r,maxIterationsNum);
            break;
        end
        bOld=b;                          
        a2=1;
        u=1;
        t=X'*b;
        g=t.*exp(-a2*t.^2/2);
        dg=(1-a2*t.^2).*exp(-a2*t.^2/2);
        b=((1-u)*t'*g*b+u*X*g)/SampleNum-mean(dg)*b;
                                           % 核心公式，参见理论部分公式2.52
        b=b-B*B'*b;                        % 对b正交化
        b=b/norm(b); 
        if abs(abs(b'*bOld)-1)<1e-9        % 如果收敛，则
             B(:,r)=b;                     % 保存所得向量b
             break;
         end
        i=i+1;        
    end
%    B(:,r)=b;                                % 保存所得向量b
end

%%%%%%%%%%%%%%%%%%%%%%%%%%  ICA计算的数据复原并构图  %%%%%%%%%%%%%%%%%%%%%%%%%
ICAedS=B'*Q*MixedS_bak;                     % 计算ICA后的矩阵

% 将混合矩阵重新排列并输出
figure
res1=sum(ICAedS(1:end,:),1);
plot(res1);
title('ICA重构信号')

%%%%%%%%%%%%%%%%%%%%%%%%%%  PCA计算并构图  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V,D]=eig(MixedS_cov);                          % 协方差矩阵对角化
Vtmp=zeros(size(V,1),1);
for j=1:2                                       % 选择最大的主元向量并排序
    for i=1:2
        if D(i,i)<D(i+1,i+1)
            tmp=D(i,i);Vtmp=V(:,i);
            D(i,i)=D(i+1,i+1);V(:,i)=V(:,i+1);
            D(i+1,i+1)=tmp;V(:,i+1)=Vtmp;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%  PCA求主元并显示  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
res2=sum(MixedS'*V(1:end,:),1);
plot(res2);
title('PCA重构信号')

