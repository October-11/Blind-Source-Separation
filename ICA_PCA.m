
%%%%%%%%%%%%%%%%%%%%%%%%%%  初始化  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all;close all;

%%%%%%%%%%%%%%  读入原始图像，混合，并输出混合图像  %%%%%%%%%%%%%%%%%%

% 读入混合前的原始图片并显示
I1=audioread ('wash5.wav')';
I2=audioread ('dragen5.wav')';
I3=audioread ('music5.wav')';
subplot(4,3,1),plot(I1),title('输入声音1'),
subplot(4,3,2),plot(I2),title('输入声音2'),
subplot(4,3,3),plot(I3),title('输入声音3'),

% 将其组成矩阵
S=[I1;I2;I3];                          % 图片个数即为变量数，图片的像素数即为采样数
                                       % 因此S_all是一个变量个数＊采样个数的矩阵
Sweight=rand(size(S,1));               % 取一随机矩阵，作为信号混合的权矩阵
MixedS=Sweight*S;                      % 得到三个图像的混合信号矩阵

% 将混合矩阵重新排列并输出
subplot(4,3,4),plot(MixedS(1,:)),title('混合声音1'),
subplot(4,3,5),plot(MixedS(2,:)),title('混合声音2'),
subplot(4,3,6),plot(MixedS(3,:)),title('混合声音3'),

MixedS_bak=MixedS;                         % 将混合后的数据备份，以便在恢复时直接调用
%%%%%%%%%%%%%%%%%%%%%%%%%%  标准化  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MixedS_mean=zeros(3,1);
for i=1:3
    MixedS_mean(i)=mean(MixedS(i,:));
end                                        % 计算MixedS的均值

for i=1:3
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
subplot(4,3,7),plot(ICAedS(1,:)),title('ICA解混声音1'),
subplot(4,3,8),plot(ICAedS(2,:)),title('ICA解混声音2'),
subplot(4,3,9),plot(ICAedS(3,:)),title('ICA解混声音3'),

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
t1=(MixedS'*V(:,1))';
t2=(MixedS'*V(:,2))';   
t3=(MixedS'*V(:,3))';

subplot(4,3,10),plot(t1),title('PCA解混声音1'),
subplot(4,3,11),plot(t2),title('PCA解混声音2'),
subplot(4,3,12),plot(t3),title('PCA解混声音3'),