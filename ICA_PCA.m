
%%%%%%%%%%%%%%%%%%%%%%%%%%  ��ʼ��  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all;close all;

%%%%%%%%%%%%%%  ����ԭʼͼ�񣬻�ϣ���������ͼ��  %%%%%%%%%%%%%%%%%%

% ������ǰ��ԭʼͼƬ����ʾ
I1=audioread ('wash5.wav')';
I2=audioread ('dragen5.wav')';
I3=audioread ('music5.wav')';
subplot(4,3,1),plot(I1),title('��������1'),
subplot(4,3,2),plot(I2),title('��������2'),
subplot(4,3,3),plot(I3),title('��������3'),

% ������ɾ���
S=[I1;I2;I3];                          % ͼƬ������Ϊ��������ͼƬ����������Ϊ������
                                       % ���S_all��һ���������������������ľ���
Sweight=rand(size(S,1));               % ȡһ���������Ϊ�źŻ�ϵ�Ȩ����
MixedS=Sweight*S;                      % �õ�����ͼ��Ļ���źž���

% ����Ͼ����������в����
subplot(4,3,4),plot(MixedS(1,:)),title('�������1'),
subplot(4,3,5),plot(MixedS(2,:)),title('�������2'),
subplot(4,3,6),plot(MixedS(3,:)),title('�������3'),

MixedS_bak=MixedS;                         % ����Ϻ�����ݱ��ݣ��Ա��ڻָ�ʱֱ�ӵ���
%%%%%%%%%%%%%%%%%%%%%%%%%%  ��׼��  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MixedS_mean=zeros(3,1);
for i=1:3
    MixedS_mean(i)=mean(MixedS(i,:));
end                                        % ����MixedS�ľ�ֵ

for i=1:3
    for j=1:size(MixedS,2)
        MixedS(i,j)=MixedS(i,j)-MixedS_mean(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%  �׻�  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MixedS_cov=cov(MixedS');                    % covΪ��Э����ĺ���
[E,D]=eig(MixedS_cov);                      % ��ͼƬ�����Э�������������ֵ�ֽ�
Q=inv(sqrt(D))*(E)';                        % QΪ�׻�����
MixedS_white=Q*MixedS;                      % MixedS_whiteΪ�׻����ͼƬ����
IsI=cov(MixedS_white');                     % IsIӦΪ��λ��            

%%%%%%%%%%%%%%%%%%%%%%%%��FASTICA�㷨  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=MixedS_white;                            % �����㷨����X���в���
[VariableNum,SampleNum]=size(X);
numofIC=VariableNum;                       % �ڴ�Ӧ���У�����Ԫ�������ڱ�������
B=zeros(numofIC,VariableNum);              % ��ʼ��������w�ļĴ����,B=[b1  b2  ...   bd]
for r=1:numofIC
    i=1;maxIterationsNum=100;               % ����������������������ÿ�������������Ե������������˴�����
    IterationsNum=0;
    b=rand(numofIC,1)-.5;                  % �������b��ֵ
    b=b/norm(b);                           % ��b��׼�� norm(b):����Ԫ��ƽ���Ϳ�����
    while i<=maxIterationsNum+1
        if i == maxIterationsNum           % ѭ����������
            fprintf('\n��%d������%d�ε����ڲ���������', r,maxIterationsNum);
            break;
        end
        bOld=b;                          
        a2=1;
        u=1;
        t=X'*b;
        g=t.*exp(-a2*t.^2/2);
        dg=(1-a2*t.^2).*exp(-a2*t.^2/2);
        b=((1-u)*t'*g*b+u*X*g)/SampleNum-mean(dg)*b;
                                           % ���Ĺ�ʽ���μ����۲��ֹ�ʽ2.52
        b=b-B*B'*b;                        % ��b������
        b=b/norm(b); 
        if abs(abs(b'*bOld)-1)<1e-9        % �����������
             B(:,r)=b;                     % ������������b
             break;
         end
        i=i+1;        
    end
%    B(:,r)=b;                                % ������������b
end

%%%%%%%%%%%%%%%%%%%%%%%%%%  ICA��������ݸ�ԭ����ͼ  %%%%%%%%%%%%%%%%%%%%%%%%%
ICAedS=B'*Q*MixedS_bak;                     % ����ICA��ľ���

% ����Ͼ����������в����
subplot(4,3,7),plot(ICAedS(1,:)),title('ICA�������1'),
subplot(4,3,8),plot(ICAedS(2,:)),title('ICA�������2'),
subplot(4,3,9),plot(ICAedS(3,:)),title('ICA�������3'),

%%%%%%%%%%%%%%%%%%%%%%%%%%  PCA���㲢��ͼ  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V,D]=eig(MixedS_cov);                          % Э�������Խǻ�
Vtmp=zeros(size(V,1),1);
for j=1:2                                       % ѡ��������Ԫ����������
    for i=1:2
        if D(i,i)<D(i+1,i+1)
            tmp=D(i,i);Vtmp=V(:,i);
            D(i,i)=D(i+1,i+1);V(:,i)=V(:,i+1);
            D(i+1,i+1)=tmp;V(:,i+1)=Vtmp;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%  PCA����Ԫ����ʾ  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1=(MixedS'*V(:,1))';
t2=(MixedS'*V(:,2))';   
t3=(MixedS'*V(:,3))';

subplot(4,3,10),plot(t1),title('PCA�������1'),
subplot(4,3,11),plot(t2),title('PCA�������2'),
subplot(4,3,12),plot(t3),title('PCA�������3'),