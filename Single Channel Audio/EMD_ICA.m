
%%%%%%%%%%%%%%%%%%%%%%%%%%  ��ʼ��  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all;close all;

%%%%%%%%%%%%%%  ����ԭʼͼ�񣬻�ϣ���������ͼ��  %%%%%%%%%%%%%%%%%%

% ������ǰ��ԭʼͼƬ����ʾ
I1=audioread ('washmix3.wav')';
figure
plot(I1),title('��������1'),

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ����EMD�ֽ�  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% ��EMD�ֽ��ķ�����ɾ���
S=[imf(1,:);imf(2,:);imf(3,:);imf(4,:);imf(5,:);imf(6,:);imf(7,:);imf(8,:);imf(9,:);imf(10,:);imf(11,:);imf(12,:);imf(13,:);imf(14,:)]; 
                                       % ���S_all��һ���������������������ľ���
Sweight=rand(size(S,1));               % ȡһ���������Ϊ�źŻ�ϵ�Ȩ����
MixedS=Sweight*S;                      % �õ�����ͼ��Ļ���źž���
MixedS_bak=MixedS;                     % ����Ϻ�����ݱ��ݣ��Ա��ڻָ�ʱֱ�ӵ���
%%%%%%%%%%%%%%%%%%%%%%%%%%  ��׼��  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MixedS_mean=zeros(14,1);
for i=1:14
    MixedS_mean(i)=mean(MixedS(i,:));
end                                        % ����MixedS�ľ�ֵ

for i=1:14
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
figure
res1=sum(ICAedS(1:end,:),1);
plot(res1);
title('ICA�ع��ź�')

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
figure(4)
res2=sum(MixedS'*V(1:end,:),1);
plot(res2);
title('PCA�ع��ź�')

