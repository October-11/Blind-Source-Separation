
%%%%%%%%%%%%%%%%%%%%%%%%%%  ��ʼ��  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
clc;clear all;close all;

 
%%%%%%%%%%%%%%  ����ԭʼ��Ƶ�ļ�����ϣ�����������Ƶ  %%%%%%%%%%%%%%%%%% 

[S]=open();

d=size(S);
Sweight=randn(size(S,1));             % ȡһ���������Ϊ�źŻ�ϵ�Ȩ���� 
MixedS=Sweight*S;                      % �õ���������źž��� 
 
% ����Ͼ����������в���� 
subplot(3,3,4),plot(MixedS(1,:)),title('����ź�1'), axis([0,50000,-0.5,0.5]);
subplot(3,3,5),plot(MixedS(2,:)),title('����ź�2'), axis([0,50000,-0.5,0.5]);
subplot(3,3,6),plot(MixedS(3,:)),title('����ź�3'), axis([0,50000,-0.5,0.5]);

% figure,plot(MixedS(1,:)),title('����ź�1'),% axis([0,1000,-4,4]);
% figure,plot(MixedS(2,:)),title('����ź�2'), %axis([0,1000,-4,4]);
% figure,plot(MixedS(3,:)),title('����ź�3'), %axis([0,1000,-4,4]);

fs=44100;
audiowrite('123.WAV',MixedS(1,:),fs);  
audiowrite('223.WAV',MixedS(2,:),fs); 
audiowrite('323.WAV',MixedS(3,:),fs);

 
MixedS_bak=MixedS;                      % ����Ϻ�����ݱ��ݣ��Ա��ڻָ�ʱֱ�ӵ��� 


%%%%%%%%%%%%%%%%%%%%%%%%%%  ��׼��  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
MixedS_mean=zeros(3,1); 
for i=1:3 
    MixedS_mean(i)=mean(MixedS(i,:)); 
end                                        % ����MixedS�ľ�ֵ 
 
for i=1:3 
    for j=1:size(MixedS,2) 
        MixedS(i,j)=MixedS(i,j)-MixedS_mean(i);    %ȥ��ֵ
    end 
end                                        
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%  �׻�  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
MixedS_cov=cov(MixedS');                    % covΪ��Э����ĺ��� 
[E,D]=eig(MixedS_cov);                      % ��ͼƬ�����Э�������������ֵ�ֽ� E����������D����ֵ
Q=inv(sqrt(D))*(E)';                        % QΪ�׻����� 
MixedS_white=Q*MixedS;                      % MixedS_whiteΪ�׻������Ƶ��Ϣ���� 
IsI=cov(MixedS_white');                     % IsIӦΪ��λ��  
% figure;
% subplot(131);plot(MixedS_white(1,:)),title('�׻�����ź�1'),  %axis([0,1000,-2,2]);
% subplot(132);plot(MixedS_white(2,:)),title('�׻�����ź�2'),   axis([0,50000,-5,5]);
% subplot(133);plot(MixedS_white(3,:)),title('�׻�����ź�3'),  %axis([0,1000,-2,2]);
 
%%%%%%%%%%%%%%%%%%%%%%%%��FASTICA�㷨  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
X=MixedS_white;                            % �����㷨����X���в��� 
[VariableNum,SampleNum]=size(X); 
numofIC=VariableNum;                       % �ڴ�Ӧ���У�����Ԫ�������ڱ������� 
B=zeros(numofIC,VariableNum) ;            % ��ʼ��������w�ļĴ����,B=[b1  b2  ...   bd] 
for r=1:numofIC                            % ������ȡÿһ������Ԫ 
    i=1;
    maxIterationsNum=100;                  % ����������������������ÿ�������������Ե������������˴����� 
   % IterationsNum=0; 
    b=rand(numofIC,1);                     % �������b��ֵ 
   
    b=b/norm(b);                           % ��b��׼�� 
    while i<=maxIterationsNum+1 
        if i == maxIterationsNum           % ѭ���������� 
            fprintf('\n��%d������%d�ε����ڲ���������', r,maxIterationsNum); 
            break; 
        end 
        bOld=b;                        
        t=X'*b; 
        g=(exp(2.*t)-1)./(exp(2.*t)+1); 
        dg=4*exp(2.*t)./(exp(2.*t)+1).^2;   %����
        b=(X*g)/SampleNum-mean(dg)*b; 
                                            % ���Ĺ�ʽ
        b=b-B*B'*b;                         % ��b������ 
        b=b/norm(b); 
        
        if abs(abs(b'*bOld)-1)<1e-9         % ����������򱣴�b 
               B(:,r)=b; 
               break; 
        end 
        i=i+1;         
     end 
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%  ���ݸ�ԭ  %%%%%%%%%%%%%%%%%%%%%%%% 
 
ICAedS=B\X;                    % x(t)=B*s(t);s(t)=B\x(t); (���)
 
% ����Ͼ����������в���� 
subplot(3,3,7),plot(ICAedS(1,:)),title('ICA����ź�1'),  %axis([0,1000,-2,2]);
subplot(3,3,8),plot(ICAedS(2,:)),title('ICA����ź�2'),   axis([0,50000,-5,5]);
subplot(3,3,9),plot(ICAedS(3,:)),title('ICA����ź�3'),  %axis([0,1000,-2,2]);

% figure,plot(ICAedS(1,:)),title('ICA����ź�1'),  %axis([0,1000,-2,2]);
% figure,plot(ICAedS(2,:)),title('ICA����ź�2'),   axis([0,50000,-5,5]);
% figure,plot(ICAedS(3,:)),title('ICA����ź�3'),  %axis([0,1000,-2,2]);

%fs=44100;
audiowrite('12345.WAV',ICAedS(1,:),fs);  
audiowrite('22345.WAV',ICAedS(2,:),fs); 
audiowrite('32345.WAV',ICAedS(3,:),fs);
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%  PCA���㲢��ͼ  %%%%%%%%%%%%%%%%%%%%%%%% 