function [y] = fastICA_Complex_lly_jw(x,ds)

	M = size(x,1); % Number of antenna 返回矩阵行数
	if (exist('ds','var')~=1)||isempty(ds)
		ds = M; % Number of Sources
	end

	SigLen = size(x,2);   % Signal length 返回x矩阵列数

	% % centering
	% x = x-repmat(mean(x,2),1,SigLen);

	% whitening
	Rx = x*x'./SigLen;
	[ux,dx] = eig(Rx);%求Rx的特征值，并构成特征对角矩阵
	% White_x = ux*inv(sqrt(dx))*ux'*x;
	if ~isempty(find(inv(sqrt(dx))==Inf, 1))%sqrt 平方根 inv 求逆矩阵 find 返回第一个非零元素1的索引值
		y=x;
		return;
	end
	WF = inv(sqrt(dx))*ux';  
	White_x = WF*x;

	% sig = [];
	% for ii = 1:M
	%     sig = [sig;real(White_x(ii,:));imag(White_x(ii,:))];
	% end
	sig = reshape([shiftdim(real(White_x),-1);shiftdim(imag(White_x),-1)],...
			2*M,SigLen);%shiftdim(A,1)使A的维号左移1位

	% initialize the diagonal filters
	% for ii=1:2*M
	%     W(ii,ii) = 1;
	% end
	W = eye(2*ds,2*M);%返回对角线矩阵

	NumIter = 200;
	convergence = zeros(1,NumIter);%0矩阵
	count = 1;
	% while abs(convergence-0)>1e-4
	while count<=NumIter
		Wk = W;
		y = Wk*sig;
		W = ((y.^3)*sig')./SigLen - diag(sum(3*(y.^2),2)./SigLen)*Wk;
	%     W = (sig*(y.^3)')./SigLen - 3.*W;
		W = W/norm(W,'fro');
		if sum(isnan(W(:)))~=0
			y=x;
			return;
		end
		W = ((W*W')^(-0.5))*W;
		W_tmp = mat2cell(W,2*ones(1,ds),2*ones(1,M));
		Wc = cell(size(W_tmp));
		for ii = 1:size(W_tmp,1)
			for jj = 1:size(W_tmp,2)
				c = W_tmp{ii,jj};
				W1 = c - c';
				W2 = diag([sum(diag(c)),sum(diag(c))]);
				Wc{ii,jj} = W1+W2;
			end
		end
		W = 0.5*cell2mat(Wc);
	%     W2 = diag(diag(W))+[0,1,0,0;1,0,0,0;0,0,0,1;0,0,1,0]*diag(diag(W))*[0,1,0,0;1,0,0,0;0,0,0,1;0,0,1,0];
	%     W = 0.5*(W1+W2);
		det_W = W - Wk;
		convergence(count) = trace(W'*Wk);
		count = count+1;
	%     if count>1000  break; end
	end
	%  figure(4);plot(1:count-1,convergence,'-*')
	 y_tmp = W*sig;
	 y_tmp = reshape(y_tmp,2,ds,SigLen);
	 y = y_tmp(1,:,:) + 1j*y_tmp(2,:,:);
	 y = shiftdim(y,1);

	 Wc = mat2cell(W,2*ones(1,ds),2*ones(1,M));
	 for ii = 1:numel(Wc)
	   Wc{ii} = Wc{ii}(1,1)+1j*Wc{ii}(2,1);
	 end
	 Wc = cell2mat(Wc);

	 W_Demix = Wc*WF;
end
% figure; plot(convergence)
    