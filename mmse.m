clear all
P=4; L=4; iter=100000
P=4; L=2; iter=1
%P=256; L=256; iter=100
H=randn(P,L) + 1j*randn(P,L);
%sigma=0.05;
%sigma=1;
sigma=0;
I=eye(L);
x=randn(L,1) + 1j*randn(L,1);
y=H*x;
clear x

%%(3a) csiを求めるなら、以下。(ZF/MMSE)
pause(.5); tic
for i=1:iter
HHH=(H'*H+sigma*I);% _HERK % C = A'*A+ beta*C % エルミート(正定値?)
INV = (HHH)\I;     % CHESV % A*x = b % (HHH)*(INV) = (I)
G = INV*H';        % _GEMV % C = A*B'
%csi3a = 1./real(diag(INV)); % isreal(csi)=1 : heavy!!
csi3a = 1./(diag(INV)); % isreal(csi)=0
xhat3a = G * y;    % BLAS  % C = A * B
end
t3a=toc
clear HHH INV G

%%(3b) csiを求めるなら、以下。(ZF/MMSE)
pause(.5); tic
for i=1:iter
HHH=(H'*H+sigma*I);% _HERK % C = A'*A+ beta*C % エルミート(正定値?)
INV = (HHH)\I;     % CHESV % A*x = b % (HHH)*(INV) = (I)
HY = H'*y;         % _GEMV % C = A'*B
%csi3b = 1./real(diag(INV)); % isreal(csi)=1 : heavy!!
csi3b = 1./(diag(INV));
xhat3b = INV * HY; % BLAS  % C = A * B
end
t3b=toc
clear HHH HY INV

%%(3c) csiを求めるなら、以下。(ZF/MMSE)
pause(.5); tic
for i=1:iter
HHH=(H'*H+sigma*I);% _HERK % C = A'*A+ beta*C % エルミート(正定値?)
G = (HHH)\H';      % CHESV % A*X = B % (HHH)*(G) = (H')
%csi3c = 1./real(diag(G*G')); % heavy
%csi3c = 1./sum(real(conj(G).*G),2); % heavy
%csi3c = 1./sum((conj(G).*G),2); % heavy a little
csi3c = 1./diag(G*G');
xhat3c = G * y;    % _GEMV
end
t3c=toc
clear HHH G GH

%%(3c) csiを求めるなら、以下。(ZF/MMSE)
pause(.5); tic
for i=1:iter
HHH=(H'*H+sigma*I);% _HERK % C = A'*A+ beta*C % エルミート(正定値?)

%G=HHH\H'; % CHESV % A*X = B % (HHH)*(G) = (H')
%csi3c = 1./sum((conj(G).*G),2);
%xhat3c = G * y ;  % _GEMV % C = A*B

GH=H/HHH; % CHESV % X*A = B % (GH)*(HHH) = (H)
%csi3c = 1./diag(GH'*GH);
csi3c = 1./sum(conj(GH).*GH);
xhat3c = GH' * y ;% _GEMV % C = A'*B

end
t3c=toc
clear HHH G GH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff_csi = norm(csi3a -csi3b)
diff_xhat= norm(xhat3a-xhat3b)
diff_csi = norm(csi3a -csi3c)
diff_xhat= norm(xhat3a-xhat3c)
disp('   csi3a    csi3b    csi3c (top10)')
disp(real([csi3a csi3b csi3c](1:min(10,end),:)))
