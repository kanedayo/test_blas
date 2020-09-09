%octave
clear all

%% Base-Serialization-Base : RowMajor
SeqToColM=@(X,ldX) reshape(X(:),ldX,[]);
%SeqToRowM=@(X,ldX) reshape(X(:),ldX,[]).';
ColM2Mat=@(X,ldX) SeqToColM(X,ldX);
RowM2Mat=@(X,ldX) SeqToColM(X,ldX).';
Mat2ColM=@(X) X(:).';
Mat2RowM=@(X) X.'(:).';

ColMajor=@(X,ldX,T) T(ColM2Mat(X,ldX));
RowMajor=@(X,ldX,T) T(RowM2Mat(X,ldX));

NT=@(X) (X);
TR=@(X) transpose(X);

M=3; N=4; K=2;
A=[1 2;1 2;1 2];
B=[1 2 3 4;1 2 3 4];
C=zeros(M,N);
%C=[1 2 3 4;1 2 3 4;1 2 3 4];

%% Matrix To ColM/RowM
%AcolM=Mat2ColM(A);ldAcolM=size(A,1); clear A; % Mat2ColM(A)
%BcolM=Mat2ColM(B);ldBcolM=size(B,1); clear B; % Mat2ColM(B)
%CcolM=Mat2ColM(C);ldCcolM=size(C,1); clear C; % Mat2ColM(C)
ArowM=Mat2RowM(A);ldArowM=size(A,2); clear A; % Mat2RowM(A)
BrowM=Mat2RowM(B);ldBrowM=size(B,2); clear B; % Mat2RowM(B)
CrowM=Mat2RowM(C);ldCrowM=size(C,2); clear C; % Mat2RowM(C)
%debug_A = RowMajor(ArowM,ldArowM,NT)
%debug_B = RowMajor(BrowM,ldBrowM,NT)
%debug_C = RowMajor(CrowM,ldCrowM,NT)

%% C = A*B + C %% RowMajor
ord=RowMajor;
C = ord(ArowM,ldArowM,NT) * ord(BrowM,ldBrowM,NT)  + ord(CrowM,ldCrowM,NT)
CrowM=Mat2RowM(C);ldCrowM=size(C,2); clear C; % Mat2RowM(C)


%% C = tr(tr(A)) %% RowMajor -> ColMajor
ord=ColMajor;
E = eye(N);
EcolM=Mat2ColM(E);ldEcolM=size(E,1); clear E; % Mat2ColM(E)
C = ord(CrowM,ldCrowM,TR) * ord(EcolM,ldEcolM,NT)
CcolM=Mat2ColM(C);ldCcolM=size(C,1); clear C; % Mat2ColM(C)

%% C = A*B + C %% ColMajor
ord=ColMajor;
C = ord(ArowM,ldArowM,TR) * ord(BrowM,ldBrowM,TR)  + ord(CcolM,ldCcolM,NT)
CcolM=Mat2ColM(C);ldCcolM=size(C,1); clear C;

%=============================

%% CHECK %% ColMajor / RowMajor
ord=ColMajor;
Cc= ord(ArowM,ldArowM,TR) * ord(BrowM,ldBrowM,TR)  + ord(CcolM,ldCcolM,NT);

ord=RowMajor;
Cr= ord(ArowM,ldArowM,NT) * ord(BrowM,ldBrowM,NT)  + ord(CcolM,ldCcolM,TR);

assert(Cc==Cr) % Ccol==Crow
