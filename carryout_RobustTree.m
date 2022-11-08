clear

%% load data
load tree_300.mat
% Normalization
XX = (X-min(min(X)))./(max(max(X))-min(min(X)));
D = XX;
[m,n]=size(D);

%%  If there are no missing values
omega=1:(m*n);
omegaC=[];

%% If there exists missing entries
% miss=2;
% omega=find(D~=miss);
% omegaC=find(D==miss);

%% Adaptive parameter selections
lambda=1/sqrt(max(m,n))*(1+3*length(omegaC)/(m*n));
 
varD=var(D);
sigma= (sum(varD))/max(m,n);
gamma=m/(n*sqrt(n));
theta=sqrt(m)/(n*sqrt(n));

%% Run RPCAtree
tol = 1e-10;
maxIter = 1000;
tic
[A1, E1, C, R, B] = RPCAtree(D, omega, lambda, theta, gamma, sigma, tol, maxIter);
toc

%% Plot tree for 2D data
plot(C(:,1),C(:,2),'ko');
hold on
for i=1:size(C,1)
 for j=1:size(C,1)
if B(i,j)>0 
   plot([C(i,1),C(j,1)],[C(i,2),C(j,2)],'ko-','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k' );
end
end
 end
plot(D(:,1),D(:,2),'.b','MarkerSize',8);
hold off


