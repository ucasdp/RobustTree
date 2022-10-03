clear
tic
load ./tree_300.mat
XX = (X-min(min(X)))./(max(max(X))-min(min(X)));
D = XX;
[m,n]=size(D);
omega=1:(m*n);
omegaC=[];
% miss=2;
% omega=find(D~=miss);
% omegaC=find(D==miss);
lambda=1/sqrt(max(m,n))*(1+3*length(omegaC)/(m*n));
 
varD=var(D);
 sigma= (sum(varD))/max(m,n);
gamma=m/(n*sqrt(n));
theta=sqrt(m)/(n*sqrt(n));

tol = 1e-10;
maxIter = 1000;

% initialize
Y = D;
norm_two = svds(Y, 1);
%norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

A1 = zeros(m, n);
A1(omega) = D(omega);
E1 = zeros(m, n);

K = m; 
%[idx,C1]= kmeans(A1,K);
C1= A1;
R = zeros( m, K);

mu = 1.25/norm_two;
mu_bar = mu * 1e7;
r = 1.1;
r1 = 1.1;
d_norm = norm(D, 'fro');

iter = 0;
total_svd = 0;
converged = false;
stopCon = 1;
sv = 10;


while ~converged       
    iter = iter + 1;
    
    dis = zeros(K, K);
    for i = 1:size(C1,1)
        dis(i,:) = diag((repmat(C1(i,:), K, 1) - C1) * (repmat(C1(i, :), K, 1) - C1)');
    end

   dis(find(dis==0))=10^(-5);
   
   for i = 1:size(C1,1)
        dis(i,i) = 0;
   end
    
   dis = sparse(dis');
   [ST, pred] = graphminspantree(dis);
   B= zeros(K,K);
   for i=1: K
       if pred(i)~=0
           B(i,pred(i))=1;
       end   
   end
   B = B + B';
   L = diag(B * ones(K,1)) - B;

   for i = 1:m
       for j = 1:K
           R(i, j) = exp(- sum((A1(i,:)-C1(j,:)).^2)/sigma);
       end
   end
   rr=var(R');
   for i=1:m
       if rr(1,i)==0
           R(i,:)=1/m*ones(1,K);
       end
       if rr(1,i)~=0
           R(i,:)=R(i,:)./sum(R(i,:));
       end
   end
   
    %R = R./repmat(sum(R, 2), 1, K);
   
    Gamma = diag(ones(1,m) * R); 
    temp_T = D - A1 + (1/mu)*Y;
    E1 = temp_T;
    EE = max(temp_T - lambda/mu, 0);
    EE = EE+min(temp_T + lambda/mu, 0);
    E1(omega) = EE(omega);
    
    Sq = R * (theta / gamma * L + Gamma) ^ (-1/2);
    Lq = norm((mu+2*gamma) .* eye(m) - 2* gamma * Sq * Sq');
    Gk = A1 -  1/Lq * (-Y + mu * (A1 - D + E1) + 2 * gamma * A1 - 2 * gamma *  Sq * Sq' * A1);

    if choosvd(n, sv) == 1
        [U S V] = svds(Gk, sv);
        %[U S V] = lansvd(D - E1 + (1/mu)*Y, sv, 'L');
    else
        [U S V] = svd(Gk, 'econ');
    end
    diagS = diag(S);
    lenS = length(find(diagS > 1/ Lq));
    if lenS < sv
        sv = min(lenS + 1, n);
    else
        sv = min(lenS + round(0.05*n), n);
    end  

    A1 = U(:, 1:lenS) * diag(diagS(1:lenS) - 1 / Lq) * V(:, 1:lenS)';    
    C1 =  (theta / gamma * L + Gamma) ^ (-1) * R' * A1;
    %A1(find(A1<0))=0;
    total_svd = total_svd + 1;
    Z = D - A1 - E1;
    
    Y = Y + mu*Z;
    mu = min(mu*r, mu_bar);
    %sigma = sigma/r1;
        
    %% stop condition  
    stopCon = norm(Z, 'fro') / d_norm;
    if stopCon < tol
        converged = true;
    end    
      
    if ~converged && iter >= maxIter
        disp('Reach maximum iterations') ;
        converged = 1;       
    end
end

[max_R,cluster_index]=max(R,[],2);
C1=roundn(C1,-2);
CGTM =unique(C1(unique(cluster_index),:),'rows');

nc=size(CGTM,1);

CC_dis=zeros(nc,nc);
for i=1:nc
     CC_dis(i,:) = diag((repmat(CGTM(i,:), nc, 1) - CGTM) * (repmat(CGTM(i, :), nc, 1) - CGTM)');
end

CC_dis = sparse(CC_dis');
[ST, pred] = graphminspantree(CC_dis);
B= zeros(nc,nc);
 for i=1: nc
    if pred(i)~=0
        B(i,pred(i))=1;
    end   
end
B = B + B';
toc

plot(CGTM(:,1),CGTM(:,2),'ko');
hold on
for i=1:nc
 for j=1:nc
if B(i,j)>0 
   plot([CGTM(i,1),CGTM(j,1)],[CGTM(i,2),CGTM(j,2)],...
       'ko-','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k' );
end
end
 end
plot(D(:,1),D(:,2),'.b','MarkerSize',8);
hold off


