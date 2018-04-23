
W = [0 0 1 2 0 0;
     1 0 0 0 0 0;
     7 0 0 0 0 3;
     0 1 4 0 0 0;
     0 0 1 0 2 0;
     0 0 0 0 1 0];
 
D = diag(sum(W,2));
 
P = D\W;

%% a) Calculate dominant eigenvector of P'
[V,D] = eig(P');
lambda = diag(D);
pi = V(:,lambda==max(lambda));

% Normalize using L1 norm
pi = pi/sum(pi);

%% b) Random walk
close all

n = 1000.*[1 2 5 10];
cumprob = cumsum(P,2);

for k = 1:size(n,2)
    node = zeros(size(n,1),1);
    node(1) = 1;
    
    % Random walking
    for m = 1:n(k)+1
        pr = rand;
        c = cumprob(node(m),:);
        node(m+1) = find(c>pr,1);
    end
    
    figure
    plot(1:20,node(1:20),'d-')
    piest = zeros(6,1);
    for m = 1:6
        piest(m) = sum(node==m)/size(node,2);
    end
    fprintf('n = %1.0f \n',n(k))
    fprintf('Norm difference (L1): %1.3f \n',norm(pi-piest))
    fprintf('%8s %8s \n','pi','piest')
    fprintf('%8.3f %8.3f \n',[pi piest]')
    fprintf('\n')
end

%% c) Random walk (return time)
n = 10000;
cumprob = cumsum(P,2);
return_est = zeros(6,1);

for k = 1:6
    node = zeros(size(n,1),1);
    node(1) = k;
    
    % Random walking
    for m = 1:n
        pr = rand;
        c = cumprob(node(m),:);
        node(m+1) = find(c>pr,1);
    end
    
    hits = find(node==k,sum(node==k));
    return_est(k) = mean(diff(hits));
end

fprintf('%8s %8s \n','1/pi','return_est')
fprintf('%8.3f %8.3f \n',[1./pi return_est]')
fprintf('\n')

%% d) Hitting time for nodes [2 5]
S = [2 5];
nd = setdiff(1:6,S);
Phat = P(nd,nd);
x = (eye(length(nd))-Phat)\ones(length(nd),1);
Ts = zeros(6,1);
Ts(nd) = x





