%% Prep
clearvars
close all

Lambda = [
0 1/30 1/15 1/60
1/60 0 1/10 1/100
1/25 1/10 0 1/50
1/100 1/10 1/10 0];

%% a) Invariant prob distr

w = sum(Lambda,2);
ws = max(w);

% Non diagonal elements of Q
Q = Lambda/ws;

% Diagonal elements of Q
Q(1:size(Q,1)+1:end) = (ones(size(Q,1),1)-sum(Q,2)+Q(1:size(Q,1)+1:end)')';

[V,D] = eig(Q');
lambda = diag(D);
pibar = V(:,lambda==max(lambda));

% Normalize using L1 norm
pibar = pibar/sum(pibar);

%% b) Jump chain
close all

n = 10000;
cumprob = cumsum(Q,2);

node = zeros(n,1);
node(1) = 1;

% Timestamp - particle arrives at position node(k) at time t(k) and leaves
% node(k) at t(k+1) (assume instantaneous travel to the next node)
t = zeros(n+1,1);

% Random walking
for k = 2:n
    pr = rand;
    c = cumprob(node(k-1),:);
    node(k) = find(c>pr,1);
    t(k) = t(k-1)-log(rand())/ws;
end

% Time when the particle leaves node(n)
t(end) = t(n)-log(rand())/ws;

% Plotting
figure
plot(1:20,node(1:20),'d-')
pibarest = zeros(size(Lambda,1),1);

% Time at each of the nodes
time = diff(t);

for k = 1:size(pibarest,1)
    pibarest(k) = sum(time(node==k))/t(end);
end
fprintf('Norm difference (L1): %1.3f \n',norm(pibar-pibarest))
fprintf('%8s %8s \n','pibar','pibarest')
fprintf('%8.3f %8.3f \n',[pibar pibarest]')
fprintf('\n')

%% c) Rate w_i

P = diag(w)\Lambda;
n = 10000;
cumprob = cumsum(P,2);

node = zeros(n,1);
node(1) = 1;

% Timestamp - particle arrives at position node(k) at time t(k) and leaves
% node(k) at t(k+1) (assume instantaneous travel to the next node)
t = zeros(n+1,1);

% Random walking
for k = 2:n
    pr = rand;
    c = cumprob(node(k-1),:);
    node(k) = find(c>pr,1);
    t(k) = t(k-1)-log(rand())/w(node(k-1));
end

% Time when the particle leaves node(n)
t(end) = t(n)-log(rand())/w(node(n));

% Plotting
figure
plot(1:20,node(1:20),'d-')
pibarest = zeros(size(Lambda,1),1);

% Time at each of the nodes
time = diff(t);

for k = 1:size(pibarest,1)
    pibarest(k) = sum(time(node==k))/t(end);
end
fprintf('Norm difference (L1): %1.3f \n',norm(pibar-pibarest))
fprintf('%8s %8s \n','pibar','pibarest')
fprintf('%8.3f %8.3f \n',[pibar pibarest]')
fprintf('\n')

%% d) Verify relation

% Find pi
[V,D] = eig(P');
lambda = diag(D);
Pi = V(:,lambda==max(lambda));

% Normalize using L1 norm
Pi = Pi/sum(Pi);

pibar2 = Pi./w/sum(Pi./w);

fprintf('Norm difference (L1): %1.3f \n',norm(pibar-pibar2))
fprintf('%8s %8s \n','pibar','pibar (calculated)')
fprintf('%8.3f %8.3f \n',[pibar pibar2]')
fprintf('\n')

%% Ice cream profits

P = diag(w)\Lambda;
n = 10000;
cumprob = cumsum(P,2);

% Profits per unit time for [sunny rainy cloudy snowy]'
f = [10 2 1 0]';

node = zeros(n,1);
node(1) = 1;

% Timestamp - particle arrives at position node(k) at time t(k) and leaves
% node(k) at t(k+1) (assume instantaneous travel to the next node)
t = zeros(n+1,1);

tnext = -log(rand())/w(1);

% Profits at timestamps
profits = zeros(n+1,1);

% Random walking
for k = 2:n
    pr = rand;
    c = cumprob(node(k-1),:);
    node(k) = find(c>pr,1);
    t(k) = t(k-1) + tnext;
    profits(k) = f(node(k))*tnext;
    
    tnext = -log(rand())/w(node(k-1));
end

% Time when the particle leaves node(n)
t(end) = t(n)+tnext;

profits(end) = f(node(end))*tnext;

% Plotting
figure
plot(t,cumsum(profits))

% Average profits and theoretical value
avgprofits = sum(profits)/t(end);
thprofits = sum(pibar.*f);

fprintf('%28s %28s \n','Theoretical value','Simulated value')
fprintf('%28.3f %28.3f \n',thprofits, avgprofits)










