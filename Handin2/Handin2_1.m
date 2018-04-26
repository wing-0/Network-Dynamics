%% Initialize
%--------------------------------------------------------------------------

clearvars
close all

% Define transition rate matrix
Lambda = [
0 2/5 1/5 0 0;
0 0 3/4 1/4 0;
1/2 0 0 1/2 0;
0 0 1/3 0 2/3;
0 1/3 0 1/3 0];

% Weights and normalized weight matrix
w = sum(Lambda,2);
P = diag(w)\Lambda;

%% a) Return time

% Steps to run the simulations for
n = 50000;

% Cumulative sum of the transition probs for a node
cumprob = cumsum(P,2);

% Array for saving the nodes at each step m = 1:n. First node is node 1,
% but due to the size of the simulation and the memoryless property of the
% Markov chain the result should be similar regardless of the starting node
node = zeros(n,1);
node(1) = 1;

% Timestamp - particle arrives at node(n) at time t(n) and leaves
% node(n) at t(n+1) (assume instantaneous travel to the next node)
t = zeros(n+1,1);

% Time the particle should stay at the starting node - using weight of
% the starting node as the rate
tnext = -log(rand())/w(node(1));

% Random walking
for m = 2:n
    
    % Cumulative transition probs for the previous node
    c = cumprob(node(m-1),:);
    
    % Use a uniform random number to select the current node
    node(m) = find(c > rand(),1);
    
    % Timestamp for current node (previous timestamp + time at previous
    % node)
    t(m) = t(m-1) + tnext;
    
    % Time the particle should stay at current node - using weight of
    % current node as the rate
    tnext = -log(rand())/w(node(m));
end

% Timestamp for when the particle leaves the final node
t(n+1) = t(n) + tnext;

% Finding return time estimates for all nodes
treturn = zeros(size(Lambda,1),1);
for k = 1:size(Lambda,1)
    
    % Timestamps for when the particle is at node k
    tk = t(node==k);
    
    % Time between the timestamps, averaging for return time estimate
    treturn(k) = mean(diff(tk));
end

fprintf('%14s %24s \n', 'Starting node', 'Return time estimate')
fprintf('%14.0f %24.2f \n',[ (1:size(Lambda,1))' treturn]')

%% b) Theoretical return time

% Find stationary probability distribution (normalized leading eigenvector 
% of P')
[V,D] = eig(P');
lambda = diag(D);
Pi = V(:,lambda==max(lambda));

% Normalize using L1 norm
Pi = Pi/sum(Pi);

% Find invariant probability vector using equation on p. 82 in lecture
% notes
pibar = Pi./w/sum(Pi./w);

% Theoretical return time 1/pibar_i/w_i
treturn_th = 1./pibar./w;

fprintf('%14s %24s %24s \n', 'Starting node', 'Return time estimate', 'Theoretical return time')
fprintf('%14.0f %24.2f %24.2f \n',[ (1:size(Lambda,1))' treturn treturn_th]')

%% c) Hitting time (using simulation from a))

% Step indices for when the particle is at node o (index 1)
ind_o = find(node==1);

% Step indices for when the particle is at node d (index 5)
ind_d = find(node==5);

% Time between each o-hit and the following d-hit
t_od = nan(size(ind_o,1),1);

for k = 1:size(ind_o,1)-1
    % Index for the d-hits between o-hit k and o-hit k+1
    hits_d = ind_d(ind_d > ind_o(k) & ind_d < ind_o(k+1));
    
    % Save the time difference between the first d-hit from above, and
    % o-hit k. NaN if no d-hits between o-hit k and o-hit k+1
    if(isempty(hits_d))
        t_od(k) = NaN;
    else
        t_od(k) = t(hits_d(1)) - t(ind_o(k));
    end
end

% Index for the d-hits after the last o-hit
hits_d = ind_d(ind_d > ind_o(end));

% Save the time difference between the first d-hit from above, and the last
% o-hit . NaN if no d-hits after the last o-hit
if(isempty(hits_d))
    t_od(end) = NaN;
else
    t_od(end) = t(hits_d(1)) - t(ind_o(end));
end

% Mean value of the non-NaN times giving the hitting time estimate
t_hit = mean(t_od(~isnan(t_od)));

%% c) Theoretical hitting time
lb = zeros(5,1);

ub = Inf*ones(5,1);
ub(5) = 0;

t_hit_th = lsqlin(eye(5)-P,1./w,[],[],[],[],lb,ub);

% Correcter solution for node 1 was given for parameters:
% lb = -Inf*ones(5,1); lb(5) = 0
% Although, this gave negative values for nodes 2-4



