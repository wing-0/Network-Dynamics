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

%% a) Particle perspective

% Start in node a (index 2)

% Steps to run each simulation for
n = 50000;

% Cumulative sum of the transition probs for a node
cumprob = cumsum(P,2);

% Array for saving the nodes for particle k at each step m = 1:n
node = zeros(n,k);
node(1,k) = 2*ones(1,k);

% Timestamp - particle k arrives at node(m,k) at time t(m,k) and leaves
% node(m,k) at t(m+1,k) (assume instantaneous travel to the next node)
t = zeros(n+1,k);

for k = 1:100
    % Time the particle should stay at the starting node - using weight of
    % the starting node as the rate
    tnext = -log(rand())/w(node(1,k));
    
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
end
