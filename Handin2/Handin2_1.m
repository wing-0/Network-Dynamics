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
n = 10000;

% Cumulative sum of the transition probs for a node
cumprob = cumsum(P,2);

% Iterate over all nodes to find a return time for each
for k = 1:size(Lambda,1)
    
    % Array for saving the nodes at each step m = 1:n
    node = zeros(n,1);
    node(1) = k;
    
    % Timestamp - particle arrives at node(n) at time t(n) and leaves
    % node(n) at t(n+1) (assume instantaneous travel to the next node)
    t = zeros(n+1,1);
    
    % Time the particle should stay at the starting node - using weight of 
    % the starting node as the rate
    tnext = -log(rand())/w(k);
    
    % Random walking
    for m = 2:n
        
        % Cumulative transition probs for the previous node
        c = cumprob(node(m-1),:);
        
        % Use a uniform random number to select this node
        node(m) = find(c > rand(),1);
        
        % Timestamp for this node (previous timestamp + time at previous
        % node)
        t(m) = t(m-1) + tnext;
        
        % Time the particle should stay at this node - using weight of this
        % node as the rate
        tnext = -log(rand())/w(node(m));
    end
    
end



