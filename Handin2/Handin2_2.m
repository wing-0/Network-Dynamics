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

% Number of particles to simulate
part = 100;

% Steps to run each simulation for
n = 20000;

% Cumulative sum of the transition probs for a node
cumprob = cumsum(P,2);

% Array for saving the nodes for particle k at each step m = 1:n
node = zeros(n,part);
node(1,:) = 2*ones(1,part);

% Timestamp - particle k arrives at node(m,k) at time t(m,k) and leaves
% node(m,k) at t(m+1,k) (assume instantaneous travel to the next node)
t = zeros(n+1,part);

% Array for saving average return times of all particles
treturn = zeros(1,part);

% Time the particles should stay at the starting node - using weight of
% the starting node as the rate
tnext = -log(rand(1,part))./w(node(1,:))';

% Random walking
for m = 2:n
    
    % Cumulative transition probs for the previous nodes
    c = cumprob(node(m-1,:),:)';
    
    % Use a uniform random number to select the current nodes for all
    % particles. First, all nodes with cumprob larger than a random number 
    % is found for every node (individual rand() for the nodes). The
    % indices for these are saved in 'row','col'. The index of the first 
    % occurence of a particle index in col is saved, and the corresponding
    % nodes are found in row
    [row,col] = find(c > rand(1,100));
    [~,Irow,~] = unique(col);
    node(m,:) = row(Irow)';

    % Timestamp for current nodes (previous timestamp + time at previous
    % nodes)
    t(m,:) = t(m-1,:) + tnext;
    
    % Time the particle should stay at current node - using weight of
    % current node as the rate
    tnext = -log(rand(1,part))./w(node(m,:))';
end

% Timestamp for when the particles leave the final node
t(n+1,:) = t(n,:) + tnext;

% Iterate over particles, calculating return time for each
for k = 1:part
    % Finding return time estimate for node a (index 2)
    % Timestamps for when the particle is at node 2
    tk = t(:,k);
    tka = tk(node(:,k)==2);
    
    % Time between the timestamps, averaging for return time estimate
    treturn(k) = mean(diff(tka));
end

% Calculating average of the return times of all nodes
treturn_avg = mean(treturn);

fprintf('Average return time (node a): %1.2f \n',treturn_avg)

%% b) Node perspective






