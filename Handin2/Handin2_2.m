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

% Number of time units to run simulation for
tmax = 60;

% Cumulative sum of the transition probs for a node
cumprob = cumsum(P,2);

% Vector for the times of Poisson ticks (any node). Preallocating for 10000
% ticks, but actual number may vary
t = NaN(10000,1);
t(1) = 0;

% Number of particles in each node at each timestep, initialized with 100
% particles in node o (index 1) at time t = 0
nump = NaN(size(t,1),size(P,2));
nump(1,:) = zeros(1,size(P,2));
nump(1,1) = 100;

% Index corresponding to time
k = 1;

% Set time until the next tick for all nodes
tnext = -log(rand(1,size(P,2)))./nump(k,:)./w';

% Time stepping
while t(k) + min(tnext) < tmax
    
    % Increase time index and set time for the tick to come
    k = k+1;
    t(k) = t(k-1) + min(tnext);
    
    % Set the number of particles in nodes to be equal to the previous time
    nump(k,:) = nump(k-1,:);
    
    % Find the node with the least time until tick
    ticknode = find(tnext==min(tnext));
    
    % Cumulative transition probs for the ticking node
    c = cumprob(ticknode,:);
    
    % Use a uniform random number to select the node to move a particle to
    movenode = find(c > rand(),1);
    
    % Change particle numbers for the two selected nodes
    nump(k,ticknode) = nump(k,ticknode)-1;
    nump(k,movenode) = nump(k,movenode)+1;
    
    % Draw new time until the next tick for all nodes. Exponential
    % distribution is memoryless, so the time until next tick is still
    % exponentially distributed (though the rate could be changed by moving
    % particles)
    tnext = -log(rand(1,size(P,2)))./nump(k,:)./w';
    
end

% Truncate matrices so they only have the size of the actual number of
% Poisson ticks
t = t(1:k);
nump = nump(1:k,:);

% Calculate average number of particles in nodes at the end of the
% simulation (last 5% of the ticks)
nump_end = nump(ceil(0.95*k):k,:);
nump_avg = mean(nump_end,1);
nump_std = std(nump_end,1);

fprintf('%10s %20s \n','Node index','Avg. particles in node at end')
fprintf('%10.0f %20.2f \n',[1:5; nump_avg])

% Plotting
plot(t,nump)
xlabel('Time')
ylabel('Number of particles')
legend(['Node o';'Node a';'Node b';'Node c';'Node d'])




