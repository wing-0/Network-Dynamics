%% Initialize
clearvars
close all

%% a) - Line graph coloring

% Construct line graph adjacency matrix
nbrNodes = 10;
W = diag(ones(nbrNodes-1,1),1)+diag(ones(nbrNodes-1,1),-1);

% Graph coordinates
xy = [(1:nbrNodes)' zeros(nbrNodes,1)];

% Iterations (timesteps) to simulate
iter = 1000;

% Timestep vector
t = (1:iter)';

% Inverse noise parameter
eta = t/100;

% Colors to choose from (1: red, 2: green)
colors = [1 2];

% State distribution initialized to every node red
x = ones(nbrNodes,1);

% Cost function so that c(1,1)=1, c(2,2)=1, others=0 (1: red, 2: green)
c = eye(2);

% Potential in each timestep
potential = zeros(iter,1);

% Calculate initial potential
for m = 1:nbrNodes
    potential(1) = potential(1) + 1/2*W(m,:)*c(x(m),x)';
end

% Time stepping
for k = 2:iter
    
    % Select a node randomly
    node = randi(nbrNodes,1);
    
    % Vector for probability of updating to either color
    prob = zeros(size(colors,2),1);
        
    % Iterate over possible colors
    for col = colors
        
        % Calculate total cost for node 'node' and color 'col'
        % (weigh times cost summed over all out-neighbors)
        cost = W(node,:)*c(col,x)';
        
        % Cacluate probability of updating to current color
        prob(col) = exp(-eta(k)*cost);
    end
    
    % Normalize probability
    prob = prob/sum(prob);
    
    % Cumulative probability vector
    cu = cumsum(prob,1);
    
    % Update color of node 'node' according to calculated probability
    x(node) = find(cu > rand(1),1);
    
    % Calculate potential
    for m = 1:nbrNodes
        potential(k) = potential(k) + 1/2*W(m,:)*c(x(m),x)';
    end
    
    % Plot links and nodes
    hold off
    gplot(W,xy,'k')
    hold on
    scatter(xy(x==1,1),xy(x==1,2),'r','filled')
    scatter(xy(x==2,1),xy(x==2,2),'g','filled')
    pause(0.01)
end

% Plot potential
figure
plot(t,potential)
xlabel('Time')
ylabel('Potential')

%% b) - Wifi channel assignment

% Import adjacency matrix and coordinates
W = load('wifi.mat','-ASCII');
xy = load('coords.mat','-ASCII');
nbrNodes = size(W,1);

% Iterations (timesteps) to simulate
iter = 2000;

% Timestep vector
t = (1:iter)';

% Inverse noise parameter
eta = t/100;
% eta = 2*log(t);
% eta = exp(t/700);

% Colors to choose from (1: red; 2: green; 3: blue; 4: yellow; 5: magenta; 
% 6 :cyan; 7: white; 8: black)
colors = 1:8;

% Custom colormap for plotting nodes with color corresponding to numbers
map = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 1 1; 0 0 0];

% State distribution for each time step. Initialized randomly
x = zeros(iter,nbrNodes);
x(1,:) = randi(size(colors,2),1,nbrNodes);

% Cost function so that c(s,s) = 2, c(s,t) = 1 if |s-t|=1 and 0 otherwise.
% This is a diagonal matrix with 2 on the main diagonal, 1 om the first 
% sub- and superdiagonals and 0 elsewhere
c = 2.*eye(8) + diag(ones(7,1),1) + diag(ones(7,1),-1);

% Potential in each timestep
potential = zeros(iter,1);

% Calculate initial potential
for m = 1:nbrNodes
    potential(1) = potential(1) + 1/2*W(m,:)*c(x(1,m),x(1,:))';
end

% Plot graph and nodes
figure
hold off
gplot(W,xy,'k')
hold on
h = scatter(xy(:,1),xy(:,2),30,map(x(1,:),:),'filled');
set(h,'MarkerEdgeColor','k')

% Time stepping
for k = 2:iter
    
    % Carry over states from previous timestep
    x(k,:) = x(k-1,:);
    
    % Select a node randomly
    node = randi(nbrNodes,1);
    
    % Vector for probability of updating to either color
    prob = zeros(size(colors,2),1);
    
    % Iterate over possible colors
    for col = colors
        
        % Calculate total cost for node 'node' and color 'col'
        % (weigh times cost summed over all out-neighbors)
        cost = W(node,:)*c(col,x(k-1,:))';
        
        % Cacluate probability of updating to current color
        prob(col) = exp(-eta(k-1)*cost);
    end
    
    % Normalize probability
    prob = prob/sum(prob);
    
    % Cumulative probability vector
    cu = cumsum(prob,1);
    
    % Update color of node 'node' according to calculated probability
    x(k,node) = find(cu > rand(1),1);
    
    % Calculate potential
    for m = 1:nbrNodes
        potential(k) = potential(k) + 1/2*W(m,:)*c(x(k,m),x(k,:))';
    end
    
    % Update node colors
    set(h,'CData',map(x(k,:),:))
    pause(0.001)
end

% Plot potential
figure
plot(t,potential)
xlabel('Time')
ylabel('Potential')

% Find timesteps with minimum potential and select one of them randomly
t_min = find(potential==min(potential));
ts = t_min(randperm(length(t_min),1));

% Plot links and nodes for minimum potential solution
figure
gplot(W,xy,'k')
hold on
h = scatter(xy(:,1),xy(:,2),30,map(x(ts,:),:),'filled');
set(h,'MarkerEdgeColor','k')








