%% a) and b)
clearvars
close all

% Iterations to run for
iter = 100;

% Utility values
% a)
a = 1; b = 1; c = 0; d = 0;
% b)
% a = 1; b = 0.5; c = 0; d = 0;

% Utility function
phi = [a d; c b];

% Potential function
pot = [a-c 0; 0 b-d];

nbrNodes = 5;
actions = [1 2];
nbrActions = 2;

% Initialize actions randomly
x = randi(2,nbrNodes,1);

% Inverse noise
eta = 2;

% Adjacency matrix
W = zeros(nbrNodes,nbrNodes);
W(1,2) = 1; W(1,3) = 1; W(1,5) = 1;
W(2,1) = 1; W(2,3) = 1;
W(3,1) = 1; W(3,2) = 1; W(3,4) = 1; W(3,5) = 1;
W(4,3) = 1; W(4,5) = 1;
W(5,1) = 1; W(5,3) = 1; W(5,4) = 1;

coord = [0 0; 1 0; 1 1; 0.5 1.5; 0 1]; % node coordinates

% Potential in each iteration
potential = zeros(iter,1);

% Draw links and nodes ----------------------------------------------------
for i = 1:nbrNodes
    for j = 1:nbrNodes
        if W(i,j) ~= 0
            plot(coord([i j],1),coord([i j],2),'k'); % draw link
            hold on
        end
    end
end
axis equal
for n = 1:nbrNodes
    if x(n) == 1
        rectangle('Position', ...
            [coord(n,1)-0.1 coord(n,2)-0.1 0.2 0.2], ...
            'Curvature', [1 1], 'Facecolor', 'r')
        % draw red node (action 1)
    else % x(n) == 2
        rectangle('Position', ...
            [coord(n,1)-0.1 coord(n,2)-0.1 0.2 0.2], ...
            'Curvature', [1 1], 'Facecolor', 'b')
        % draw blue node (action 2)
    end
end
% -------------------------------------------------------------------------

% Iterate over all timesteps
for k = 1:iter
    
    % Iterate over all nodes and update actions
    for m = 1:nbrNodes
        prob = zeros(nbrActions,1);
        
        % Iterate over possible actions
        for action = actions
            % Calculate utility for node 'm' and action 'action'
            util = W(m,:)*phi(action,x)';
            % Cacluate probability of current action
            prob(action) = exp(eta*util);
        end
        % Normalize probability
        prob = prob/sum(prob);
        
        % Cumulative probability vector
        c = cumsum(prob,1);
        
        % Update action of node 'm' according to calculated probability
        x(m) = find(c > rand(1),1);
    end
    
    % Calculate potential
    for m = 1:nbrNodes
        potential(k) = potential(k) + 1/2*W(m,:)*pot(x(m),x)';
    end
    


    % Draw nodes ----------------------------------------------------------
    pause(0.05)
    for n = 1:nbrNodes
        if x(n) == 1
            rectangle('Position', ...
                [coord(n,1)-0.1 coord(n,2)-0.1 0.2 0.2], ...
                'Curvature', [1 1], 'Facecolor', 'r')
            % draw red node (action 1)
        else % x(n) == 2
            rectangle('Position', ...
                [coord(n,1)-0.1 coord(n,2)-0.1 0.2 0.2], ...
                'Curvature', [1 1], 'Facecolor', 'b')
            % draw blue node (action 2)
        end
    end
    % ---------------------------------------------------------------------
end

figure
plot(potential,'-o')
xlabel('Time')
ylabel('Potential')






