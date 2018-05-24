clearvars
close all

% Load data
citation = load('citation.mat','-ascii');
W = spconvert(citation);

% Count number of nodes (articles) and links (citations)
nodes = length(W);
links = length(find(W));

%% a) - Erdos-Renyi

% Probability for a link between any two nodes i and j
p = links/(nchoosek(nodes,2)*2);

% Generate a random square matrix, only keep values smaller than p
% This is the generated weight matrix (also remove diagonal to avoid self
% loops)
W_er = rand(nodes) < p;
W_er = sparse(W_er);
W_er = W_er - diag(diag(W_er));

%% b) - Configuration model

W_cm = zeros(nodes);

% Out- and in-degree vectors
w_out = sum(W,2);
w_in = sum(W,1);

% Iterate over all nodes
for k = 1:nodes
    % Connect one link for each unit of out-degree
    while w_out(k) > 0
        
        % Identify valid targets
        targets = 1:nodes;
        targets = targets(w_in > 0);
        
        % Select one of the valid targets at random and connect a link
        t = targets(randi(length(targets)));
        W_cm(k,t) = 1;
        
        % Remove one from out-degree and in-degree respectively
        w_out(k) = w_out(k) - 1;
        w_in(t) = w_in(t) - 1;
    end
end

W_cm = sparse(W_cm);

%% Measures
Ws = {W, W_er, W_cm};
Wss = {'W','W_er','W_cm'};

for k = 1:length(Ws)
%     figure
%     spy(Ws{k})
    
    % Find shortest paths between all nodes
    DIST = graphallshortestpaths(Ws{k},'Directed',true);
    % Non-infinite distances
    D_noinf = DIST(DIST ~= Inf);
    % Diameter
    diam = max(D_noinf);
    % Avg distance
    adist = sum(D_noinf)/nodes^2;
    
    fprintf('%8s \n', Wss{k})
    fprintf('Diameter: %1.5f , Avg.dist: %1.5f \n', diam, adist)
    
    % Degree distributions
    win = sum(Ws{k},1);
    indist = zeros(max(win),1);
    for m = 1:max(win)
        indist(m) = sum(ismember(win,m));
    end
    indist = indist/nodes;
    
    figure
    plot(1:max(win),indist)
    
    wout = sum(Ws{k},2);
    outdist = zeros(max(wout),1);
    for m = 1:max(wout)
        outdist(m) = sum(ismember(wout,m));
    end
    outdist = outdist/nodes;
    figure
    plot(1:max(wout),outdist)
    
end










