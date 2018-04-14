%% Initialize
%--------------------------------------------------------------------------

clearvars
close all

% Load data from file
load IOdownload.mat

% Select data from sets 'can2000', swe2000' and store in a cell.
% The index of the name in 'sets' is the same as for the weight
% matrix in 'W'
sets = econ([6 34]);
W = {getfield(io,sets{1}), getfield(io,sets{2})};

%--------------------------------------------------------------------------
%% a) Rank sectors using in-degree and out-degree centrality
%--------------------------------------------------------------------------

% Iterate over the three datasets
for k = 1:length(sets)
    
    fprintf('Dataset ''%s''\n',sets{k})
    
    % Calculate in-degree centrality and extract the indices for the
    % three most central sectors
    ic = sum(W{k},1);
    [~,max3] = maxk(ic,3);
    
    % Print results
    fprintf('\tMost central sectors (in-degree)\n')
    fprintf('\t\t%-70s %1.0f\n',name{max3(1)},ic(max3(1)))
    fprintf('\t\t%-70s %1.0f\n',name{max3(2)},ic(max3(2)))
    fprintf('\t\t%-70s %1.0f\n',name{max3(3)},ic(max3(3)))
    
    % Calculate out-degree centrality and extract the indices for the
    % three most central sectors
    oc = sum(W{k},2);
    [~,max3] = maxk(oc,3);
    
    % Print results
    fprintf('\tMost central sectors (out-degree)\n')
    fprintf('\t\t%-70s %1.0f\n',name{max3(1)},oc(max3(1)))
    fprintf('\t\t%-70s %1.0f\n',name{max3(2)},oc(max3(2)))
    fprintf('\t\t%-70s %1.0f\n',name{max3(3)},oc(max3(3)))
    
    fprintf('\n')
    
end

%--------------------------------------------------------------------------
%% b) Rank sectors using eigenvector centrality
%--------------------------------------------------------------------------

% Iterate over the three datasets
for k = 1:length(sets)
    fprintf('Dataset ''%s''\n',sets{k})
    
    % Calculate the dominant eigenvalue of weight matrix
    [~,D] = eig(W{k});
    lambda = max(diag(D));
    
    % Calculate eigenvector centrality as the eigenvector to lambda^(-1)*W'
    % corresponding to the eigenvalue 1. Index for eigenvalue 1 is found by
    % finding where abs(diag(D)-1) is sufficiently close to zero
    [V,D] = eig(1/lambda*W{k}');
    ec = V(:,abs(diag(D)-1)<1e-9);
    
    % Normalize centrality vector (eigenvector) according to L1 norm
    ec = ec./sum(ec);
    
    %Extract the indices for the three most central sectors
    [~,max3] = maxk(abs(ec),3);
    
    % Print results
    fprintf('\tMost central sectors (eigenvector)\n')
    fprintf('\t\t%-70s %1.3f\n',name{max3(1)},ec(max3(1)))
    fprintf('\t\t%-70s %1.3f\n',name{max3(2)},ec(max3(2)))
    fprintf('\t\t%-70s %1.3f\n',name{max3(3)},ec(max3(3)))
    
    fprintf('\n')
end

%--------------------------------------------------------------------------
%% c) Rank sectors using Katz centrality
%--------------------------------------------------------------------------

% Iterate over the three datasets
for k = 1:length(sets)
    fprintf('Dataset ''%s''\n',sets{k})
    
    % Set first values for beta and mu
    beta = 0.15;
    mu = ones(size(W{k},1),1);
    
    % Calculate the dominant eigenvalue of weight matrix
    [~,D] = eig(W{k});
    lambda = max(diag(D));
    
    % Calculate Katz centrality using the inversion described in the
    % lecture notes
    kc = (eye(size(W{k}))-1/lambda*(1-beta)*W{k}')\mu*beta;
    
    %Extract the indices for the three most central sectors
    [~,max3] = maxk(abs(kc),3);
    
    % Print results
    fprintf('\tMost central sectors (Katz w/ mu = 1 everywhere)\n')
    fprintf('\t\t%-70s %1.3f\n',name{max3(1)},kc(max3(1)))
    fprintf('\t\t%-70s %1.3f\n',name{max3(2)},kc(max3(2)))
    fprintf('\t\t%-70s %1.3f\n',name{max3(3)},kc(max3(3)))
    
    % Change mu so that it is 1 for "Wholesale & retail trade; repairs"
    % and 0 for all others
    mu = zeros(size(W{k},1),1);
    mu(contains(name,'Wholesale & retail trade; repairs')) = 1;
    
    % Calculate Katz centrality using the inversion described in the
    % lecture notes
    kc = (eye(size(W{k}))-1/lambda*(1-beta)*W{k}')\mu*beta;
    
    %Extract the indices for the three most central sectors
    [~,max3] = maxk(abs(kc),3);
    
    % Print results
    fprintf('\tMost central sectors (Katz w/ mu = 1 for Wholesale...)\n')
    fprintf('\t\t%-70s %1.3f\n',name{max3(1)},kc(max3(1)))
    fprintf('\t\t%-70s %1.3f\n',name{max3(2)},kc(max3(2)))
    fprintf('\t\t%-70s %1.3f\n',name{max3(3)},kc(max3(3)))
    
    fprintf('\n')
    
end






