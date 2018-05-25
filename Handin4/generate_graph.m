function W = generate_graph(k,n)
% GENERATE_GRAPH Generate a random graph
%   W = generate_graph(k,n) uses preferential attachment to generate a
%   random graph with n nodes and average degree k. W is a sparse weight
%   matrix

% Initial graph (complete graph with k0 nodes)
k0 = k + 1;
W = ones(k0) - diag(ones(k0,1));
W = sparse(W);

% Add nodes from k0+1 to n
for m = k0+1:n
    
    % Set degree c of the new node. To ensure avg degree of k even for odd
    % k, c is set to floor or ceil of k/2 every other iteration
    if mod(m,2)==0
        c = floor(k/2);
    else
        c = ceil(k/2);
    end
    
    % Out-degree vector
    w = sum(W,2);
    
    % Probability of adding links
    p = w./sum(w);
    
    % Select c neighbors
    for j = 1:c
        
        % Select neighbor, and remove that neighbor from the population to
        % choose from
        neigh = randsample(k0,1,true,full(p));
        p(neigh) = 0;
        % Add link (both directions)
        W(k0+1,neigh) = 1;
        W(neigh,k0+1) = 1;
    end
    
    % Increase k0
    k0 = k0 + 1;
end
end

