clearvars
close all

% Properties of the final graph
n = 100;
k = 2;

% Initial graph (complete graph with k nodes)
k0 = k;
W = ones(k) - diag(ones(k,1));
W = sparse(W);

for j = (k0+1):n
    
    if mod(j,2)==0 % If j even number
        c = floor(k/2); % Degree of new node
    else % If j odd number
        c = ceil(k/2); % Degree of new node
    end
    
    % Out-degree vector
    w = sum(W,2);
    
    % Probability of adding links
    p = w./sum(w);
    
    % Select c neighbors
    for m = 1:c
        % Select neighbor, and set prob to 0 for future selections
        neigh = randsample(k0,1,true,full(p));
        p(neigh) = 0;
        % Add link (both directions)
        W(k0+1,neigh) = 1;
        W(neigh,k0+1) = 1;
    end
    k0 = k0 + 1;
end

G = graph(W);
plot(G,'Layout','force')





