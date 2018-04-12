clearvars
close all

% Weight matrix for the intermarriages.
W = [0 1 1 1 0 0 0 0 0 0 0 0 0 0 0; % 1. Castellani
	1 0 1 0 0 0 1 0 0 0 0 0 0 0 0; % 2. Peruzzi
	1 1 0 0 1 0 1 0 0 0 0 0 0 0 0; % 3. Strozzi
	1 0 0 0 0 0 0 0 1 0 0 0 0 0 0; % 4. Barbadori
	0 0 1 0 0 0 0 0 1 0 0 0 0 0 0; % 5. Ridolfi
	0 0 0 0 0 0 0 0 1 0 0 0 0 0 0; % 6. Acciaiuoli
	0 1 1 0 0 0 0 0 0 1 0 0 0 0 0; % 7. Bischeri
	0 0 0 0 0 0 0 0 1 1 0 0 0 0 0; % 8. Tornabuoni
	0 0 0 1 1 1 0 1 0 0 1 1 0 0 0; % 9. Medici
	0 0 0 0 0 0 1 1 0 0 1 0 1 0 0; % 10. Guadagni
	0 0 0 0 0 0 0 0 1 1 0 0 0 1 0; % 11. Albizzi
	0 0 0 0 0 0 0 0 1 0 0 0 0 0 1; % 12. Salviati
	0 0 0 0 0 0 0 0 0 1 0 0 0 0 0; % 13. Lamberteschi
	0 0 0 0 0 0 0 0 0 0 1 0 0 0 0; % 14. Ginori
	0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]; % 15. Pazzi

% Calculating the invariant probability distribution P = D\W = diag(w)\W.
% w is the out-degree of each node, so the sum of each row
P = diag(sum(W,2))\W;

%% Calculating invariant distribution centrality
% This is equivalent to finding the eigenevctor to P' corresponding to the
% eigenvalue 1 (equation is pi = P'*pi)
[V,D] = eig(P');

% Find index of the eigenvalue 1
J = find(abs(diag(D)-1)<1e-6);

% Find eigenvector and normalize according to L1 norm (taxi norm)
pi = V(:,J);
pi = pi/sum(pi);

% Centrality for Medici (9), Strozzi (3), Tornabuoni (8)
fprintf('Medici: %0.3f\n', pi(9))
fprintf('Strozzi: %0.3f\n', pi(3))
fprintf('Tornabuoni: %0.3f\n', pi(8))

%% Calculating closeness centrality

% Defining distance matrix as NaN with 0 on diag (self-distance = 0)
D = NaN(15);
D(1:16:end) = 0;

Wa = W;
len = 1;

% Iterates over distance (by multiplying W together with itself). Sets the
% distance to D if the weight matrix is nonzero and has not been changed
% from NaN yet
while(sum(sum(isnan(D))) > 0)
    D(Wa>0 & isnan(D)) = len;
    len = len + 1;
    Wa = Wa*W;
end

% Calculates closeness centrality
pi = 15./sum(D,2);

% Centrality for Medici (9), Strozzi (3), Tornabuoni (8)
fprintf('Medici: %0.3f\n', pi(9))
fprintf('Strozzi: %0.3f\n', pi(3))
fprintf('Tornabuoni: %0.3f\n', pi(8))

%% Calculating decay centrality

% Defining distance matrix as NaN with 0 on diag (self-distance = 0)
D = NaN(15);
D(1:16:end) = 0;

Wa = W;
len = 1;

% Iterates over distance (by multiplying W together with itself). Sets the
% distance to D if the weight matrix is nonzero and has not been changed
% from NaN yet
while sum(sum(isnan(D))) > 0
    D(Wa>0 & isnan(D)) = len;
    len = len + 1;
    Wa = Wa*W;
end

delta = [0.25 0.5 0.75];

for d = delta
    
    % Calculates centrality for all nodes
    pi = sum(d.^D,2) - ones(15,1);
    
    fprintf('delta = %0.2f\n',d)
    
    % Centrality for Medici (9), Strozzi (3), Tornabuoni (8)
    fprintf('Medici: %0.3f\n', pi(9))
    fprintf('Strozzi: %0.3f\n', pi(3))
    fprintf('Tornabuoni: %0.3f\n\n', pi(8))
end

%% Calculating PageRank centrality

beta = 0.15;
mu = ones(15,1);

term = beta*mu;
pi = term;
k = 1;
while term > 1e-6
    term = beta*(1-beta)^k*(P')^k*mu;
    pi = pi + term;
    k = k + 1;
end

% Centrality for Medici (9), Strozzi (3), Tornabuoni (8)
fprintf('Medici: %0.3f\n', pi(9))
fprintf('Strozzi: %0.3f\n', pi(3))
fprintf('Tornabuoni: %0.3f\n\n', pi(8))






