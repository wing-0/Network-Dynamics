%% Part 1 - Epidemic on a symmetric k-regular graph
clearvars
close all

% Generate k-regular graph with n=500, k=4
n = 500;
W = zeros(n);
W = W + diag(ones(n-1,1),1); % add ones on the +1 off-diagonal
W = W + diag(ones(n-1,1),-1); % add ones on the -1 off-diagonal
W = W + diag(ones(n-2,1),2); % add ones on the +2 off-diagonal
W = W + diag(ones(n-2,1),-2); % add ones on the -2 off-diagonal
W = W + diag(ones(1,1),n-1); % add ones on the +n-1 off-diagonal
W = W + diag(ones(1,1),1-n); % add ones on the -n+1 off-diagonal
W = W + diag(ones(2,1),n-2); % add ones on the +n-2 off-diagonal
W = W + diag(ones(2,1),2-n); % add ones on the -n+2 off-diagonal
W = sparse(W); % transform it into a sparse matrix

% Define states Susceptible, Infected, Recovered
S = 0;
I = 1;
R = 2;

% Define infection and recovery probabilities
beta = 0.3;
rho = 0.7;

% Select patient zero's. This is done by randomly selecting 10 individuals
pz = randperm(n,10);

% Number of epidemics to simulate, and how many weeks each
n_iter = 100;
n_epi = 15;

% Variables for mean values
m_susc = zeros(1,n_epi);
m_infc = zeros(1,n_epi);
m_rec = zeros(1,n_epi);
% Simulate 100 epidemics
for k = 1:n_iter
    
    % Initialize epidemic (pz's have state I, others S)
    X = zeros(n,n_epi);
    X(pz,1) = I;
    
    % Simulate 15 week epidemic
    for t = 1:n_epi-1
        
        % Carry over from previous timestep
        X(:,t+1) = X(:,t);
        
        % Logical matrices for susceptible and infected individuals
        % (previous timestep)
        susc = X(:,t) == S;
        infc = X(:,t) == I;
        
        % Number of infected neighbors for all individuals
        infc_neigh = sum(W(:,infc),2);
        
        % Probability of infection for all individuals
        prob = 1-(1-beta).^infc_neigh;
        
        % If a random number if smaller than the infection probability for
        % a susceptible individual, change its state from S to I
        change = rand(n,1) < prob;
        X(change & susc,t+1) = I;
        
        % If a random number is smaller than the recovery probability for
        % an infected individual, change its state from I to R. Individuals
        % infected in this timestep are not considered
        change = rand(n,1) < rho;
        X(change & infc,t+1) = R;
    end
    
    % Add to the mean values (which are still just total values)
    m_susc = m_susc + sum(X==S,1);
    m_infc = m_infc + sum(X==I,1);
    m_rec = m_rec + sum(X==R,1);
end

% Calculate mean values by dividing by the number of simulations
m_susc = m_susc./n_iter;
m_infc = m_infc./n_iter;
m_rec = m_rec./n_iter;

% Plot avg number of newly infected each week
figure
plot(1:n_epi,[10 -diff(m_susc)])
xlabel('Week')
ylabel('Number of individuals')
xlim([1 n_epi])
xticks(1:n_epi)

% Plot avg total number of susceptible, infected and recovered each week
figure
hold on
plot(1:n_epi,m_susc,'b')
plot(1:n_epi,m_infc,'r')
plot(1:n_epi,m_rec,'c')
legend('Susceptible','Infected','Recovered','Location','East')
xlabel('Week')
ylabel('Number of individuals')
xlim([1 n_epi])
xticks(1:n_epi)

%% Part 2 - Generate random graph
clearvars
close all

k = 3;







