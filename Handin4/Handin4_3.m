clearvars
close all

% Generate a random graph with avg degree 6 and 500 nodes
n = 500;
W = generate_graph(6,n);

% Plot the graph
figure
plot(graph(W),'Layout','force')

% Define states Susceptible, Infected, Recovered, Vaccinated
S = 0;
I = 1;
R = 2;
V = 3;

% Define infection and recovery probabilities
beta = 0.3;
rho = 0.7;

% Vaccination scheme
vacc = [0; 5; 15; 25; 35; 45; 55; 60; 60; 60; 60; 60; 60; 60; 60];

% Number of individuals to vaccinate in the beginning of each week
vacc_n = [0; diff(vacc)]/100*n;

% Select patient zero's. This is done by randomly selecting 10 individuals
pz = randperm(n,10);

% Number of epidemics to simulate, and how many weeks each
n_iter = 100;
n_epi = 15;

% Variables for mean values
m_susc = zeros(1,n_epi);
m_infc = zeros(1,n_epi);
m_rec = zeros(1,n_epi);
m_vacc = zeros(1,n_epi);

% Simulate 100 epidemics
for k = 1:n_iter
    
    % Initialize epidemic (pz's have state I, others S)
    X = zeros(n,n_epi);
    X(pz,1) = I;
    
    % Simulate 15 week epidemic
    for t = 1:n_epi-1
        
        % Perform vaccination according to scheme. Individuals are selected
        % randomly from the population that has not yet been vaccinated.
        % Vaccinated individuals are no longer in susceptible, infected or 
        % recovered state but instead in a vaccinated state which they
        % cannot leave
        novacc = find(X(:,t) ~= V);
        vaccinate = randperm(length(novacc),vacc_n(t));
        X(novacc(vaccinate),t) = V;
        
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
    m_vacc = m_vacc + sum(X==V,1);
end

% Calculate mean values by dividing by the number of simulations
m_susc = m_susc./n_iter;
m_infc = m_infc./n_iter;
m_rec = m_rec./n_iter;
m_vacc = m_vacc./n_iter;

% Line color definition
c_S = [167,0,255]./255;
c_I = [255,39,0]./255;
c_R = [0,189,0]./255;
c_V = [0,214,255]./255;

% Plot avg number of newly infected and newly vaccinated each week
figure
hold on
plot(1:n_epi,[10 -diff(m_susc)],'Color',c_I)
plot(1:n_epi,[0 diff(m_vacc)],'Color',c_V)
legend('Infected','Vaccinated')
xlabel('Week')
ylabel('Number of individuals')
xlim([1 n_epi])
xticks(1:n_epi)

% Plot avg total number of susceptible, infected and recovered each week
figure
hold on
plot(1:n_epi,m_susc,'Color',c_S)
plot(1:n_epi,m_infc,'Color',c_I)
plot(1:n_epi,m_rec,'Color',c_R)
plot(1:n_epi,m_vacc,'Color',c_V)
legend('Susceptible','Infected','Recovered','Vaccinated','Location','East')
xlabel('Week')
ylabel('Number of individuals')
xlim([1 n_epi])
xticks(1:n_epi)