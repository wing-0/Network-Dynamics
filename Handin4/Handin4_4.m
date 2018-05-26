%% Initial definitions
clearvars
close all

% Line color definition
c_S = [167,0,255]./255;
c_I = [255,39,0]./255;
c_R = [0,189,0]./255;
c_V = [0,214,255]./255;

% Population (scaled)
n = 934;

% Week numbers
week = [42:52 1:5]';

% Reference number for infected individuals each week (scaled)
I0 = [1; 1; 3; 5; 9; 17; 32; 32; 17; 5; 2; 1; 0; 0; 0; 0];

% Define states Susceptible, Infected, Recovered, Vaccinated
S = 0;
I = 1;
R = 2;
V = 3;

% Vaccination scheme
vacc = [5; 9; 16; 24; 32; 40; 47; 54; 59; 60; 60; 60; 60; 60; 60; 60];

% Number of individuals to vaccinate in the beginning of each week
vacc_n = [0; diff(vacc)]/100*n;

% Initial values and deltas for parameters
k0 = 10;
d_k = 1;
beta0 = 0.3;
d_beta = 0.1;
rho0 = 0.6;
d_rho = 0.1;

% Values for parameters in previous iteration. NaN here since no previous
% iteration has been done
k0_old = nan;
beta0_old = nan;
rho0_old = nan;

%% Simulation

% Keeping track of the number of iterations
iter = 0;

% The number of times the algorithm has decreased the values in the deltas.
% This happens each time the algorithm cannot find a better set of
% parameters
refined = 0;

% Loop until the delta values have been refined twice and a solution has
% been found
while refined < 3
    
    % Update values for parameters in previous iteration
    k0_old = k0;
    beta0_old = beta0;
    rho0_old = rho0;
    
    % Parameter spaces for all parameters. Boundaries for what values they 
    % can take on are included such that:
    % 2 <= k <= n-1
    % 0 <= beta <= 1
    % 0 <= rho <= 1
    k = [2 + (k0-d_k > 0)*(k0-d_k - 2); 
         k0;
         n-1 + (k0+d_k < n)*(k0+d_k - n+1)];
    beta = [(beta0-d_beta > 0)*(beta0-d_beta);
            beta0;
            (beta0+d_beta < 1)*(beta0+d_beta)];
    rho = [(rho0-d_rho > 0)*(rho0-d_rho);
           rho0;
           (rho0+d_rho < 1)*(rho0+d_rho)];
    
    % Parameter vectors for all combinations of parameters in the
    % parameter space. In total there are 27 combinations
    k = [k(1)*ones(9,1); k(2)*ones(9,1); k(3)*ones(9,1)];
    beta = repmat([beta(1)*ones(3,1); beta(2)*ones(3,1); beta(3)*ones(3,1)],3,1);
    rho = repmat(rho,9,1);
    
    % Number of epidemics to simulate, and how many weeks each
    n_iter = 10;
    n_epi = 15;
    
    % Mean values for all parameter simulations
    m_susc = zeros(length(k),n_epi);
    m_infc = zeros(length(k),n_epi);
    m_rec = zeros(length(k),n_epi);
    m_vacc = zeros(length(k),n_epi);
    
    % Newly infected each week (mean)
    n_infc = zeros(length(k),n_epi);
    
    % Iterate over possible parameters
    for m = 1:length(k)
        
        % Generate a random graph with avg degree k(m) and n nodes
        W = generate_graph(k(m),n);
        
        % Select patient zero's. This is done by randomly selecting I0(1)
        % individuals
        pz = randperm(n,I0(1));
        
        % Simulate 10 epidemics
        for j = 1:n_iter
            
            % Initialize epidemic (pz's have state I, others S)
            X = zeros(n,n_epi);
            X(pz,1) = I;
            
            % Add the first infected to mean
            n_infc(m,1) = n_infc(m,1) + length(pz);
            
            % Initial vaccination from vacc(1). Random selection, but make
            % sure that the initial infected are not vaccinated at this
            % time (that would lead to the epidemic not happening at all)
            noinfc = find(X(:,1) ~= I);
            vaccinate = randperm(length(noinfc),round(vacc(1)*n/100));
            X(noinfc(vaccinate),1) = V;
            
            % Simulate 15 week epidemic
            for t = 1:n_epi-1
                
                % Perform vaccination according to scheme. Individuals are
                % selected randomly from the population that has not yet
                % been vaccinated. Vaccinated individuals are no longer in
                % susceptible, infected or recovered state but instead in a
                % vaccinated state which they cannot leave
                novacc = find(X(:,t) ~= V);
                vaccinate = randperm(length(novacc),round(vacc_n(t)));
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
                prob = 1-(1-beta(m)).^infc_neigh;
                
                % If a random number if smaller than the infection
                % probability for a susceptible individual, change its
                % state from S to I
                change = rand(n,1) < prob;
                X(change & susc,t+1) = I;
                
                % Add to the mean value of newly infected
                n_infc(m,t+1) = n_infc(m,t+1) + sum(change & susc);
                
                % If a random number is smaller than the recovery
                % probability for an infected individual, change its state
                % from I to R. Individuals infected in this timestep are
                % not considered
                change = rand(n,1) < rho(m);
                X(change & infc,t+1) = R;
            end
            
            % Add to the mean values (which are still just total values)
            m_susc(m,:) = m_susc(m,:) + sum(X==S,1);
            m_infc(m,:) = m_infc(m,:) + sum(X==I,1);
            m_rec(m,:) = m_rec(m,:) + sum(X==R,1);
            m_vacc(m,:) = m_vacc(m,:) + sum(X==V,1);
        end
        
        % Calculate mean values by dividing by the number of simulations
        m_susc(m,:) = m_susc(m,:)./n_iter;
        m_infc(m,:) = m_infc(m,:)./n_iter;
        m_rec(m,:) = m_rec(m,:)./n_iter;
        m_vacc(m,:) = m_vacc(m,:)./n_iter;
        n_infc(m,:) = n_infc(m,:)./n_iter;
    end
    
    % Reference for newly infected per week
    ref_infc = ones(length(k),1)*I0(2:end)';
    
    % Root mean square error between the results of each parameter
    % combination and the real epidemic
    RMSE = sqrt(sum((n_infc - ref_infc).^2, 2)/15);
    
    % Update k0, beta0 and rho0 to the combination resulting in minimum
    % RMSE
    k0 = k(RMSE==min(RMSE));
    beta0 = beta(RMSE==min(RMSE));
    rho0 = rho(RMSE==min(RMSE));
    
    iter = iter+1;
    fprintf('%1.0f iterations complete, minimum RMSE = %1.2f \n',...
        iter,min(RMSE))
    
    % If k0, beta0 and rho0 are the same as in the previous iteration (i.e
    % when the algorithm is about to complete) the size of deltas is
    % decreased so that the algorithm can find a more accurate solution
    % Not done for d_k since it needs to be an integer
    if(k0 == k0_old && beta0 == beta0_old && rho0 == rho0_old)
        d_beta = d_beta/2;
        d_rho = d_rho/2;
        refined = refined + 1;
        fprintf('Refined! #%1.0f \n',refined)
    end
    
    % Only save the mean values that are related to the minimum RMSE
    % solution
    m_susc = m_susc(RMSE==min(RMSE),:);
    m_infc = m_infc(RMSE==min(RMSE),:);
    m_rec = m_rec(RMSE==min(RMSE),:);
    m_vacc = m_vacc(RMSE==min(RMSE),:);
    n_infc = n_infc(RMSE==min(RMSE),:);
    
    hold off
    plot(1:n_epi,I0(2:end),'b')
    hold on
    plot(1:n_epi,n_infc,'r')
    legend('True','Simulation')
    xlabel('Week')
    ylabel('Number of infected individuals')
    xlim([1 n_epi])
    xticks(1:n_epi)
    xticklabels(week(2:end))
    pause(0.1)
end

% Plot avg number of newly infected (simulated) and the true value each
% week
figure
hold on
plot(1:n_epi,I0(2:end),'b')
plot(1:n_epi,n_infc,'r')
legend('True','Simulation')
xlabel('Week')
ylabel('Number of infected individuals')
xlim([1 n_epi])
xticks(1:n_epi)
xticklabels(week(2:end))

% Plot avg total number of susceptible, infected, recovered and vaccinated
% each week
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













