%% Initialize

clearvars
close all

% Load data from file
load -ascii twitter
load -ascii users
Wa = spconvert(twitter);

% Number of nodes
n = max(size(Wa));

% Make W into a square matrix. 'users' is of size 6893x1, which indicates
% that W should be 6893x6893 instead of the given 6893x6881. By inspecting
% the 'twitter' matrix it can be found that the highest index for a tail
% node is 6893. However, the highest index for a head node is 6881. From
% this it can be concluded that the nodes with indices between 6881 and
% 6893 have no links heading towards them, and have therefore been excluded
% from W due to the sparse format. To include them and make W square as it
% should be, all zero columns should be added at the right end of W to make
% it 6893x6893.
W = sparse(n,n);
W(1:size(Wa,1),1:size(Wa,2)) = Wa;

% Calculate normalized weight matrix P
D = diag(sum(W,2));
P = D\W;

% Change all NaN's in P to 0. P has the same number of non-nan, non-zero
% entries as W, so NaN's are numerical errors and should really be zeros
P(isnan(P)) = 0;

%% a) PageRank

% Define parameters according to PageRank standard
beta = 0.15;
mu = ones(n,1);

% Initial values for y, y0
y = zeros(n,1);
yold = inf*ones(n,1);

% Iterate until y has converged sufficiently
while find(abs(y-yold) > 1e-6)
    
    % Update y according to the reformulation given by (3.24) in the
    % lecture notes. Also update yold to previous y
    yold = y;
    y = (1-beta)*P'*y + beta*mu;
    
end

%Extract the indices for the five most central sectors
[~,max5] = maxk(y,5);

fprintf('Five most central users:\n')
fprintf('%-20s %-10s\n','User ID','PageRank centrality')
fprintf('%-20.0f %-10.2f\n',[users(max5) y(max5)]')

%% b) Discrete-time consensus algorithm

close all

% Set stubborn and regular nodes.
s0 = 22;
s1 = 270;
stubborn = [s0 s1];
regular = setdiff(1:n, stubborn);

% Create matrices Q and R. P = [Q R; zeros S]
Q = P(regular,regular);
R = P(regular,stubborn);

% Set number of time steps to run for
steps = 500;

% Set opinion values for nodes. s0 = 0, s1 = 1, others 0
x = zeros(n,steps);
x(stubborn,1) = [0 1];

% Iterate for the prescribed number of time steps
for k = 2:steps
    
    % Update regular nodes according to opinion dynamics with stubborn
    % nodes using P split up into Q, R and S.
    x(regular,k) = Q*x(regular,k-1) + R*x(stubborn,k-1);
    x(stubborn,k) = x(stubborn,k-1);
end

% Split the PageRank range into three equal bins (size of each bin is
% (max(y)-min(y))/3. Randomly pick three regular nodes from each bin and 
%store their indices
ysplit = (max(y)-min(y))/3;

high = find(y > max(y)-ysplit);
high = high(randperm(length(high),3));
high = setdiff(high,stubborn);

med = find(y > min(y)+ysplit & y < max(y)-ysplit & y~=s0 & y~=s1);
med = med(randperm(length(med),3));
med = setdiff(med,stubborn);

low = find(y < min(y) + ysplit & y~=s0 & y~=s1);
low = low(randperm(length(low),3));
low = setdiff(low,stubborn);

% Plot the opinions over time. Use different lines for different bins
figure
hold on
plot(x(high,:)','-');
plot(x(med,:)','--');
plot(x(low,:)',':');
leg = legend(string([high;med;low]));
title(leg,'Node index')
title(sprintf('Opinions over time\nStubborn nodes: %4.0f (value 0), %4.0f (value 1)',...
    s0,s1))
xlabel('Timestep')
ylabel('Opinion (0-1)')

%% b) Discrete-time consensus algorithm with PageRank considerations

close all

% Split the PageRank range into three equal bins (size of each bin is
% (max(y)-min(y))/3. Randomly pick two nodes from each bin and store their
% indices. 2 nodes each so that simulations can be run for two stubborn
% nodes from the same bin
ysplit = (max(y)-min(y))/3;
high = find(y > max(y)-ysplit);
high = high(randperm(length(high),2));
med = find(y > min(y)+ysplit & y < max(y)-ysplit);
med = med(randperm(length(med),2));
low = find(y < min(y) + ysplit);
low = low(randperm(length(low),2));



% Set pairs of stubborn nodes. s0 = 0, s1 = 1
% Scenarios: combination of high PR, med PR, low PR
s0 = [high(1) high(1) high(1) med(1) med(1) med(1) low(1) low(1) low(1)];
s1 = [high(2) med(1) low(1) high(1) med(2) low(1) high(1) med(1) low(2)];


% Iterate over all PageRank scenarios
for k = 1:size(s0,2)
    
    % Set stubborn and regular nodes.
    stubborn = [s0(k) s1(k)];
    regular = setdiff(1:n, stubborn);
    
    % Create matrices Q and R. P = [Q R; zeros S]
    Q = P(regular,regular);
    R = P(regular,stubborn);
    
    % Set opinion values for nodes. s0 = 0, s1 = 1, others 0
    x = zeros(n,1);
    x(stubborn) = [0 1];
    xold = inf*ones(n,1);
    
    % Iterate until opinions have converged sufficiently
    while (sum(abs(x-xold) > 1e-6) > 0)
        
        % Update regular nodes according to opinion dynamics with stubborn
        % nodes using P split up into Q, R and S.
        xold = x;
        x(regular) = Q*xold(regular) + R*xold(stubborn);
        x(stubborn) = xold(stubborn);
    end
    
    % Plot opinion distribution as a histogram
    figure
    hist(x,10)
    ylim([0 7000])
    title(sprintf('Opinion distribution\nPageRank of stubborn nodes: %3.2f (value 0), %3.2f (value 1)',...
        y(s0(k)),y(s1(k))))
    xlabel('Opinion (0-1)')
    ylabel('Number of nodes')

end



