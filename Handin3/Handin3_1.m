%% Initialize
clearvars
close all

% Node link incidence matrix
B = load('traffic.mat','-ASCII');

% Link maximum flow capacities
Cap = load('capacities.mat','-ASCII');

% Link minimum travel times
l = load('traveltime.mat','-ASCII');

% Construct weight matrices with travel times and capacities as weight
Wl = zeros(size(B,1));
Wcap = Wl;

% Iterate over all links in B and set W('tail','head') to the
% minimum travel time/maximum flow capacity of each link
for k = 1:size(B,2)
    Wl(B(:,k)==1,B(:,k)==-1) = l(k);
    Wcap(B(:,k)==1,B(:,k)==-1) = Cap(k);
end

%% a) and b) - calculate shortest path and maximum flow between node 1-17
[DIST,~,~] = graphshortestpath(sparse(Wl),1,17);
[M,~,~] = graphmaxflow(sparse(Wcap),1,17);

fprintf('Shortest path node 1-17: %0.2f h \n',DIST)
fprintf('Maximum flow node 1-17: %0.0f vehicles/h \n',M)

%% c) - External inflow and outflow
f = load('flow.mat','-ASCII');

% Cacluate external flows, and extract positive and negative flows as in-
% and outflows
ext = B*f;
ext_in = ext(ext>0);
ext_out = -ext(ext<0);

fprintf('--------------------------------------------------------------\n')
fprintf('External inflow \n')
fprintf('\t%-6s %-8s \n','Node','Inflow')
fprintf('\t%-6.0f %-8.0f \n',[find(ext>0) ext_in]')

fprintf('External outflow \n')
fprintf('\t%-6s %-8s \n','Node','Outflow')
fprintf('\t%-6.0f %-8.0f \n',[find(ext<0) ext_out]')
fprintf('--------------------------------------------------------------\n')

%% d) - Social optimum

% Exogenous inflows
lambda = zeros(size(B,1),1);
lambda(1) = ext(1);

% Exogenous outflows
mu = zeros(size(B,1),1);
mu(end) = ext(1);

% Net exogenous flow
nu = lambda - mu;

% Optimize using cvx
cvx_begin
    variable f_so(size(B,2))
    minimize sum( l.*Cap.*inv_pos( ones(size(B,2),1) - f_so./Cap ) - l.*Cap )
    subject to
        B*f_so == nu
        f_so >= 0
cvx_end

fprintf('--------------------------------------------------------------\n')
fprintf('Social optimum flows \n')
fprintf('\t%-6s %-8s \n','Link','Flow')
fprintf('\t%-6.0f %-8.0f \n',[(1:size(B,2))' f_so]')
fprintf('--------------------------------------------------------------\n')

%% e) - Wardrop equilibrium

% Primitive of the delay function is -C*l*log(c-f)

% Optimize using cvx
cvx_begin
    variable f_uo(size(B,2))
    minimize sum( -Cap.*l.*log(Cap-f_uo) )
    subject to
        B*f_uo == nu
        f_uo >= 0
cvx_end

fprintf('--------------------------------------------------------------\n')
fprintf('\t%-6s %-10s %-10s \n','Link','SO flow', 'UO flow')
fprintf('\t%-6.0f %-10.0f %-10.0f \n',[(1:size(B,2))' f_so f_uo]')
fprintf('--------------------------------------------------------------\n')

%% f) - Tolls

% Derivative of the delay function is C*l/(c-f)^2

% Toll vector
w = f_so.*Cap.*l./(Cap-f_so).^2;

% Optimize using cvx
cvx_begin
    variable f_t(size(B,2))
    minimize sum( -Cap.*l.*log(Cap-f_t) + w.*f_t )
    subject to
        B*f_t == nu
        f_t >= 0
cvx_end

fprintf('--------------------------------------------------------------\n')
fprintf('\t%-6s %-10s %-10s %-10s \n','Link','SO flow', 'UO flow', 'UO w/tolls')
fprintf('\t%-6.0f %-10.0f %-10.0f %-10.0f \n',[(1:size(B,2))' f_so f_uo f_t]')
fprintf('--------------------------------------------------------------\n')

%% g) - Total additional delay






