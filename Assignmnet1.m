%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
rng(5);
%setting up parameters

beta = 0.96;
gamma = 1.3;
r = 0.04;    
rho = 0.9;    
sigma = 0.04;

%Setting values for rows and col:

row_dim = 30;
col_dim = 5;

%setting up V0 matrix

V0 = zeros(row_dim,col_dim);


%Doing Tauchen function

[y, P] = Tauchen(col_dim, 0, rho, sigma, 2);

%setting up a_min and a_max

a_min = min(y(1))/r;
%a_min = (-1)*y(1)/r;
a_max = 10;

a_grid = linspace(a_min, a_max, row_dim);

%Other matrices

V1 = ones(row_dim,col_dim);
policy = zeros(row_dim,col_dim);
consumption = zeros(row_dim,col_dim);

%Setting up tolerance:

e = 1e-9;
diff = 1;
it = 0;
max_it = 10000;

% VAlue function:

while diff > e %&& it < max_it
    it = it + 1;

    for j = 1:col_dim
        EV = V0 * P(j,:)';    % 30x1

        for i = 1:row_dim
            c = a_grid(i) + exp(y(j)) - a_grid/(1+r);

            u_c = -Inf(size(c));
            pos_cons = c > 0;
            u_c(pos_cons) = u(c(pos_cons), gamma);

            value = u_c + beta * EV';   % 1x30

            [V1(i,j), best_a] = max(value);
            policy(i,j) = a_grid(best_a);
            consumption(i,j) = c(best_a);
        end
    end
diff = max(max(abs(V1 -V0)));
V0=V1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Graphs for the Valu Function for each a
figure;
hold on;

colors = lines(length(y));  % different color for each income state

for j = 1:length(y)
    plot(a_grid, V1(:, j), 'Color', colors(j,:));
end

xlabel('Assets');
ylabel('Value function');
title('Value Function for all income');

grid on;
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Part C

T = 1000;
forget = 500;  % drop 500

eps = randn(T,1);   % normal shocks

sim_y = zeros(T,1);         
sim_a = zeros(T,1);         
sim_c = zeros(T,1);         

% Income AR(1)
for t = 2:T
    sim_y(t) = rho * sim_y(t-1) + eps(t);
end

% Drop first 500
sim_y = sim_y(forget+1:T);
sim_a = sim_a(forget+1:T);
sim_c = sim_c(forget+1:T);

T_sim = length(sim_y);

% Extracting y difference
for t = 1:T_sim
    [y_period, y_dif(t)] = min(abs(sim_y(t) - y));
end

% Simulate assets and consumption using policy functions
for t = 1:T_sim-1
    % Assets
    [a_period, a_dif] = min(abs(sim_a(t) - a_grid));
    
    
    % Apply policy function
    sim_a(t+1) = policy(a_dif, y_dif(t));
    sim_c(t) = consumption(a_dif, y_dif(t));
end

%%% Heatmap
data = [sim_y, sim_a, sim_c];

figure;
tiledlayout(3,1);
nexttile; plot (sim_y'); colorbar; title('Y');
nexttile; plot (sim_a'); colorbar; title('A');
nexttile; plot (sim_c'); colorbar; title('C');


%figure;
%h = heatmap(data');
%h.XLabel = 'Time';
%h.YLabel = 'Variable';
%h.YDisplayLabels = {'y', 'a''', 'c'};
%h.Title = 'Simulated y, a'' and c';
%h.Colormap = parula;   % optional: set colormap


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Part D

%Calculating standard deviation
std_c = std(sim_c);  %Matlab Function

%SD Code:
N = length(sim_c);           % number of observations
mean_c = sum(sim_c) / N;     % mean

sq_diff = (sim_c - mean_c).^2;

variance_c = sum(sq_diff) / (N - 1);

std_c_calc = variance_c^(1/2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Part D responses

%a) Increases but really low (there may be no change): from 0.1843 to 0.1844

%b) No change: from 0.1843 to 0.1843

%c) Increases: from 0.1843 to 0.2342

%d) Increases: from 0.1843 to 0.1863

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Utility function:

function u_func = u(c, gamma)
    if gamma == 1
        u_func = log(c);
    else
        u_func = (c.^(1 - gamma)) ./ (1 - gamma);
    end
end
