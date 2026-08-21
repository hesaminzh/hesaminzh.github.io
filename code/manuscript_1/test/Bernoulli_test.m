%% Parameters
p = 0.73;     % LoS probability 0 ≤ p ≤ 1 so 1-p would be blockage probability
N = 1e5;      % number of Monte-Carlo snapshots


%% Bernoulli generation
B = rand(N,1) < p; % B \sim Bernoulli(p), p here is LoS probability
%B = (1 - rand(N,1)) > 1-p; % B \sim Bernoulli(1-p), p here is blockage probability
B = double(B);     % optional: cast to double if you need maths later

%% Sanity check
fprintf('Empirical mean ≈ %.4f (should be %.2f)\n', mean(B), p);


%%% Another example
N = 1e6;
p_avail = 0.8;               % 80% chance link is available

% (a) Direct approach: compare against p_avail
isAvail1 = (rand(N,1) < p_avail);
fprintf('Empirical avail ≈ %.3f\n', mean(isAvail1));

% (b) If you only have blockage probability q = 1-p
q_block = 1 - p_avail;
isBlocked = (rand(N,1) < q_block);
%isAvail2  = ~isBlocked;
% or
isAvail2  = rand(N,1) > q_block;
fprintf('Empirical avail ≈ %.3f\n', mean(isAvail2));
