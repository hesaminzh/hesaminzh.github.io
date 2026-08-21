nMonte = 10;
channels = cell(nMonte, 1);

% Fill channels with random 2x2 matrices
for i = 1:nMonte
    channels{i} = randn(2,2);
end




parpool('local', 2);
wsr_all_direct = zeros(nMonte, 1);
parfor i = 1:nMonte
    H = channels{i};
    wsr_all_direct(i) = norm(H, 'fro');  % Just compute Frobenius norm
end

channels_const = parallel.pool.Constant(channels);  % One copy shared across workers

wsr_all_const = zeros(nMonte, 1);
parfor i = 1:nMonte
    H = channels_const.Value{i};
    wsr_all_const(i) = norm(H, 'fro');
end

delete(gcp('nocreate'));

isequal(wsr_all_const, wsr_all_direct)  % Should return true

max(abs(wsr_all_const - wsr_all_direct)) < 1e-12  % Should return true
