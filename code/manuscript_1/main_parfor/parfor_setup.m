% Access the local cluster object
cluster = parcluster('local');

% Override the NumWorkers property
cluster.NumWorkers = 75; % Set this to your desired number (max 512)

% Save the updated profile
saveAsProfile(cluster, 'local');

% Start the pool with new setting
parpool(cluster, 75);
