% Requires must installation
% https://www.biomecardio.com/MUST/index.html
% Run this script inside the must directory

transducer_name = 'L11-5V';
transducer_params = getparam(transducer_name);

% Define scatterers
scatterers.x = [ 0 0 ].* 1e-2; % x coords of each of two scatterers
scatterers.z = [ 1 2 ].* 1e-2; % z coords of each of two scatterers
scatterers.RC = [ 1 1 ]; % reflection coefficients of each of two scatterers

% Define some additional properties
transmit_delays = zeros(1, param.Nelements); % transmit delays - all zero == plane wave
transducer_params.fs = 4*transducer_params.fc; % sampling frequency (choose 4x tx frequency for great glory)
opt.ElementSplitting = 1; % to make simulations faster

% Generate RF simulation
RF = simus(scatterers.x, scatterers.z, scatterers.RC, transmit_delays, transducer_params, opt);

% Export to disk as json
export.transducer_name = transducer_name;
export.transducer_params = transducer_params;
export.scatterers = scatterers;
export.transmit_delays = transmit_delays;
export.RF = RF;

json_export = jsonencode(export);

fid = fopen(transducer_name + ".json", 'w');
fprintf(fid, '%s', json_export);
fclose(fid);

