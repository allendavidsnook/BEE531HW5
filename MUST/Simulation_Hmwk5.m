% Requires must installation
% https://www.biomecardio.com/MUST/index.html
% Run this script inside the must directory
%% Make movie
param = getparam('L11-5V');

% Define scatterers
xs = [ 0 0 ].* 1e-2; % x coords of each of two scatterers
zs = [ 1 2 ].* 1e-2; % z coords of each of two scatterers
RC = [ 1 1 ]; % reflection coefficients of each of two scatterers

txdel = zeros(1,128); % transmit delays - all zero == plane wave
param.fs = 4*param.fc; % sampling frequency (choose 4x tx frequency for great glory
opt.ElementSplitting = 1; % to make simulations faster

% Save time - don't generate the movie over and over again
make_movie = 0;

if make_movie == 1
    param.movie = [5 3];
    [F,info] = mkmovie(xs,zs,RC,txdel,param,'l11-5v_2.gif');

    figure;colormap([1-hot(128); hot(128)]);
    framenum = round(linspace(1,size(F,3),9));
    file = 'L11-5v_.gif'
    for k = 1:1000 % k is the frame number in the roughly 47 sec long movie this creates
        image(info.Xgrid*100,info.Zgrid*100,1.0*F(:,:,k))
        hold on
        scatter(xs*100,zs*100,10,'w','filled')
        hold off
        axis equal ij tight
        title([int2str(info.TimeStep*k*1e6) ' \mus'])
        ylabel('[cm]')
        xlabel('[cm]')
        set(gca,'box','off')
        snapnow
        frame(k) = getframe(gcf);
        %pause(.1)
    end

    vid = VideoWriter(['test_2.avi']);
    vid.FrameRate = 20;
    open(vid)
    writeVideo(vid,frame);
    close(vid);
end

%% Beamform
%% Generate RF
%% simus takes
%%% the xs and zs position of the scatterers,
%%% the RC?
%%% the transmit delays (all zeros for this plane wave)
%%% the transducer parameters (e.g. center frequency, pitch, number of
%%% elements
%%% any other options (e.g. enable element splitting)
RF = simus(xs,zs,RC,txdel,param,opt);

% This will generate a 1100r x 128c RF matrix
% Let's export it



xc = ((0:127)-63-0.5).*param.pitch; % -0.019 to +0.019 m
fs = param.fs; % 30,400,000 samp/s
c = 1540; % m/s
dz = c/(2*fs); % dz = 25 um

line_64_rf_shifted = zeros(1100, 128);

line_64_delay_matrix = zeros(1100,128);

% copy the RF data (1100x128) and pad it out to 1900x128
rf_pad     = [RF; zeros(800,128) ]; % extra pad for shifts at bottom of image

for ln = 1:128 % reconstruct a line per element
  x_ca = xc(ln); % x-axis (along transducer face) location of reconstructed line

  % we'll calculate a delays matrix (temp) for the current rx element
  % we need to calculate the correct delay for each depth (1100 rows) from each
  % possible tx element (128 columns)
  % so, temp contains 1100 rows and 128 columns
  temp = zeros(1100,128);

  % for the line we are receiving on
  % calculate the additional delay we need to shift by
  for id = 0:(1100-1)
      % z will hold the depth that we are evaluating, which is a function
      % of id
      z = id .* dz; % z will vary from 0 when id=0 to 0.0279 m when id=1100
      % Calculate num samples to shift from position of tstart
      % temp(depth, element) contains the delay for that depth
      % calculate the delay matrix
      % as a row vector with one column per transmitting element
      for tx_element_index = 1:128
          % we are listening on rx element ln at position x_ca
          % the echo from depth z below tx element tx_element_index will be
          % delayed to the rx element ln by an amount we need to calculate

          % first, calculate the distance between the rx and tx elements
          x_sep_tx_rx = abs(x_ca - xc(tx_element_index));
          % then calculate the total distance from that depth below the tx
          % element all the way back to the rx element
          total_distance = sqrt(x_sep_tx_rx * x_sep_tx_rx + z * z);
          % subtract out the z depth so we just get the delta we need
          delta_distance = total_distance - z;
          % discretize it to index into the pad correctly
          delay = delta_distance / dz / 2;
          temp(id+1,tx_element_index) = delay; % save the delta for the given depth
      end
  end

  if ln == 64
      line_64_delay_matrix = temp;
  end

  % apply the delay to the received signal
  for dep = 1:1100 % for each discrete depth
    shift = temp(dep,:); % grab the delays for this depth
    % iterate over each transmit element
    for ele=1:128 % first_piezo:last_piezo
      % calculate time shifted rf at this depth contributed by each transmitting element
      shifted_depth = dep + round(shift(ele));
      rf_shifted(dep, ele) = rf_pad(shifted_depth, ele);
      % NOT USED shift_see(dep,ele)  = shift(ele);
    end
  end

  if ln == 64
      line_64_rf_shifted = rf_shifted
  end

  % lastly, sum all the contributions for this receive line
  % sum(x, 2) returns a column vector containing the sum of each row
  RF_beamformed(:,ln) = sum(rf_shifted,2);
end

%% FIGURE X - Delay Matrix for Line 64
figure(4);
imagesc(line_64_delay_matrix);
title('Line 64 Delay Matrix');
colorbar;

%% FIGURE 5 - RF unshifted
% complex -> magnitude
env = abs(RF);
% log compress it and display it
env_log_compressed = 20*log10(env);
figure(5);
imagesc(env_log_compressed);
colormap(gray);
caxis([30 78])
title('Figure 5 - RF signal, unflattened');
xlabel("Transducer Element");
ylabeltext = sprintf('Depth / %0.2e m', dz);
ylabel(ylabeltext);

%% FIGURE 6 - RF shifted
% complex -> magnitude
line_64_env_shifted = abs(line_64_rf_shifted);
% log compress it and display it
line_64_env_shifted_log_compressed = 20*log10(line_64_env_shifted);
figure(6);
imagesc(line_64_env_shifted_log_compressed);
colormap(gray);
caxis([30 78])
title('Figure 6 - RF signal, flattened');
xlabel("Transducer Element");
ylabeltext = sprintf('Depth / %0.2e m', dz);
ylabel(ylabeltext);

%% FIGURE 7 - Line 64 Beamformed
% complex -> magnitude
line_64_env_beamformed = abs(RF_beamformed(:, 64));
% log compress it and display it
line_64_env_beamformed_log_compressed = 20*log10(line_64_env_beamformed);
figure(7);
plot(line_64_env_beamformed_log_compressed);
xlabeltext = sprintf('Depth / %0.2e m', dz);
xlabel(xlabeltext);
ylabel("Envelope, dB");
title("Figure 7 - Line 64 Bmode");

%% FIGURE 8 - Beamformed Bmode
% complex -> magnitude
env_beamformed = abs(RF_beamformed);
% log compress it and display it
env_beamformed_log_compressed = 20*log10(env_beamformed);
figure(8);
imagesc(env_beamformed_log_compressed);
colormap(gray);
caxis([30 78])
title("Figure 8 - Bmode Image");
xlabel("Transducer Element");
ylabeltext = sprintf('Depth / %0.2e m', dz);
ylabel(ylabeltext);
