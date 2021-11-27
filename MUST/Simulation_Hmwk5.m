% Requires must installation
% https://www.biomecardio.com/MUST/index.html
% Run this script inside the must directory
%% Make movie
param = getparam('L11-5V');

% Define scatterers
xs = [ 0 0 ].* 1e-2;
zs = [ 1 2 ].* 1e-2;
RC = [ 1 1 ];

txdel = zeros(1,128); % transmit delays
param.fs = 4*param.fc; % sampling frequency
opt.ElementSplitting = 1; % to make simulations faster

param.movie = [5 3];
[F,info] = mkmovie(xs,zs,RC,txdel,param,'l11-5v_2.gif');

figure;colormap([1-hot(128); hot(128)]);
framenum = round(linspace(1,size(F,3),9));
file = 'L11-5v_.gif'
for k = 1:1000
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

%% Beamform
%% Generate RF
RF = simus(xs,zs,RC,txdel,param,opt);

xc = ((0:127)-63-0.5).*param.pitch;
fs = param.fs;
c = 1540;
dz = c/(2*fs);

for ln = 1:128 % reconstruct a line per element
  x_ca = xc(ln); % x-axis (along transducer face) location of reconstructed line
  rf_pad     = [RF; zeros(500,128) ]; % extra pad for shifts at bottom of image

  for id = 0:(1100-1)
      z = id .* dz;
      % Calculate num samples to shift from position of tstart
      temp(id+1,:) % = you insert code here to calculate the delay matrix
  end

  for dep = 1:1100
    shift = temp(dep,:); % take delays for depth depth
    for ele=1:128 % first_piezo:last_piezo
      rf_shift(dep,ele)   = rf_pad(dep+round(shift(ele)),ele);
      shift_see(dep,ele)  = shift(ele);
    end
  end
  bmode(:,ln) = sum(rf_shift,2);
end

%%%% figure(36);imagesc(20*log10(abs(%insert envelope signal here% )));colormap(gray);
caxis([30 78])
