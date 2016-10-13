clear all;clc; close all;
%% Validation file :
% Here, again the goal is get the envelope of an damping signal
% Now, I identify the peaks of the signal, interpolate them using spline
% The remaining code is equal to 'test_get_envelope.m'

%% Get some original data
%load spins.mat
%load tempo.mat
%plot(time_span,a,'*')

tt = 0:1e-4:1;
xx = ((exp(-(tt./0.15).^2)).*cos(2*pi*50*tt)); % "0.15" is the value that I look for

t = [-tt(end:-1:1) tt(2:end)];
x = [xx(end:-1:1) xx(2:end)];

% Duplicate the signal, to make easy to identify the initial value such as
% a peak
%t = [-time_span(end:-1:1) time_span(2:end)];
%x = abs([a(end:-1:1) a(2:end)]);

figure(1);
findpeaks(x,t);
[pks,loc] = findpeaks(x,t)
% Eliminate the valleys - they are peaks, but not the correct one
remover = find(pks < 0.35); % 0.35 could be changed
pks(remover) = [];
loc(remover) = [];

% Create the interpolation
tq = min(loc):0.001:max(loc);
interpolado = interp1(loc,pks,tq,'spline');
figure(2);
plot(t,x,'o',tq,interpolado,':.','Linewidth',2);

%% The follow code is a copy from test_get_envelope.m
half_life = interpolado(round(end/2)-1:end);
time_span0 = tq(round(end/2)-1:end);

remover = find(half_life < half_life(1)/2);
time_span0(remover) = [];
half_life(remover) = [];

figure(3); 
plot(time_span0, half_life);
hold on
plot(t,x)
xlim([0 0.15])

%% make the fit
% Here there is something to say:
% There are some values VERY close to zero, however they are negative. This
% negative values make complex number after the square root (and generate
% the Warning on the prompt). Then I must to eliminate them

remover = find(-log(half_life/half_life(1)) < 0);
time_span0(remover) = [];
half_life(remover) = [];


[c, gof] = fit(time_span0',sqrt(-log(half_life/half_life(1)))','poly1');
fprintf('decaimento = %.2f \n', 1/c.p1)

figure(4); 
plot(time_span0, sqrt(-log(half_life/half_life(1))));
title(gof.rsquare)
