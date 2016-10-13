function [ decaimento ] = damped_oscilations_coeff( time_span,spin )
% For more information look test_get_envelop_.m
% Here I only copy and paste the source code that have been worked on
% previous tests

%%  
t = [-time_span(end:-1:1) time_span(2:end)];
x = abs([spin(end:-1:1) spin(2:end)]);
[pks,loc] = findpeaks(x,t);
% Eliminate the valleys - they are peaks, but not the correct one
remover = find(pks < 0.35); % 0.35 could be changed
pks(remover) = [];
loc(remover) = [];

% Create the interpolation
tq = min(loc):0.001:max(loc);
interpolado = interp1(loc,pks,tq,'spline');
%% The follow code is a copy from test_get_envelope.m
half_life = interpolado(round(end/2)-1:end);
time_span0 = tq(round(end/2)-1:end);

remover = find(half_life < half_life(1)/2);
time_span0(remover) = [];
half_life(remover) = [];

remover = find(-log(half_life/half_life(1)) < 0);
time_span0(remover) = [];
half_life(remover) = [];
%% For debug
% plot(time_span0, half_life);
% hold on
% plot(t,x)
%hold off
%xlim([0 1])
%% The most important informatio: The fit
[c, ~] = fit(time_span0',sqrt(-log(half_life/half_life(1)))','poly1');
decaimento = 1/c.p1;

end

