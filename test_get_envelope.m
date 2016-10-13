clear all; close all; clc;
%% Validation file :
% Here I show an example of how to get the "decaying exponent" of a damping oscilating
% signal. 
% I use the hilbert-function to get the envelope of the signal, and save it
% on the variable "env". 
%However, it has issues at the begging of the signal, make a plot for
%yourself using the comments.
% My solution was 
% -> To duplicate the signal, thus I have a full modulation envelope
% -> apply the  hilbert-function
% -> then slice the result in the middle

%t = 0:1e-4:1;
%x = ((exp(-(t./0.55).^2)).*cos(2*pi*50*t));

%duplication
tt = 0:1e-4:1;
xx = ((exp(-(tt./0.15).^2)).*cos(2*pi*50*tt)); % "0.15" is the value that I look for

t = [-tt(end:-1:1) tt(2:end)];
x = [xx(end:-1:1) xx(2:end)];

% get the envelope
y = hilbert(x);
env = abs(y);

% slice in the middle
half_life = env(round(end/2)-1:end);
time_span0 = t(round(end/2)-1:end);


remover = find(half_life < half_life(1)/2);
time_span0(remover) = [];
half_life(remover) = [];

% remover = find(env < env(1)/2);
% time_span0 = t(remover);
% half_life = env(remover);


figure(1); 
plot(time_span0, half_life);
hold on
plot(t,x)
xlim([0 0.1])

%% make the fit
% Here there is something to say:
% There are some values VERY close to zero, however they are negative. This
% negative values make complex number after the square root, and generate
% the Warning on the prompt.

% They are small enough to be ignored on the computations, otherwise I will
% have to create more source code, and will lower the code performance

[c, gof] = fit(time_span0',sqrt(-log(half_life/half_life(1)))','poly1');
fprintf('decaimento = %.2f \n', 1/c.p1)

figure(2); 
plot(time_span0, sqrt(-log(half_life/half_life(1))));
title(gof.rsquare)





