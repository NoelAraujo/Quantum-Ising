%% Original  - 100% funcionando como esperado
clear all; clc; close all;
t = (0:0.1:3)';

N = 2500;
a = 0.25;
logP=0*t;
for i=1:N
    logP=logP+log(abs(cos(2*t/(min(i,N+1-i))^a)));
end

plot(t.^2,N^(2*a-1)*logP)



%% Crio os dados periodicos e aleatorios
clear all; clc; close all;
%t = (0:0.1:sqrt(10))';
N = 2000;
alpha = 0.25;

modo = 'random'; time_init = 0; time_end = sqrt(10); time_steps = 500; d=1; J=0;
[d_ij,~,~,~] = generate_data(modo, N, time_init,time_end,time_steps,alpha,d,J);

t = linspace(time_init,time_end,time_steps);
logP = 0*t;
for i=1:N
    %logP = logP+log(abs(cos(2*t/(min(i,N+1-i))^alpha)));
    logP = logP+log(abs(cos(2*t/(d_ij(1,i).^alpha))));
end
plot(t.^2,N^(2*alpha-1)*logP,'o-')

%% Calculo valor do spin para 1 particula
clear all; clc; close all;
%t = (0:0.1:sqrt(10))';
N = 2000;
alpha = 0.25;

modo = 'random'; time_init = 0; time_end = sqrt(10); time_steps = 500; d=1; J=0;h=0;
[d_ij,~,spin0,~] = generate_data(modo, N, time_init,time_end,time_steps,alpha,d,J);

t = linspace(time_init,time_end,time_steps);
logP = 0*t;
for j=1:N
    logP = logP + log(abs(cos(2*t/(d_ij(1,j).^alpha))));
end
plot(t.^2,N^(2*alpha-1)*logP,'o-')

spin = spin0(1).*cos(2*h*t).*exp(logP);
figure()
plot(t,spin);
xlim([0 0.1])

%% Calculo valor do spin para N particula
clear all; clc; close all;
%t = (0:0.1:sqrt(10))';
N = 2000;
alpha = 0.25;

modo = 'random';  time_init = 0; time_end = sqrt(10); time_steps = 500; d=1; J=1;h=0;
[d_ij,~,spin0,~] = generate_data(modo, N, time_init,time_end,time_steps,alpha,d,J);

t = linspace(time_init,time_end,time_steps);
for i=1:N
logP = 0*t;
for j=1:N
    logP = logP + log(abs(cos(2*t/(d_ij(i,j).^alpha))));
end
%plot(t.^2,N^(2*alpha-1)*logP,'o-')

spin = spin0(i).*cos(2*h*t).*exp(logP);
%figure()
plot(t,spin);
hold on
xlim([0 0.1])
drawnow
end

%% Calculo valor do spin para N particula com alpha > 0.5 E exibo a curva teorica
clear all; clc; close all;
N = 2000;
alpha = 0.85;

modo = 'random'; time_init = 0; time_end = sqrt(10); time_steps = 500; d=1; J=1; h=0;
[d_ij,y_t,spin0,time_span] = generate_data(modo, N, time_init,time_end,time_steps,alpha,d,J);


hold on
for i=1:N
    logP = 0*time_span;
    for j=1:N
        logP = logP + log(abs(cos(2*time_span/(d_ij(i,j).^alpha))));
    end
    spin(i,:) = spin0(i).*cos(2*h*time_span).*exp(logP);
    
    plot(time_span,abs(spin(i,:)));
    plot(time_span,y_t,'k','LineWidth',5)
    hold on
    xlim([0 1])
    drawnow
end

%% Calculo valor do spin para N particula com alpha < 0.5 exibindo a curva de validacao

clear all; clc; close all;
N = 5000;
alpha = 0.25;

modo = 'static';  time_init = 0; time_end = sqrt(10); time_steps = 500; d=1; J=1; h=0;
[d_ij,~,spin0,time_span] = generate_data(modo, N, time_init,time_end,time_steps,alpha,d,J);


hold on
for i=1:N
    logP = 0*time_span;
    for j=1:N
        logP = logP + log(abs(cos(2*time_span/(d_ij(i,j).^alpha))));
    end
    spin(i,:) = spin0(i).*cos(2*h*time_span).*exp(logP);
    plot(time_span.^2,N^(2*alpha-1)*logP)
    xlim([0 0.014])
    drawnow
    %plot(time_span,abs(spin(i,:)));
    %plot(time_span,y_t,'k','LineWidth',5)
    %hold on
    %xlim([0 1])
    %drawnow
end


%% Calculo valor do spin para N particula com alpha < 0.1
% exibindo o spin de meia vida 
clear all; clc; close all;
N = 250;
alpha = 0.25;

modo = 'static'; nps = N; time_init = 1e-5; time_end = sqrt(1); time_steps = 1000; d=1; J=1; h=0;
[d_ij,y_t,spin0,time_span] = generate_data(modo, N, time_init,time_end,time_steps,alpha,d,J);


hold on
tic;
for i=1:N
    logP = 0*time_span;
    for j=1:N
        logP = logP + log(abs(cos(2*time_span/(d_ij(i,j).^alpha))));
    end
    spin(i,:) = spin0(i).*cos(2*h*time_span).*exp(logP);
       
    % spin' = spin0*cos()*exp(logP)
    % spin'/spin0*cos() = exp(logP)
    %---> log(spin'/spin0*cos()) = logP <---
    
    %plot(time_span.^2,N^(2*alpha-1)*logP)
    
    time_span0 = time_span;
    half_life = abs(spin(i,:));
    remover = find(half_life < 0.1);
    
    time_span0(remover) = [];
    half_life(remover) = [];
    logP_1 = log(half_life./(0.5.*cos(2*h*time_span0)));
    %plot(time_span0.^2, N^(2*alpha-1)*logP_1);
    %xlim([0 0.004])
    %drawnow
end
toc;

%% estimativa de tempo em horas
clear all; close all; clc;
tempo = 0;
% tempo_std = tempo para rodar N = 1000 particulas, tempo =
% linspace(0,1,100),alpha = 0.25

tempo_std = 13.72; 
intervalo = 1000:1000:15000;
for ii=intervalo
    tempo = tempo + (tempo_std*(ii/1000)^2)/3600;
end
disp(tempo)
%% eu estimo o tempo da simulacao com o seguinte codigo
clear all; clc; close all;
t = linspace(0,1,100);
N = 1000;
a = 0.25;
logP = 0*t;
tic;
for ii=1:N
    logP = 0;
    for i=1:N
        logP = logP + log(abs(cos(2*t/(min(i,N+1-i))^a)));
    end
    %plot(t.^2,N^(2*a-1)*logP)
    %drawnow
end
toc;


