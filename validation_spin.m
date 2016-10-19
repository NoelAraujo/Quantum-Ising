clear all;close all;clc;

d = 1; % dimension
J = 1; 
h = 0; % external field
alpha = [0.6 0.7 0.8 0.9];
%alpha = [0.1 0.2 0.3 0.4];

N = 50000;

time_init = 1e-9;
time_end = 0.25;
time_steps = 500;
time_span = logspace(log10(time_init),log10(time_end),time_steps);

pack_colors = ['r','b','k','g'];

for aa=1:length(alpha)    
    y_t = 0;
    logP = 0;
    
    tic;
    for i=1:N
        logP = logP + log(abs(cos(2*time_span/(min(i,N+1-i))^alpha(aa))));
    end
    toc;
    spin = 1.*cos(2*h*time_span).*exp(logP);
    



    if alpha(aa) < d/2
            y_t = 0.*time_span -2^(5+2*alpha(aa)-d);
            y_t = y_t*pi^(d/2  -2 );
            y_t = y_t/(  (d-2*alpha(aa))*gamma(d/2) );
            y_t = y_t*(J^2).*time_span.^2;
            y_t = y_t*(N.^(1  -2*alpha(aa)/d ));
            y_t = exp(y_t);
            to_show(aa,:) = y_t;
    else    
            y_t = -2^(2 - 2*alpha(aa));
            y_t = y_t*sqrt(pi);
            y_t = y_t/(  (2*alpha(aa) - 1)*gamma(1/2) );
            y_t = y_t.* abs(4*J.*time_span./pi).^(1/alpha(aa));
            
%             y_t = -2^(1 + d - 2*alpha(aa));
%             y_t = y_t*pi^(d/2);
%             y_t = y_t/(  (2*alpha(aa)-d)*gamma(d/2) );
%             y_t = y_t.* abs(4*J.*time_span./pi).^(d/alpha(aa));
            y_t = exp(y_t);
            to_show(aa,:) = y_t;
    end
    
    %hh = plot(time_span,y_t,'--',time_span,spin,'o')
    %%hh = plot(time_span,y_t,'--',time_span,spin,strcat(pack_colors(aa),'o'));
    hh = plot(time_span,spin,strcat(pack_colors(aa),'o'));
    Legenda{aa} = strcat('\alpha = ',sprintf('%.2f',alpha(aa)));
    set(hh, 'MarkerFaceColor', get(hh(1), 'Color'));
    
    %plot(time_span,y_t,'--');
    hold on
    
end
legend(Legenda,'Location','best')

for aa=1:4
    plot(time_span,to_show(aa,:),strcat(pack_colors(aa),'--'))
    hold on
end
xlabel('Time')
ylabel('Spin')