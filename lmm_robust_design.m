clear all
close all
clc

% Add path of functions
addpath('functions')

return


%% SYSTEM PARAMTERS

% Primary system
m1=1;
k1=1E5;
c1=500;
wn_s=sqrt(k1/m1)/(2*pi);

% Absorber system
m2=0.05*m1;
k2=(k1*m2)/m1;
c2=500;

% Frequency range
omega=0:0.1:100*2*pi;

%% FRF INITIAL DESIGN

freq=omega/(2*pi);

Hs=zeros(1,length(omega));
Hc_init=zeros(1,length(omega));

% Receptance of primary system
for i=1:length(omega)
    Hs(i)=1/(k1+1j*c1-omega(i)^2*m1);
end

% Receptance of coupled system (inital design)
Hc_init=recept2MCK(omega,m1,k1,c1,m2,k2,c2);
    
% DISPLAY RESULTS
figure()
semilogy(freq,abs(Hs)); hold on
semilogy(freq,abs(Hc_init)); hold on
set(gca,'FontSize',12,'TickLabelInterpreter','latex')
ylabel('$ \left | H (\lambda) \right | $','interpreter','latex')
xlabel('$ \lambda  $','interpreter','latex')



%% DESIGN OPTIMIZATION - DETERMINISTIC

% Analytical solution of optimal values

mu=m2/m1;
gamma_opt_ana=1/(1+mu);
ka_opt_ana=(gamma_opt_ana*sqrt(k1/m1)*sqrt(m2))^2;

eta_opt_ana=sqrt(3*mu/(8*(1+mu)));
ca_opt_ana=eta_opt_ana*2*sqrt(m2*ka_opt_ana);

% Numerical solution of the problem

% x0 = [1e3,1];
% sol_det_num = fminsearch(@(x)recept2MCK_opt(omega,m1,k1,c1,m2,x(1),x(2)),x0)

x0 = [1e3];
sol_det_num = fminsearch(@(x)recept2MCK_opt(omega,m1,k1,c1,m2,x(1),c2),x0)


% recept_abs_damp_opt(omega,ms,ks,ma,x(1),x(2)),x0);

ka_opt_num=sol_det_num(1);
% ca_opt_num=sol_det_num(2);
ca_opt_num=c2;


% Frequency response curve of optimal design
Hc_opt=recept2MCK(omega,m1,k1,c1,m2,ka_opt_ana,ca_opt_ana);
Hc_opt_num=recept2MCK(omega,m1,k1,c1,m2,ka_opt_num,ca_opt_num);

% Propgation of uncertainty using Monte Carlo
Nm=1000;
k_dist=zeros(1,Nm);
h=zeros(Nm,length(omega));
h_norm=zeros(1,Nm);

for i=1:Nm
    k_dist(i)=unifrnd(0.6*k1,1.4*k1);
    h(i,:)=abs(recept2MCK(omega,m1,k_dist(i),c1,m2,ka_opt_num,ca_opt_num));
    h_norm(i)=norm(h(i,:),Inf);
end


% Display envelope of samples H
lamb_vec = [omega/wn_s, fliplr(omega/wn_s)];
env = [min(h,[],1) fliplr(max(h,[],1))];


figure()
semilogy(freq,abs(Hc_opt_num)); hold on
fill([freq, fliplr(freq)], env,'k','edgecolor','none','FaceAlpha',0.1); hold on
set(gca,'FontSize',12,'TickLabelInterpreter','latex')
ylabel('$|G (\omega)|$ [m/N]','interpreter','latex')
xlabel('Frequency [Hz]','interpreter','latex')
% legend({'LTVA - Analytical','Envelope of samples'},'Location','Best')

[~,h_legend] = legend({'LTVA Deterministic','Envelope of samples'},'Location','Best');
% PatchInLegend = findobj(h_legend, 'type', 'patch');
% set(PatchInLegend(1), 'FaceAlpha', 0.1);


%% ROBUST DESIGN OPTIMIZATION

% Risk level vector
e_vec=[0.01 0.1:0.1:1];
% e_vec=[0.01];

sol_rob_vec=zeros(length(e_vec),2);

% Create the subset of your random variable
for i=1:length(e_vec)
    i
    
    e=e_vec(i);          % Risk level
    beta=1e-10;     % Confidence level
    d = 3;
    Nm=fix((2/e)*(d-log(beta)))+1;
    ks_rnd=unifrnd(k1-0.2*(1-e)*k1,k1+0.2*(1-e)*k1,[1,Nm]);

    % Setup of optimization problem
    options = optimoptions('fmincon','Display','iter','PlotFcns',@optimplotfval,'Algorithm','sqp');
    problem.options = options;
    problem.solver = 'fmincon';
    problem.objective = @(q)objFcn_abs_damp(q);
    problem.x0 = [5e3 5];
    problem.lb = [0 0];
    problem.nonlcon = @(q)cnstFcn_abs_damp_mod(q,m1,ks_rnd,c1,m2,omega,c2);

    % Run optimization
    % Estimated time : 86s
    tic
    sol_rob = fmincon(problem)
    toc
    
    sol_rob_vec(i,:)=sol_rob;
    
end


%% Propgation of uncertainty using Monte Carlo
clear env_rob h_rob

Nm=10000;
% Solve the optimization above or run the saved results frol
% DESIGN_ROBUST.mat
sol= sol_rob_vec(1,:)
e=0.01;

k_dist=unifrnd(0.8*k1,1.2*k1,[1,Nm]);

for i=1:Nm   

    h_rob(i,:)=abs(recept2MCK(omega,m1,k_dist(i),c1,m2,sol(1),c2));
end

Hc_rob=abs(recept2MCK(omega,m1,k1,c1,m2,sol(1),c2));

env_rob = [min(h_rob,[],1), fliplr(max(h_rob,[],1))];


co=lines(7)

figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])

fill([freq, fliplr(freq)],env,co(1,:),'edgecolor','none','FaceAlpha',0.1); hold on
plot(freq,abs(Hc_opt_num),'Color',co(1,:),'linewidth',2); hold on

fill([freq, fliplr(freq)], env_rob,co(2,:),'edgecolor','none','FaceAlpha',0.1); hold on
plot(freq,abs(Hc_rob),'Color',co(2,:),'linewidth',2); hold on

set(gca,'FontSize',28)

[~,h_legend] = legend({'Equal-peak design','Mean sample - deterministic','Robust equal-peak design','Mean sample - robust'},'Location','Best','interpreter','latex');
PatchInLegend = findobj(h_legend, 'type', 'patch');

set(PatchInLegend(1), 'FaceAlpha',0.1);
set(PatchInLegend(2), 'FaceAlpha',0.1);
% set(gca,'YScale','log')
ylabel('$|G (\omega)|$ [m/N]','interpreter','latex')
xlabel('Frequency [Hz]','interpreter','latex')
set(gca,'FontSize',28,'TickLabelInterpreter','latex')

set(gcf,'PaperOrientation','landscape');;
set(gcf,'PaperSize',[38,21])
print(gcf, '-dpdf','-fillpage', 'figures/H_robust_v2.pdf');

% TRADE-OFF PERFORMANCE x ROBUSTNESS

co=lines(7);

figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
plot(e_vec*100,sol_rob_vec(:,2)*1e-4,'o-','linewidth',2,'color',co(1,:),'markerfacecolor',co(1,:))
set(gca,'FontSize',28,'TickLabelInterpreter','latex')
ylabel('$g^*$','interpreter','latex')
xlabel('$ \varepsilon$ $[\%]$','interpreter','latex')
xlim([0 100])

set(gcf,'PaperOrientation','landscape');;
set(gcf,'PaperSize',[38,21])
print(gcf, '-dpdf','-fillpage', 'figures/tradeoff_v2.pdf');

