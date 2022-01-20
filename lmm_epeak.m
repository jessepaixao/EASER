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

% Natural frequency
wn_s=sqrt(k1/m1)/(2*pi);

% Absorber system
m2=0.05*m1;
% k2=(wn_s*(0.05/(1+0.05)))^2*m2;
k2=(k1*m2)/m1;
c2=500;

% Natural frequency
w2_a=sqrt(k2/m2)/(2*pi);

omega=0:0.1:100*2*pi;
% omega=omega./(2*p;

%% OPTIMAL SOLUTION - DETERMINISTIC

% Numerical solution of the problem

x0 = [1e3];
k2_opt = fminsearch(@(x)recept2MCK_opt(omega,m1,k1,c1,m2,x(1),c2),x0)


%% Frequency response curve of inital design

freq=omega/(2*pi);
Hs=zeros(1,length(omega));
Hc_init=zeros(1,length(omega));

% Receptance of primary system
for i=1:length(omega)
    Hs(i)=1/(k1+1j*c1-omega(i)^2*m1);
    
end

% Receptance of coupled system (inital design)
Hc1=recept2MCK(omega,m1,k1,c1,m2,k2*1.2,c2);
Hc2=recept2MCK(omega,m1,k1,c1,m2,k2*0.8,c2);
Hc_opt=recept2MCK(omega,m1,k1,c1,m2,k2_opt,c2);

% DISPLAY RESULTS
figure1=figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
semilogy(freq,abs(Hs),'linewidth',2); hold on
semilogy(freq,abs(Hc1),'linewidth',2); hold on
semilogy(freq,abs(Hc2),'linewidth',2); hold on
semilogy(freq,abs(Hc_opt),'linewidth',2); hold on
set(gca,'FontSize',28,'TickLabelInterpreter','latex')
ylabel('$|G (\omega)|$ [m/N]','interpreter','latex')
xlabel('Frequency [Hz]','interpreter','latex')
ylim([3E-6,0.5E-2])

% Create textbox
annotation(figure1,'textbox',...
    [0.63564668769716 0.554130744087957 0.180280776376529 0.0827824790525397],...
    'String',{'System with TMD','(Optimally Tuned)'},...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create arrow
annotation(figure1,'arrow',[0.532071503680337 0.584647739221872],...
    [0.755198347107438 0.807851239669422],'LineWidth',2);

% Create arrow
annotation(figure1,'arrow',[0.439537329127234 0.387749737118822],...
    [0.430818181818181 0.471590909090909],'LineWidth',2);

% Create ellipse
annotation(figure1,'ellipse',...
    [0.517692954784437 0.823605371900827 0.00491482649842279 0.00929752066115674],...
    'FaceColor',[0 0 0]);

% Create ellipse
annotation(figure1,'ellipse',...
    [0.471688748685594 0.505940082644626 0.00491482649842279 0.00929752066115674],...
    'FaceColor',[0 0 0]);

% Create ellipse
annotation(figure1,'ellipse',...
    [0.453287066246057 0.436208677685949 0.00491482649842279 0.00929752066115674],...
    'FaceColor',[0 0 0]);

% Create ellipse
annotation(figure1,'ellipse',...
    [0.488118822292324 0.619576446280991 0.00491482649842279 0.00929752066115674],...
    'FaceColor',[0 0 0]);

% Create ellipse
annotation(figure1,'ellipse',...
    [0.542009463722397 0.56017561983471 0.00491482649842279 0.00929752066115674],...
    'FaceColor',[0 0 0]);

% Create ellipse
annotation(figure1,'ellipse',...
    [0.557125131440589 0.504648760330577 0.00491482649842279 0.00929752066115674],...
    'FaceColor',[0 0 0]);

% Create ellipse
annotation(figure1,'ellipse',...
    [0.593271293375394 0.427169421487602 0.00491482649842279 0.00929752066115674],...
    'FaceColor',[0 0 0]);

% Create arrow
annotation(figure1,'arrow',[0.572029442691903 0.62434279705573],...
    [0.494834710743802 0.558367768595042],'LineWidth',2);

% Create arrow
annotation(figure1,'arrow',[0.482386961093585 0.430599369085173],...
    [0.612119834710743 0.652892561983471],'LineWidth',2);

% Create textbox
annotation(figure1,'textbox',...
    [0.258391167192429 0.448242314335891 0.168396700962241 0.0827824790525397],...
    'String',{'System with TMD','(Undertuned)'},...
    'Interpreter','latex',...
    'FontSize',20,...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.256813880126182 0.644523306071428 0.168396700962241 0.0827824790525397],...
    'String',{'System with TMD','(Overtuned)'},...
    'Interpreter','latex',...
    'FontSize',20,...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.593585699263933 0.80995041233646 0.181179814457568 0.0501033066718046],...
    'String',{'System without TMD'},...
    'Interpreter','latex',...
    'FontSize',20,...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.50446898002103 0.835743801652893 0.0297055730809676 0.0392892561983476],...
    'String',{'PQ'},...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.481598317560462 0.639462809917355 0.0178759200841223 0.0377396694214882],...
    'String','P',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.457150368033649 0.515495867768595 0.0197160883280759 0.0372231404958684],...
    'String','P''',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.434279705573081 0.450929752066116 0.0197160883280759 0.0372231404958684],...
    'String','P''''',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.557045215562566 0.525826446280992 0.019716088328076 0.0372231404958684],...
    'String','Q''',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.588328075709779 0.443698347107438 0.019716088328076 0.0372231404958684],...
    'String','Q',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.539957939011567 0.576962809917355 0.019716088328076 0.0372231404958684],...
    'String','Q''''',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);


%% AMPLITUDE VS TUNING RATIO AND PQ SPACE

omega=0:0.1:1000*2*pi;

% Tuning ratio
gamma_vec=0.8:.01:2;

k2_vec=gamma_vec.^2*(k1/m1)*m2;


% gamma_vec=(vec.*2*pi).^2*m2;

% Initialize vectors
P=zeros(1,length(gamma_vec));
Q=zeros(1,length(gamma_vec));

for i=1:length(gamma_vec)
      
    Hc=recept2MCK(omega,m1,k1,c1,m2,k2_vec(i),c2);
    
    val=findpeaks(abs(Hc));
    
    P(i)=val(1);
    Q(i)=val(2);
     
end
    
% DISPLAY RESULTS
figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
semilogy(gamma_vec,P,'linewidth',2); hold on
semilogy(gamma_vec,Q,'linewidth',2); hold on
set(gca,'FontSize',32,'TickLabelInterpreter','latex')
xlabel('Tunning ratio $\gamma$','interpreter','latex')
legend({'$ \left | G (\omega_P)  \right | $','$ \left | G (\omega_Q)  \right | $'},'interpreter','latex')

figure()
set(gcf,'units','normalized','outerposition',[0 0 1 1])
plot(P,Q,'-o','linewidth',2,'MarkerFaceColor',lines(1)); hold on
plot([0 1E-3],[0 1E-3],'k--')
set(gca,'FontSize',32,'TickLabelInterpreter','latex')
xlabel('P','interpreter','latex')
ylabel('Q','interpreter','latex')
xlim([0 2E-4])
ylim([0 2E-4])

%%

% DISPLAY RESULTS
figure1=figure()
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% semilogy(freq,abs(Hs),'linewidth',2); hold on
semilogy(freq,abs(Hc1),'r','linewidth',2); hold on
% semilogy(freq,abs(Hc2),'b','linewidth',2); hold on
% semilogy(freq,abs(Hc_opt),'linewidth',2); hold on
set(gca,'FontSize',28,'TickLabelInterpreter','latex')
ylabel('$|H (\omega)|$','interpreter','latex')
xlabel('$\omega$','interpreter','latex')
ylim([3E-6,0.5E-3])

