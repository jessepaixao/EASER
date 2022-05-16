close all
clear all
clc

addpath functions/

return


%% LOAD EXPERIMENTAL DATA - STARTUP

% LOAD DATA
load('data\PB_INI\L1')

% Acquisiion frequency
Fs=10240;

% Create frequency vector
N=size(sData{1}.Time,1);
f=Fs*(0:(N/2))/N;

% CUMPUTE FRF
for i=1:size(sData,2)

    % Store variables
    data=sData{i}.Variables;
    time=seconds(sData{i}.Time);

    input=data(:,2);
    output=data(:,1);

    % Compute FFT
    X=fft(input);
    Y=fft(output); 
    
    % Compute spectral densities
    sxx=conj(X).*X;
    sxy=conj(X).*Y;

    % Estimate FRF (H1 Estimator)
    G_win(:,i)=sxy./sxx;
end

% Select interval for LMM identification
[val loc1]=max(f==145);
[val loc2]=max(f==400);

global G

% Mean of FRF (modulus)
G_exp=(mean(transpose(G_win(1:N/2+1,:))));

G=G_exp(loc1:loc2);
f=f(loc1:loc2);
G = G.';
f = f.';

G = G./(1i*2*pi*f);

% FRF fitting

f0 = [(2*f(1)+f(end))/2 ; (2*f(1)+2*f(end))/3];
maxiter = 10;
[q,err] = Fitting(G,f,f0,maxiter);

u = q(1);
v = q(2);
r = q(3:4);
% r = 1i*abs(q(3:4));
s = q(5:6);
Hs = zeros(size(G));
W = 2*pi*f;
for k=1:length(r)
    Hs = Hs + r(k)./(1i*W-s(k));
end
Hs = Hs+u+1i*(W-mean(W))*v;

font_size=24;

figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
clf
plot(real(G),imag(G),'b','linewidth',2)
hold on
axis('equal')
plot(real(Hs),imag(Hs),'or','linewidth',2)
set(gca,'FontSize',font_size,'TickLabelInterpreter','latex')
ylabel('$ Imag[ H (\omega)] $','interpreter','latex')
xlabel('$ Real[ H (\omega)] $','interpreter','latex')
legend({'Experimental','Identified'},'interpreter','latex','location','best')

% 
figure(2)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
clf
semilogy(W/2/pi,abs(G).*W.^2,'b','linewidth',2)
hold on
semilogy(W/2/pi,abs(Hs).*W.^2,'r','linewidth',2)
set(gca,'FontSize',font_size,'TickLabelInterpreter','latex')
ylabel('$ \left | H (\omega) \right | [\frac{mm}{s^2V}]$','interpreter','latex')
xlabel('Frequency [Hz]','interpreter','latex') 
xlim([145,400])
legend({'Experimental','Identified'},'interpreter','latex','location','best')


wn=abs(s)/2/pi;
damp=-real(s)./abs(s)*100;
