clc
clear all
close all

addpath functions/

return


%% COMSOL INTIALIZATION
status = system(['C:\APPLICATIONS-FEMTO-ST\Comsol56\bin\win64\comsolmphserver.exe','&']);

pause(3)
addpath C:\APPLICATIONS-FEMTO-ST\Comsol56\mli

clear P
mphstart
import com.comsol.model.* 
import com.comsol.model.util.*

% Load Model
model = mphload('model\plate_beam.mph');


%% EXPERIMENTAL DATA

% LOAD DATA
% Sample PB_INI
load('data\PB_INI\L1')

% Acquisiion frequency
Fs=10240;

% Create frequency vector
N=size(sData{1}.Time,1);
f=Fs*(0:(N/2))/N;

clear G_win

for i=1:size(sData,2)

    data=sData{i}.Variables;
    time=seconds(sData{i}.Time);

    input=data(:,2);
    output=data(:,1);   
    
    window=hanning(2^14);
    overlap=length(window)*(3/4);
    L=2^nextpow2(N);
    nfft=L;

    % PSD of Input
    [pxx f]=cpsd(input,input,window,overlap,nfft,Fs);

    % Cross PSD of Input-Output
    [pxy f]=cpsd(input,output,window,overlap,nfft,Fs);

    % H1 Estimator
    G_win(:,i)=pxy./pxx;
end

% Accelerance
H_win_mean=abs(mean(transpose(G_win)).*(2*1j*pi.*f'));

[val loc1]=max(f==50);
[val loc2]=max(f==400);

% Select points for model updating
l_upt=loc1:12:loc2;
f_upt=f(loc1:12:loc2);

% FRF normalization by maximum amplitude
H_exp=H_win_mean(l_upt);

%% FEM MODEL PARAMTERS

% Define parameters
par.Hp = 150 ; % Heigth of the plate
par.Wp = 100 ;  % Width of the plate
par.Tp = 3 ;  % Thickness of the plate
par.Hb = 46.47; % Heigth of the beam S5 
par.Wb=10; % Width of the beam
par.Of=1; % Offset of the beam
par.Ey=2.174E9; % Young's modulus
par.Nu=0.35; % Poisson's Ratio
par.Rho=1120; % Density
par.k_f=1E13; % Foundation spring
par.damp1=0.0193; % Damping Ratio - Mode 1
par.damp2=0.02105; % Damping Ratio - Mode 2

% Setup parameters
model.param.set('Hp', [num2str(par.Hp,10),'[mm]']);
model.param.set('Wp', [num2str(par.Wp,10),'[mm]']);
model.param.set('Tp', [num2str(par.Tp,10),'[mm]']');
model.param.set('Wb', [num2str(par.Wb,10),'[mm]']);
model.param.set('Hb', [num2str(par.Hb,10),'[mm]']);
model.param.set('Of', [num2str(par.Of,10),'[mm]']);
model.param.set('Ey', [num2str(par.Ey,10),'[Pa]']);
model.param.set('Nu', [num2str(par.Nu,10)]);
model.param.set('Rho', [num2str(par.Rho,10),'[kg/m^3]']);
model.param.set('k_f', [num2str(par.k_f,10),'[N/(m*m)]']);


% FEM SOLUTION
freq=f_upt;

model.sol('sol4').feature('mo1').set('dampratio', [par.damp1 par.damp2 0 0 0 0 0 0 0 0]); %S1
model.study('std3').feature('frmod').set('plist', freq);
model.study('std3').run;

data=mphplot(model,'pg7','createplot','off');
   
% Accelerance from COMSOL
H=data{1, 1}{1, 1}.d;
cal_fac=(max(H_exp)/max(H));
H_sim=H.*cal_fac;

figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
plot(freq,H_sim,'b-','linewidth',2); hold on
plot(freq,H_exp,'r-','linewidth',2); hold on
set(gca,'FontSize',32,'TickLabelInterpreter','latex')
ylabel('$ \left | H (\omega) \right | $','interpreter','latex')
xlabel('Frequency [Hz]','interpreter','latex')
legend({'FEM','Experimental'},'interpreter','latex','location','best')


%% BEAM LENGTH OPTIMIZATION

lb = [40];
ub = [50];
options = optimset('Display','iter'); 
fun = @(x) fobjHinf(x,model);

x = fminbnd(fun,lb,ub,options);

%% OPTIMAL DESIGN

% Optimal beam length
x=46.07;
model.param.set('Hb', [num2str(x,10),'[mm]']);

% Set frequency range
freq=[50:400]';
model.study('std3').feature('frmod').set('plist', freq);

% Run FEM 
model.study('std3').run;

data=mphplot(model,'pg7','createplot','off');
    
H=data{1, 1}{1, 1}.d;
H_epeak_sim=H*cal_fac;

figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
plot(freq,H_epeak_sim,'-','linewidth',2); hold on
set(gca,'FontSize',32,'TickLabelInterpreter','latex')
ylabel('$ \left | H (\omega) \right | [\frac{mm}{s^2V}]$','interpreter','latex')
xlabel('Frequency [Hz]','interpreter','latex')
