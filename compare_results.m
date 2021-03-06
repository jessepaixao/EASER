clc
clear all
close all

addpath functions/

return

%% EXPERIMENTAL DATA

% Load data
folder_data='data\';
files={'P1\L1','P2\L1','P3\L1',...
    'PB1\L1','PB2\L1','PB3\L1',...
    'PB4\L1','PB5\L1'};

%Setup colors
color=parula(length(files));
PQ={};

close all
font_size=32;

for k=1:length(files)
    
    load(strcat(folder_data,files{k}))

    % Acquisiion frequency
    Fs=10240;

    % Create frequency vector
    N=size(sData{1}.Time,1);
    f=Fs*(0:(N/2))/N;
    
    clear G_win
    
    for i=1:size(sData,2)    
        
        data=sData{i}.Variables;
        time=seconds(sData{i}.Time);

        input=data(:,2); % [mm/s2]
        output=data(:,1);   % [V]

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

    % Receptance
    H_win_mean=abs(mean(transpose(G_win))).*(2*pi.*f');

    [val loc1]=max(f==50);
    [val loc2]=max(f==400);

    % Select points for model updating
    step=14;
    l_upt=loc1:step:loc2;
    f_upt=f(loc1:step:loc2);

    % FRF normalization by maximum amplitude
    H_norm_exp=H_win_mean(l_upt);
    
    % Find natural natural frequencies from experimental data
    [pks,locs]=findpeaks(H_norm_exp,'MinPeakHeight',1e4);
    fn_exp=f_upt(locs);
    PQ{k}=pks;

    figure(1)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    plot(f_upt,H_norm_exp,'-','linewidth',2,'Color',color(k,:)); hold on
    xlim([50,400])
    ylim([0 7e4])
    set(gca,'FontSize',font_size,'TickLabelInterpreter','latex')
    ylabel('$ \left | H (\omega) \right | [\frac{mm}{s^2V}]$','interpreter','latex')
    xlabel('Frequency [Hz]','interpreter','latex') 
    
    if k==length(files)
        legend({'P1','P2','P3','PB1','PB2','PB3','PB4','PB5'},'interpreter','latex','location','best')
    end
       
    figure(2)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    if length(pks)==1
        plot(PQ{k},PQ{k},'ko','MarkerSize',12,'MarkerFaceColor',color(k,:)); hold on
    else
        plot(PQ{k}(1),PQ{k}(2),'ko','MarkerSize',12,'MarkerFaceColor',color(k,:)); hold on
    end
    set(gca,'FontSize',font_size,'TickLabelInterpreter','latex')
    ylabel('Q $[\frac{mm}{s^2V}]$','interpreter','latex')
    xlabel('P $[\frac{mm}{s^2V}]$','interpreter','latex')
    xlim([0 7e4])
    ylim([0 7e4])
    
end
plot([0 7e4],[0 7e4],'k--')
legend({'P1','P2','P3','PB1','PB2','PB3','PB4','PB5'},'interpreter','latex','location','best')

% FIGURES 7A and 7B