clc
clear all
close all

return

%% EXPERIMENTAL DATA

% LOAD DATA
% load('CHIRP_INI_0')
% load('G:\Utilisateurs\jesse.dossantos\Nextcloud\Workspace\MATLAB\3D Printer Control\Modal Anlysis\data\05_18_21_Deterministic_Design_Updated_Test\CHIRP_INI_2')
% load('G:\Utilisateurs\jesse.dossantos\Nextcloud\Workspace\MATLAB\3D Printer Control\Modal Anlysis\data\05_12_21_Deterministic_Design_Test\CHIRP_INI_0')
% load('G:\Utilisateurs\jesse.dossantos\Nextcloud\Workspace\MATLAB\3D Printer Control\Modal Anlysis\data\05_21_21_S3\CHIRP_INI_0')

folder_data='G:\Utilisateurs\jesse.dossantos\Nextcloud\Workspace\MATLAB\3D Printer Control\Modal Anlysis\data\';
% Specimen 1
% load(strcat(folder_data,'05_12_21_Deterministic_Design_Test_S1\CHIRP_INI_1'))
% Specimen 2
% load(strcat(folder_data,'05_18_21_Deterministic_Design_Test_S2\CHIRP_INI_2'))
% Specimen 3
% load(strcat(folder_data,'05_26_21_Deterministic_Design_Test_S3\CHIRP_INI_3'))
% Specimen 4
load(strcat(folder_data,'05_28_21_Deterministic_Design_Test_S4\CHIRP_INI_5'))

files={'05_12_21_Deterministic_Design_Test_S1\CHIRP_INI_1',...
    '05_18_21_Deterministic_Design_Test_S2\CHIRP_INI_2',...
    '05_26_21_Deterministic_Design_Test_S3\CHIRP_INI_3',...
    '05_28_21_Deterministic_Design_Test_S4\CHIRP_INI_5',...
    '06_04_21_Deterministic_Design_Test_S5\CHIRP_INI_6'}
% color={'k','r','b','m','y'}

color=parula(16);

for k=1:16
    
    load(['G:/Utilisateurs/jesse.dossantos/Nextcloud/Workspace/MATLAB/LMM Identification/data/chirp1_b1_',num2str(k),'_p4_22_02_21']) %Layer 0
    
%     load(strcat(folder_data,files{k}))

    % Acquisiion frequency
    Fs=10240;

    % Create frequency vector
    N=size(sData{1}.Time,1);
    f=Fs*(0:(N/2))/N;

    for i=1:size(sData,2)

        data=sData{i}.Variables;
        time=seconds(sData{i}.Time);

        input=data(:,2);
        output=data(:,1)*125;   

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
    H_win_mean=abs(mean(transpose(G_win)));
    
    

    [val loc1]=max(f==160);
    [val loc2]=max(f==400);

    % Select points for model updating
    step=14;
    l_upt=loc1:step:loc2;
    f_upt=f(loc1:step:loc2);

    % FRF normalization by maximum amplitude
    H_norm_exp=H_win_mean(l_upt);
    
    H_inf=norm(H_norm_exp,Inf);

    % Find natural natural frequencies from experimental data
    [pks,locs]=findpeaks(H_norm_exp,'MinPeakHeight',350);
    fn_exp=f_upt(locs);

    P=pks(1);
    Q=pks(2);

    % PLOT
    figure(1)
    subplot(1,2,1)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    plot(k,H_inf,'o','MarkerEdgeColor',color(k,:),'MarkerFaceColor',color(k,:)), hold on
    ylabel('$ \left | G (\omega) \right |_\infty $','interpreter','latex')
    xlabel('Layer','interpreter','latex')
    set(gca,'FontSize',32,'TickLabelInterpreter','latex')
    

%     figure(2)
    subplot(1,2,2)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    % plot(f(loc1:loc2),H_norm_exp,'-','linewidth',2); hold on
    plot(f_upt,H_norm_exp,'-','linewidth',1,'Color',color(k,:)); hold on
%     plot(f,H_win_mean,'-','linewidth',1,'Color',color(k,:)); hold on
    % plot(P,Q,'o'); hold on
    % plot([0 1000],[0 1000],'k--'); hold on
    % plot([50 400],[1/sqrt(2) 1/sqrt(2)],'k-','linewidth',2); hold on
    % plot([50 400],[0.5301/sqrt(2) 0.5301/sqrt(2)],'k-','linewidth',2); hold on
    % plot([50 400],[0.74/sqrt(2) 0.74/sqrt(2)],'k-','linewidth',2); hold on
    % xlim([0,600])
    xlim([160,400])
    % ylim([1,1e4])
    set(gca,'FontSize',28,'TickLabelInterpreter','latex')
    ylabel('$ \left | G (\omega) \right | $','interpreter','latex')
    xlabel('Frequency [Hz]','interpreter','latex')
    % legend({'No Window','Hanning Window'},'interpreter','latex')
    % legend({'Experimental'},'interpreter','latex')
    
    figure(2)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    % plot(f(loc1:loc2),H_norm_exp,'-','linewidth',2); hold on
    % plot(f_upt,H_norm_exp,'-','linewidth',2); hold on
    plot(P,Q,'o','MarkerEdgeColor',color(k,:),'MarkerFaceColor',color(k,:)); hold on
    
    % plot([50 400],[1/sqrt(2) 1/sqrt(2)],'k-','linewidth',2); hold on
    % plot([50 400],[0.5301/sqrt(2) 0.5301/sqrt(2)],'k-','linewidth',2); hold on
    % plot([50 400],[0.74/sqrt(2) 0.74/sqrt(2)],'k-','linewidth',2); hold on
    xlim([0,1200])
    ylim([0,1200])
    % ylim([1,1e4])
    set(gca,'FontSize',28,'TickLabelInterpreter','latex')
    ylabel('$ Q $','interpreter','latex')
    xlabel('$ P $','interpreter','latex')
    % legend({'No Window','Hanning Window'},'interpreter','latex')
    % legend({'Experimental'},'interpreter','latex')

end
plot([0 2000],[0 2000],'k--'); hold on

legend({'L1','L2','L3','L4','L5','L6','L7','L8','L9','L10','L11','L12','L13','L14','L15','L16'},'interpreter','latex')

