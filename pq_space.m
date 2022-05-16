clc
clear all
close all

addpath functions/

return

%% EXPERIMENTAL DATA

% LOAD DATA
folder_data='data\'

% DECISION:
% 1 - Decrease natural frequency
% 2 - Increase natural frequency
% 0 - Stop

% PB1
files_cell{1}={'PB1\L1','PB1\L2','PB1\L3',...
    'PB1\L4','PB1\L5','PB1\L6',...
    'PB1\L7','PB1\L8'}
% decision={'2','2','2','2','2','2','2','0'}
decision_cell{1}={'k','k','k','k','k','k','k','k'}
line_cell{1}={':',':',':',':',':',':',':',':'}

% PB2
files_cell{2}={'PB2\L1','PB2\L2','PB2\L3','PB2\L4'}
% decision={'2','2','2','2','2','0'}
decision_cell{2}={'k','k','k','k','k','k'}
line_cell{2}={':',':',':',':',':',':'}

% PB3
files_cell{3}={'PB3\L1','PB3\L2','PB3\L3',...
    'PB3\L4','PB3\L5','PB3\L6',...
    'PB3\L7'}
% decision={'2','2','1','1','1','2','0'}
decision_cell{3}={'k','k','k','k','k','k','k'}
line_cell{3}={':',':','-','-','-',':'}

% PB4
files_cell{4}={'PB4\L1','PB4\L2','PB4\L3',...
    'PB4\L4','PB4\L5','PB4\L6',...
    'PB4\L7'}
% decision={'1','1','1','1','1','1','0'}
decision_cell{4}={'k','k','k','k','k','k','k'}
line_cell{4}={'-','-','-','-','-','-'}

% PB5
files_cell{5}={'PB5\L1','PB5\L2','PB5\L3',...
    'PB5\L4','PB5\L5','PB5\L6',...
    'PB5\L7','PB5\L8',...
    'PB5\L9','PB5\L10','PB5\L11',...
    'PB5\L12','PB5\L13','PB5\L14'}

% decision={'1','1','1','1','2','1','1','1','1','2','1','1','1','0'}
decision_cell{5}={'k','k','k','k','k','k','k','k','k','k','k','k','k','k'}
line_cell{5}={'-','-','-','-',':','-','-','-','-',':','-','-','-'}

% Colors
color=parula(4);
line_color=parula(6);
cor=1;

PQ={};
P=[];
Q=[];
H_inf=[];

TH=550;

font_size=24;

for g=1:5
    files=files_cell{g};
    decision=decision_cell{g};
    line_type=line_cell{g};
    
    clear P Q
    l=0;
    
    xi=[];
    yi=[];
    
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

        if length(pks)==1
            P=pks;
            Q=pks;

        else
            P(k)=pks(1);
            Q(k)=pks(2);

        end
        
        H_inf(k)=max(pks);


        figure(1)
        if k==1 | k==fix(length(files)/2) | k==length(files)
            l=l+1;
            plot(f_upt,H_norm_exp,'-','linewidth',2,'Color',color(l,:)); hold on
        end
        xlim([50,400])
        set(gca,'FontSize',font_size*0.8,'TickLabelInterpreter','latex')
        ylabel('$ \left | H (\omega) \right | [\frac{mm}{s^2V}]$','interpreter','latex')
        xlabel('Frequency [Hz]','interpreter','latex')
        legend({'Step$\;$0',strcat('Step$\;$',num2str(fix(length(files)/2))-1),strcat('Step$\;$',num2str(length(files)-1))},'interpreter','latex','location','best')
       saveas(1,strcat('figures/FRF_PB',num2str(g)),'epsc')
        
        figure(2)
        xi(k)=k-1;
        yi(k)=max(pks);

        if k==length(files)
            for n=1:length(files)-1
                p=plot(xi(n:n+1),yi(n:n+1),line_type{n},'linewidth',2,'color','k'); hold on
                uistack(p,'bottom')
                p.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
            plot([-2 -1],[0 1 ],'-','linewidth',2,'color','k');
            plot([-2 -1],[0 1 ],':','linewidth',2,'color','k');
            legend({'Decision I','Decision II'},'location','southwest','interpreter','latex')
            
            
            yi(1)
            yi(end)
            (yi(1)-yi(end))/yi(1)*100
        end


        p=plot(k-1,max(pks),'ko','MarkerSize',8,'MarkerFaceColor',decision{k}); hold on
        p.Annotation.LegendInformation.IconDisplayStyle = 'off';
        set(gca,'FontSize',font_size*0.8,'TickLabelInterpreter','latex')
        ylabel('$ \left | H (\omega) \right |_{\infty} [\frac{mm}{s^2V}]$','interpreter','latex')
        xlabel('Step','interpreter','latex')
        xlim([0 length(files)-1])
        ylim([1e4,3.2e4])
        saveas(2,strcat('figures/Hinf_PB',num2str(g)),'epsc')
        

        figure(3)
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        p=plot(P(k),Q(k),'ko','MarkerSize',8,'MarkerFaceColor',line_color(g,:)); hold on
        p.Annotation.LegendInformation.IconDisplayStyle = 'off';
        p=plot(P(:),Q(:),'k-','linewidth',2,'color',line_color(g,:)); hold on
        p.Annotation.LegendInformation.IconDisplayStyle = 'off';
        uistack(p,'bottom')

        set(gca,'FontSize',font_size,'TickLabelInterpreter','latex')
        ylabel('Q $[\frac{mm}{s^2V}]$','interpreter','latex')
        xlabel('P $[\frac{mm}{s^2V}]$','interpreter','latex')
        xlim([1.e4 3.5e4])
        ylim([1.e4 3.5e4])

    end
    close(1)
    close(2)
end


% Set threshold value
TH=700
x_axis = [0 7e4*25];
x_plot =[x_axis, fliplr(x_axis)];
y_plot=[x_axis+TH, fliplr(x_axis)-TH];

p=patch(x_plot, y_plot, 1,'facecolor', 'k', 'edgecolor', 'none', 'facealpha', 0.1);hold on
uistack(p,'bottom')
p=patch('Vertices',[0 0+TH;0 7e4*25; 7e4*25 7e4*25+TH],'facecolor', 'b', 'edgecolor', 'none', 'facealpha', 0.1);
uistack(p,'bottom')
p=patch('Vertices',[0 0-TH;7e4*25 0; 7e4*25 7e4*25-TH],'facecolor', 'r', 'edgecolor', 'none', 'facealpha', 0.1);
uistack(p,'bottom')

plot([0 7e4*25],[0 7e4*25],'k-'); hold on
plot([0 7e4*25],[0+TH 7e4*25+TH],'k--'); hold on


plot([-50 -60],[-100 -110],'-','linewidth',2,'color',line_color(1,:))
plot([-50 -60],[-100 -110],'-','linewidth',2,'color',line_color(2,:))
plot([-50 -60],[-100 -110],'-','linewidth',2,'color',line_color(3,:))
plot([-50 -60],[-100 -110],'-','linewidth',2,'color',line_color(4,:))
plot([-50 -60],[-100 -110],'-','linewidth',2,'color',line_color(5,:))

plot([0 7e4*25],[0-TH 7e4*25-TH],'k--'); hold on

legend({'Decision I','Decision II','Stop','Equal-peak','Threshold','PB1','PB2','PB3','PB4','PB5'},'interpreter','latex','FontSize',font_size,'location','northwest')

% Create arrow
annotation('arrow',[0.466745531019979 0.462145110410095],...
    [0.628390495867768 0.619576446280992],'HeadStyle','plain');

% Create arrow
annotation('arrow',[0.435462670872765 0.432308096740273],...
    [0.659123966942149 0.649793388429752],'HeadStyle','plain');

% Create arrow
annotation('arrow',[0.508806519453207 0.512092534174553],...
    [0.68753305785124 0.68595041322314],'HeadStyle','plain');

% Create arrow
annotation('arrow',[0.616719242902208 0.612644584647739],...
    [0.429010330578512 0.429493801652893],'HeadStyle','plain');

% Create arrow
annotation('arrow',[0.69111461619348 0.694137749737118],...
    [0.463068181818182 0.465650826446281],'HeadStyle','plain');