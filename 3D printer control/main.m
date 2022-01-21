%% CONTROL 3D PRINTER : SELF-DESIGN MANUFACTURING PARADIGM
% 
% Functions of slicer implementation in MATLAB
% Copyright (c) 2018, Sunil Bhandari
% 
% Function stlwrite
% Copyright (c) 2018, Sven Holcombe
%
% Created by
% Copyright (c) 2020 Jesse Paixao


close all
clear all
clc

addpath('functions')
addpath('functions/Slicer')

%% PREPARE ACQUISITION

% SET SAMPLE NAME
experiment_name='PB1';

% NI ACQUISITION SETUP
NI = daq("ni");
NI.Rate = 10240; % Acquisition frequency
ch0=addinput(NI,"cDAQ1Mod1","ai0","Voltage");
ch0.Name='Vibrometer - V100';
ch1=addinput(NI,"cDAQ1Mod1","ai1","Voltage");
ch1.Name='Amplifier Output';
ch2=addinput(NI,"cDAQ1Mod1","ai2","Voltage");
ch2.Name='Amplifier Input';
scale_vibrometer=25; % Scale signal from vibrometer
scale_amplifier=1; % Scale signal from amplifier

% SETUP EXCITATION SIGNAL
fs=NI.Rate;  % sampling frequency
duration=10; % duration
t=0:1/fs:duration; % input vector time 
exc=0.2*chirp(t,10,duration,1500,'logarithmic'); % generate excitation signal
N_rept=10; % number of repetition acquisition

sData={};

% ACQUISITION
for i = 1:N_rept
    % Send excitation signal
    sound(exc,fs);
    % Start acquistion of input and output
    sData{i} = read(NI,seconds(duration));
    % Pause between acqusition
    pause(1)
end

% Save data startup
save(strcat('data/',experiment_name,'/',experiment_name,'_INI.mat'),'sData');

%% COMPUTE AND PLOT FRF INITIAL DESIGN

fs=10240; % Sampling frequency
load(strcat('data/',experiment_name,'/',experiment_name,'_INI.mat'),'sData');

N=size(sData{1}.Time,1);
freq=fs*(0:(N/2))/N;

clear G_exp

for i=1:size(sData,2)

    data=sData{i}.Variables;
    time=seconds(sData{i}.Time);

    x=data(:,2);
    y=data(:,1)*scale_vibrometer;   
    
    window=hanning(2^14);
    overlap=length(window)*(3/4);
    L=2^nextpow2(N);
    nfft=L;

    % PSD of Input
    [pxx f]=cpsd(x,x,window,overlap,N,fs);

    % Cross PSD of Input-Output
    [pxy f]=cpsd(x,y,window,overlap,N,fs);
    
    % H1 Estimator
    G_exp(:,i)=pxy./pxx;  
    
end

[val loc1]=max(f==50);
[val loc2]=max(f==400);

H_exp=(2*pi.*f').*mean(transpose(G_exp));

H_exp=H_exp(loc1:loc2);
f=f(loc1:loc2);

[pks pn]=findpeaks(abs(H_exp),'MinPeakProminence',100); 

fn_exp=f(pn);

H_inf_exp(1)=max(pks);

% PLOT : FRF INITIAL DESIGN
figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
plot(f,abs(H_exp),'-','linewidth',2); hold on
xlim([100,400])
set(gca,'FontSize',32,'TickLabelInterpreter','latex')
ylabel('$ \left | H (\omega) \right | $','interpreter','latex')
xlabel('Frequency [Hz]','interpreter','latex')


%% INITIALIZE 3D PRINTER

format short;
warning('off');
global bed_hight bed_width bed_temp nozzle_dim filament_dim filament_temp bed_center00 bed_x0 bed_y0 bed_z0 init_z
global layer_hight shell_thick top_bottom_thick infill_type top_bottom_type infill_density print_speed skirt_dis

% DEFINE 3D PRINTER PARAMETERS
% 3D  printer parameters
bed_hight=235;      % bed hight
bed_width=250;      % bed width(mm)% fprintf(fid,';end\n');
bed_temp=0;        % bed temperature(C)
nozzle_dim=0.4;     % nozzle diameter(mm)
filament_dim=1.75;  % filament diameter(mm)
filament_temp=250;  % filament temperature(C)
% bed_center00=false; % bed center coordinate at center of the bed X0.0000 Y0.0000 
bed_center00=true; % bed center coordinate
% SET REFERENCE POINT 3D PRINTER
bed_z0=66.88;
bed_y0=163;
bed_x0=156.5;
init_z=110;
% Slicer parameters
layer_hight=0.1; %layer height(mm)
shell_thick=1; %shell thickness(mm)
top_bottom_thick=0.5; %top and bottom surface thickness(mm)
infill_type='rec'; %rec and offset is optional
top_bottom_type='rec'; %rec and offset is optional
infill_density=100; %0->100 by percent
print_speed=30; %nozzle move speed(mm/s)
skirt_dis=0; %skirt distance

% Functionalities
adptive=0; %adptive:0 equal thickness, 1 adptive
% File input/output 
filename='files/cube.stl';
gfilename='files/mass.gcode';

% CREATE GEOMETRY TO PRINT
a=10; % Layer size
h=0.2; % Layer height

% DESIGN MODIFICATION - DECISION I
% Coordinates cube
X=[0 0 0;
   a 0 0;
   a a 0;
   0 a 0;
   0 0 0;
   0 0 h;
   a 0 h;
   a a h;
   0 a h;
   0 0 h;
   0 0 0];

dt=delaunayTriangulation(X);
tetramesh(dt, 'FaceColor', 'cyan');
[F,P] = freeBoundary(dt);
stlwrite('files/cube_M.stl',F,P);

a=10; % Layer size
h=0.2; % Layer height

% DESIGN MODIFICATION - DECISION II
% Coordinates cube
X=[0 0 0;
   2*a 0 0;
   2*a a 0;
   0 a 0;
   0 0 0;
   0 0 h;
   2*a 0 h;
   2*a a h;
   0 a h;
   0 0 h;
   0 0 0];

dt=delaunayTriangulation(X);
tetramesh(dt, 'FaceColor', 'cyan');
[F,P] = freeBoundary(dt);
stlwrite('files/cube_R.stl',F,P);


% START CLOSED LOOP
% NI ACQUISITION SETUP
NI = daq("ni");
NI.Rate = 10240; % Acquisition frequency
ch0=addinput(NI,"cDAQ1Mod1","ai0","Voltage");
ch0.Name='Vibrometer - V100';
ch1=addinput(NI,"cDAQ1Mod1","ai1","Voltage");
ch1.Name='Amplifier Output';
ch2=addinput(NI,"cDAQ1Mod1","ai2","Voltage");
ch2.Name='Amplifier Input';
scale_vibrometer=1; % Scale signal from vibrometer
scale_amplifier=1; % Scale signal from amplifier

% SETUP EXCITATION SIGNAL
fs=NI.Rate;  % sampling frequency
duration=10; % duration
t=0:1/fs:duration; % input vector time 
exc=0.2*chirp(t,10,duration,1500,'logarithmic'); % generate excitation signal
N_rept=10; % number of repetition acquisition

sData={};

font_size=44;

P=[];
Q=[];

str = 'Y';
mod = 'M';

z0_M=bed_z0;
z0_R=bed_z0;
X0=bed_x0;

k=0;

color=parula(30);

close all

while str == 'Y'
    
    clear H_exp H_inf
    
    pause(50)
    
    sData={};
    
    % MODAL ANALYSIS MEASUREMENTS
    
    for i = 1:N_rept
        % Send excitation signal
        sound(exc,fs);
        % Start acquistion of input and output
        sData{i} = read(NI,seconds(duration));
        % Pause between acqusition
        pause(1)
    end

    % Save data startup
    save(strcat('data/',experiment_name,'/',experiment_name,'_L',num2str(k),'.mat'),'sData');
    
    N=size(sData{1}.Time,1);
    
    freq=fs*(0:(N/2))/N;
    
    clear G_exp
    
    for i=1:size(sData,2)

        data=sData{i}.Variables;
        time=seconds(sData{i}.Time);

        x=data(:,2);
        y=data(:,1);   

        window=hanning(2^14);
        overlap=length(window)*(3/4);
        L=2^nextpow2(N);
        nfft=L;

        % PSD of Input
        [pxx f]=cpsd(x,x,window,overlap,N,fs);

        % Cross PSD of Input-Output
        [pxy f]=cpsd(x,y,window,overlap,N,fs);

        % H1 Estimator
        G_exp(:,i)=pxy./pxx;  

    end

    [val loc1]=max(f==50);
    [val loc2]=max(f==400);
    
    clear loc val
  
    
    H_exp=(2*pi.*f').*mean(transpose(G_exp));

    H_exp=H_exp(loc1:loc2);
    f=f(loc1:loc2);

    [pks pn]=findpeaks(abs(H_exp),'MinPeakProminence',100); 
    
    P(k)=pks(1);
    Q(k)=pks(2);
      
    fn_exp=f(pn(2));

    H_inf_exp(1)=max(pks);

    % PLOT 1 : FRF
    figure(1)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    plot(f,abs(H_exp),'-','Color',color(k,:),'linewidth',2); hold on
    xlim([100,400])
    set(gca,'FontSize',32,'TickLabelInterpreter','latex')
    ylabel('$ \left | H (\omega) \right | $','interpreter','latex')
    xlabel('Frequency [Hz]','interpreter','latex')
    
    
    % PLOT 2 : PQ SPACE
    figure(2)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    plot(P(k),Q(k),'ko','MarkerSize',8,'MarkerFaceColor',color(k,:)); hold on
    plot([0 1500],[0 1500],'k--')
    xlim([0,1500])
    ylim([0,1500])
    set(gca,'FontSize',32,'TickLabelInterpreter','latex')
    ylabel('P','interpreter','latex')
    xlabel('Q','interpreter','latex')
    
    

    prompt = 'Do you want print a new layer more? Y/N [Y]: ';
    str = input(prompt,'s');
    if str=='N'
        break
    end
    
    prompt = 'Do you want to increse mass or rigidity? M/R [M]: ';
    mod = input(prompt,'s');   
    if isempty(mod)
        mod = 'M';
    end
    
    if mod=='M'
        filename='files/cube_M.stl';
        bed_x0=X0-36.05;
        bed_z0=z0_M+h;
        z0_M=z0_M+h;
    else
        filename='files/cube_R.stl';
        bed_x0=X0;
        bed_z0=z0_R+h;
        z0_R=z0_R+h;
    end
    
    slicer()

    % 3D PRINTER COMMUNICATION

    filetext = readlines("files/mass.gcode")
    block_size=10;

    printer = tcpip('10.42.0.175', 23);
    printer.OutputBufferSize = 100000;
    printer.InputBufferSize = 100000;


    vec=1:10:length(filetext);

    for i=vec
       fopen(printer)

       for m=0:block_size-1 
           if i+m<length(filetext)
               fprintf(printer,filetext(i+m));
               filetext(i+m)
           end
       end

       if i==vec(end)
           break;
       end

       txt=fscanf(printer);
       if isempty(txt)==1
              txt='none';
       end

       while strcmp(txt(1:2),'ok')==0
          txt_aux=fscanf(printer);
          if isempty(txt_aux)==0
              txt=txt_aux;
          end
       end
       fclose(printer)
    end
       
    k=k+1;

end
