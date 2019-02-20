close all 
clear all 
clc
%% READ FILE PARAMETERS AND DATA MATRIX
var = 37; %number of variables being looked at
% Open the file.
fullFileName = fullfile(pwd,'Sine','2um','1mm IED','10000_Hz_-102.0_mA.txt');
fileID = fopen(fullFileName, 'rt');
% Read the first line of the file.
textLine = strtrim(fgetl(fileID));
while ischar(textLine)
	% Print this line to the command window.
	fprintf('%s\n', textLine);
	% Now parse the line
	if contains(textLine(1:2), 'dt')
		dt = str2double(textLine(4:end))
	elseif contains(textLine(1:6), 'tfinal', 'IgnoreCase', true)
		tfinal = str2double(textLine(8:end))
	elseif contains(textLine(1:6), 'fsin.f', 'IgnoreCase', true)
		f = str2double(textLine(8:end))
	elseif contains(textLine(1:8), 'fsin.amp', 'IgnoreCase', true)
		f_amp = str2double(textLine(10:end))
	elseif contains(textLine(1:4), 'nseg', 'IgnoreCase', true)
		nseg = str2double(textLine(6:end))
	elseif contains(textLine(1:2), 'L ', 'IgnoreCase', true)
		L = str2double(textLine(3:end))
	elseif contains(textLine(1:4), 'diam', 'IgnoreCase', true)
		diam = str2double(textLine(6:end))
	elseif contains(textLine(1:8), 'stim_amp ', 'IgnoreCase', true)
		stim_amp = str2double(textLine(10:end))
	elseif contains(textLine, 'batch_run from', 'IgnoreCase', true)
		% Now, after this batch_run line we're expecting 5 lines with 5 numbers per line.
		m = zeros(tfinal/dt+1, nseg*var); % Preallocate space for speed.
		for row = 1 : tfinal/dt+1
			% Read the next line.
			textLine = fgetl(fileID);
			m(row, :) = sscanf(textLine, '%f');
		end
	end
	% Read the next line.
	textLine = fgetl(fileID);
end
% All done reading all lines, so close the file.
fclose(fileID);
%% SPLIT DATA INTO SEPARATE AXON SEGMENTS 
axon = mat2cell(m, [tfinal/dt+1], repmat(var, 1, nseg));
t = [0:dt:tfinal]';

%% CALC CONDUCTION VELOCITY
%first read spike time file
ffname = fullfile(pwd,'Sine','2um','1mm IED','10000_Hz_0.0_mA_ST.txt');
fid = fopen(ffname);
A = [];
tline = fgetl(fid);
ii = 1;
jj = 1;
while ischar(tline)
  while ~isempty(tline)
      A(ii,jj) = str2num(tline);
      ii = ii + 1;
      disp(tline)
      tline = fgetl(fid);
  end
  ii = 1;
  jj = jj + 1;
  tline = fgetl(fid);
end
fclose(fid);

%% Find spikes occuring after test stim
B = A.*(A<30);

C = zeros(1,size(B,2));
D = zeros(1,size(B,2));
for i=1:size(B,2)-1
    C(1,i) = find(B(:,i)~=0, 1, 'first');
end 
for i=1:size(B,2)-1
D(1,i) = B(C(1,i),i);
end 
% Remove first 3 segments because test stim occurs in the 4th
E = D(:,4:size(B,2))
%% calculate conduction velocity 
deltax = L/nseg;
deltat = zeros(1,length(E));
cond_vel = zeros(1,length(E));

for i = 1:length(E)-1
    deltat(1,i) = E(1,i+1)-E(1,i);
end 

for i = 1:length(E)-1
    cond_vel(1,i) = deltax/deltat(1,i)/1000
end 
    
condvel_final = mean(cond_vel(1:31));
%% calculate conduction velocity (only one spike)
deltax = L/nseg;
deltat = zeros(1,length(A));
cond_vel = zeros(1,length(A));

for i = 1:length(A)-1
    deltat(1,i) = A(1,i+1)-A(1,i);
end 

for i = 1:length(A)-1
    cond_vel(1,i) = deltax/deltat(1,i)/1000
end 
    
condvel_final = mean(cond_vel(5:31));
   
  %% Reorganise for plotting 
    start_onset = find(t==0);
    end_onset = find(t==50);
    start_test = find(t==95);
    end_test = find(t==145);
    rows_r = (end_test-start_test)+(end_onset-start_onset);
  
    for i = 1:18
    axondata(:,i) = axon{1,2*i}(:,1);
    end 
    
    specaxondata = zeros(rows_r,18);
    specaxondata = axondata([start_onset:end_onset,start_test:end_test],:);
    %specaxondata = axondata([1000:2000],[6000:8000],:);
    
    %% 3D plot of Vm over axon length 
    p = figure
    hold on 
    l = L/1000;
    y = 0.5:0.5:l;
    x = 1:rows_r+2;
    z = specaxondata(:,[1:18]);
    waterfall(x', y', z'); 
    axis([0 rows_r 0 9 -120 50]);
     set(gca,...
     'XTickLabel',{'0','20','40','60','80','100'},...
     'XTick',[0 4000 8000 12000 16000 20000],...
    'YTick',[0:1:9],...
    'ZTick',[],...
    'Zcolor','none');
%     xlabel('Time, ms')
%     ylabel('Length along axon, mm')
    set(gcf, 'color','white')
    view(-35,66)
    box off
    hold off 
    
    %saveas(p,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/vm_3D.svg'],'svg');
    
   %% Params over whole time 
%    1-vm
%    2-e_ext
%    3-ik
%    4-ina
%    5-ki
%    6-nai
%    7-ko
%    8-nao
%    9-ina_hcn
%    10-ik_hcn
%    11-ina_nav1p7
%    12-ina_nav1p8
%    13-ina_nav1p9
%    14-ina_NaKpump
%    15-ik_NaKpump
%    16-ik_knatype
%    17-ik_km
%    18-ik_kdr
%    19-ik_ka
%    20-nf_hcn
%    21-ns_hcn
%    22-n_ka
%    23-h_ka
%    24-n_kdr
%    25-ns_km
%    26-nf_km
%    27-w_knatype
%    28-m_nav1p8
%    29-h_nav1p8
%    30-s_nav1p8
%    31-u_nav1p8
%    32-m_nav1p9
%    33-h_nav1p9
%    34-s_nav1p9
%    35-m_nav1p7
%    36-h_nav1p7
%    37-s_nav1p7
   
   
   c = jet(5);
    h = figure 
    hold on 
%     plot(t(:,1),axon{1,1}(:,1),'color',c(1,:));
%     plot(t(:,1),axon{1,5}(:,1),'color',c(2,:));
    plot(t(:,1),axon{1,14}(:,1),'color',c(1,:));
    plot(t(:,1),axon{1,16}(:,1),'color',c(2,:));
    plot(t(:,1),axon{1,18}(:,1),'color',c(3,:));
    plot(t(:,1),axon{1,20}(:,1),'color',c(4,:));
    plot(t(:,1),axon{1,22}(:,1),'color',c(5,:));
%     plot(t(:,1),axon{1,32}(:,1),'color',c(8,:));
%     plot(t(:,1),axon{1,36}(:,1),'color',c(9,:));
    axis([0 tfinal -120 40])
%     legend({'node 1','node 5','node 14','node 16','node 18','node 20','node 22','node 32','node 36'}, 'Location','northeast');
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northeast');
    xlabel('Time, ms');
    ylabel('Membrane potential, mV');
    set(gcf, 'color','white')
    hold off 
    %%
    %saveas(h,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/vm.svg'],'svg');
    
    for i = 2:1:37
    k = figure 
    hold on 
%     plot(t(:,1),axon{1,1}(:,i),'color',c(1,:));
%     plot(t(:,1),axon{1,5}(:,i),'color',c(2,:));
    plot(t(:,1),axon{1,14}(:,i),'color',c(1,:),'linewidth',0.8);
    plot(t(:,1),axon{1,16}(:,i),'color',c(2,:),'linewidth',0.8);
    plot(t(:,1),axon{1,18}(:,i),'color',c(3,:),'linewidth',0.8);
    plot(t(:,1),axon{1,20}(:,i),'color',c(4,:),'linewidth',0.8);
    plot(t(:,1),axon{1,22}(:,i),'color',c(5,:),'linewidth',0.8);
%     plot(t(:,1),axon{1,32}(:,i),'color',c(8,:));
%     plot(t(:,1),axon{1,36}(:,i),'color',c(9,:));
%     legend({'node 1','node 5','node 14','node 16','node 18','node 20','node 22','node 32','node 36'}, 'Location','northwest');
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northeast');
    xlim([0 tfinal]);
    xlabel('Time, ms');
    set(gcf, 'color','white')
    if i == 2 
        ylabel('e_extracellular (mV)');
        %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/e_extracellular.svg'],'svg');
%         ylim([0 1]);
    elseif i == 3
        ylabel('ik (mA/cm^2)');
        legend('Location','northeast');
        %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ik.svg'],'svg');
%         ylim([0 1]);
    elseif i == 4
        ylabel('ina (mA/cm^2) ');
        legend('Location','southeast');
        %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ina.svg'],'svg');
%         ylim([0 1]);
    elseif i == 5
         ylabel('K^+ intracellular concentration (mM)');
         %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ki.svg'],'svg');
         legend('Location','southwest');
         legend boxoff
%         ylim([0 1]);
    elseif i == 6
        ylabel('Na^+ intracellular concentration (mM)');
          % saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/nai.svg'],'svg');
         legend boxoff
%         ylim([0 1]);
    elseif i == 7
        ylabel('K^+ periaxonal space concentration (mM)');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ko.svg'],'svg');
%         ylim([0 1]);
            legend('Location','northeast');
            legend boxoff
    elseif i == 8
        ylabel('Na^+ periaxonal space concentration (mM)');
        legend('Location','southeast');
         legend boxoff   
        %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/nao.svg'],'svg');
%         ylim([0 1]);
    elseif i == 9
        ylabel('ina_hcn (mA/cm^2)');
        legend('Location','northeast');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ih_hcn.svg'],'svg');
%         ylim([0 1]);
    elseif i == 10
        ylabel('ik_hcn (mA/cm^2)');
        legend('Location','northeast');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ih_hcn.svg'],'svg');
%         ylim([0 1]);
    elseif i == 11
        ylabel('ina_nav17 (mA/cm^2)');
        legend('Location','southeast');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ina_nav17.svg'],'svg');
%         ylim([0 1]);
    elseif i == 12
        ylabel('ina_nav18 (mA/cm^2)');
        legend('Location','southwest');
          % saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ina_nav18.svg'],'svg');
%         ylim([0 1]);
    elseif i == 13
        ylabel('ina_nav19 (mA/cm^2)');
        legend('Location','northeast');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ina_nav19.svg'],'svg');
%         ylim([0 1]);
    elseif i == 14
        ylabel('ina_NaKpump (mA/cm^2)');
        legend('Location','northeast');
          % saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ink_NaKpump.svg'],'svg');
%         ylim([0 1]);
    elseif i == 15
        ylabel('ik_NaKpump (mA/cm^2)');
        legend('Location','northeast');
          % saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ink_NaKpump.svg'],'svg');
%         ylim([0 1]);
    elseif i == 16
        ylabel('ik_knatype (mA/cm^2)');
        legend('Location','northeast');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ik_knatype.svg'],'svg');
%         ylim([0 1]);
    elseif i == 17
        ylabel('ik_km (mA/cm^2)');
        legend('Location','northeast');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ik_kmtype.svg'],'svg');
%         ylim([0 1]);
    elseif i == 18
        ylabel('ik_Kdr (mA/cm^2)');
        legend('Location','northeast');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ik_Kdr.svg'],'svg');
%         ylim([0 1]);
    elseif i == 19
        ylabel('ik_ka (mA/cm^2)');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/jika_ka.svg'],'svg');
%         ylim([0 1]);
    elseif i == 20
        ylabel('nf_hcn (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/mf_hcn.svg'],'svg');
    elseif i == 21
        ylabel('ns_hcn (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/msl_hcn.svg'],'svg');
    elseif i == 22
        ylabel('n_ka (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/a_ka.svg'],'svg');
    elseif i == 23
        ylabel('h - k_a (mA/cm^2)');
        legend boxoff 
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/b_ka.svg'],'svg');
    elseif i == 24
        ylabel('n_Kdr (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/m_Kdr.svg'],'svg');
    elseif i == 25
        ylabel('ns_km (mA/cm^2)');
        ylim([0 1]);
       %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/n_kmtype.svg'],'svg');
    elseif i == 26
        ylabel('nf_km (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/n_kmtype.svg'],'svg');
    elseif i == 27
        ylabel('w_knatype (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/w_knatype.svg'],'svg');
    elseif i == 28
        ylabel('m_nav1p8 (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/m_nav1p8.svg'],'svg');
    elseif i == 29
        ylabel('h_nav1p8 (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/h_nav1p8.svg'],'svg');
    elseif i == 30
        ylabel('s_nav1p8 (mA/cm^2)');
        ylim([0 1]);   
        %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/h_nav1p8.svg'],'svg');    elseif i == 30
    elseif i == 31
        ylabel('u_nav1p8 (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/h_nav1p8.svg'],'svg');       
    elseif i == 32     
        ylabel('m_nav1p9 (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/m_nav1p9.svg'],'svg');
     elseif i == 33
        ylabel('h_nav1p9(mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/h_nav1p9.svg'],'svg');
    elseif i == 34
        ylabel('s_nav1p9(mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/h_nav1p9.svg'],'svg'); 
    elseif i == 35
        ylabel('m_nav1p7 (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/O_nav17.svg'],'svg');
     elseif i == 36
        ylabel('h_nav1p7 (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/C1_nav17.svg'],'svg');
    elseif i == 37
        ylabel('s_nav1p7 (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/C1_nav17.svg'],'svg');       
    hold off 
    end
    end
    
      %% Params at pulse 
    tstart = 110;
    tstop = 140;
    start_test = find(t==tstart);
    end_test = find(t==tstop);
   
   c = jet(5);
    h = figure 
    hold on 
%     plot(t([start_test:end_test],1),axon{1,1}([start_test:end_test],1),'color',c(1,:));
%     plot(t([start_test:end_test],1),axon{1,5}([start_test:end_test],1),'color',c(2,:));
    plot(t([start_test:end_test],1),axon{1,14}([start_test:end_test],1),'color',c(1,:),'linewidth',0.5);
    plot(t([start_test:end_test],1),axon{1,16}([start_test:end_test],1),'color',c(2,:),'linewidth',0.5);
    plot(t([start_test:end_test],1),axon{1,18}([start_test:end_test],1),'color',c(3,:),'linewidth',0.5);
    plot(t([start_test:end_test],1),axon{1,20}([start_test:end_test],1),'color',c(4,:),'linewidth',0.5);
    plot(t([start_test:end_test],1),axon{1,22}([start_test:end_test],1),'color',c(5,:),'linewidth',0.5);
%     plot(t([start_test:end_test],1),axon{1,32}([start_test:end_test],1),'color',c(8,:));
%     plot(t([start_test:end_test],1),axon{1,36}([start_test:end_test],1),'color',c(9,:));
    axis([tstart tstop -120 80])
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northeast');
    legend boxoff 
%     legend({'node 1','node 5','node 14','node 16','node 18','node 20','node 22','node 32','node 36'}, 'Location','northeast');
    xlabel('Time, ms');
    ylabel('Membrane potential, mV');
    set(gcf, 'color','white')
    hold off 

    %saveas(h,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/vm.svg'],'svg');
    
    for i = 35:1:37
    k = figure 
    hold on 
%     plot(t([start_test:end_test],1),axon{1,1}([start_test:end_test],i),'color',c(1,:));
%     plot(t([start_test:end_test],1),axon{1,5}([start_test:end_test],i),'color',c(2,:));
    plot(t([start_test:end_test],1),axon{1,14}([start_test:end_test],i),'color',c(1,:),'linewidth',0.8);
    plot(t([start_test:end_test],1),axon{1,16}([start_test:end_test],i),'color',c(2,:),'linewidth',0.8);
    plot(t([start_test:end_test],1),axon{1,18}([start_test:end_test],i),'color',c(3,:),'linewidth',0.8);
    plot(t([start_test:end_test],1),axon{1,20}([start_test:end_test],i),'color',c(4,:),'linewidth',0.8);
    plot(t([start_test:end_test],1),axon{1,22}([start_test:end_test],i),'color',c(5,:),'linewidth',0.8);
%     plot(t([start_test:end_test],1),axon{1,32}([start_test:end_test],i),'color',c(8,:));
%     plot(t([start_test:end_test],1),axon{1,36}([start_test:end_test],i),'color',c(9,:));
%     legend({'node 1','node 5','node 14','node 16','node 18','node 20','node 22','node 32','node 36'}, 'Location','northwest');
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northwest');
    xlim([tstart tstop]);
    xlabel('Time, ms');
    set(gcf, 'color','white')
    if i == 2 
        ylabel('e_extracellular (mV)');
        %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/e_extracellular.svg'],'svg');
%         ylim([0 1]);
    elseif i == 3
        ylabel('K^+ current (mA/cm^2)');
        legend('Location','northeast');
        legend boxoff 
        xlim([110 130])
        %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ik.svg'],'svg');
%         ylim([0 1]);
    elseif i == 4
        ylabel('Na^+ current (mA/cm^2)');
        legend('Location','southeast');
        legend boxoff 
        xlim([110 130])
        %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ina.svg'],'svg');
%         ylim([0 1]);
    elseif i == 5
         ylabel('ki (mM)');
         %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ki.svg'],'svg');
         legend('Location','southwest');
%         ylim([0 1]);
    elseif i == 6
        ylabel('nai (mM)');
          % saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/nai.svg'],'svg');
           %lgd2.Location('northeast');
%         ylim([0 1]);
    elseif i == 7
        ylabel('ko (mM)');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ko.svg'],'svg');
%         ylim([0 1]);
    elseif i == 8
        ylabel('nao (mM)');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/nao.svg'],'svg');
%         ylim([0 1]);
    elseif i == 9
        ylabel('ina_hcn (mA/cm^2)');
        legend('Location','northeast');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ih_hcn.svg'],'svg');
%         ylim([0 1]);
    elseif i == 10
        ylabel('ik_hcn (mA/cm^2)');
        legend('Location','northeast');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ih_hcn.svg'],'svg');
%         ylim([0 1]);
    elseif i == 11
        ylabel('i_N_a_V_1_7 (mA/cm^2)');
        legend('Location','southeast');
        legend boxoff 
        xlim([110 130])
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ina_nav17.svg'],'svg');
%         ylim([0 1]);
    elseif i == 12
       ylabel('i_N_a_V_1_8 (mA/cm^2)');
        legend('Location','southwest');
        legend boxoff 
        xlim([110 130])
          % saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ina_nav18.svg'],'svg');
%         ylim([0 1]);
    elseif i == 13
        ylabel('ina_nav19 (mA/cm^2)');
        legend('Location','northeast');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ina_nav19.svg'],'svg');
%         ylim([0 1]);
    elseif i == 14
        ylabel('ina_NaKpump (mA/cm^2)');
        legend('Location','northeast');
          % saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ink_NaKpump.svg'],'svg');
%         ylim([0 1]);
    elseif i == 15
        ylabel('ik_NaKpump (mA/cm^2)');
        legend('Location','northeast');
          % saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ink_NaKpump.svg'],'svg');
%         ylim([0 1]);
    elseif i == 16
        ylabel('ik_knatype (mA/cm^2)');
        legend('Location','northeast');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ik_knatype.svg'],'svg');
%         ylim([0 1]);
    elseif i == 17
        ylabel('ik_km (mA/cm^2)');
        legend('Location','northeast');
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ik_kmtype.svg'],'svg');
%         ylim([0 1]);
    elseif i == 18
        ylabel('i_K_d_r (mA/cm^2)');
        legend('Location','northeast');
        legend boxoff 
        xlim([110 130])
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ik_Kdr.svg'],'svg');
%         ylim([0 1]);
    elseif i == 19
        ylabel('i_K_a (mA/cm^2)');
        legend('Location','northeast');
        legend boxoff 
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/jika_ka.svg'],'svg');
%         ylim([0 1]);
    elseif i == 20
        ylabel('nf_hcn (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/mf_hcn.svg'],'svg');
    elseif i == 21
        ylabel('ns_hcn (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/msl_hcn.svg'],'svg');
    elseif i == 22
        ylabel('n - k_a (mA/cm^2)');
        legend boxoff
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/a_ka.svg'],'svg');
    elseif i == 23
        ylabel('h - k_a (mA/cm^2)');
        legend boxoff 
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/b_ka.svg'],'svg');
    elseif i == 24
        ylabel('n - K_d_r (mA/cm^2)');
        legend('Location','northeast'); 
        legend boxoff 
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/m_Kdr.svg'],'svg');
    elseif i == 25
        ylabel('ns_km (mA/cm^2)');
        ylim([0 1]);
       %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/n_kmtype.svg'],'svg');
    elseif i == 26
        ylabel('nf_km (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/n_kmtype.svg'],'svg');
    elseif i == 27
        ylabel('w_knatype (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/w_knatype.svg'],'svg');
    elseif i == 28
        ylabel('m - NaV1.8 (mA/cm^2)');
        ylim([0 1]);
        legend('Location','northeast');
        legend boxoff 
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/m_nav1p8.svg'],'svg');
    elseif i == 29
        ylabel('h - NaV1.8 (mA/cm^2)');
        ylim([0 1]);
        legend('Location','southeast');
        legend boxoff 
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/h_nav1p8.svg'],'svg');
    elseif i == 30
        ylabel('s - NaV1.8 (mA/cm^2)');
        ylim([0 1]);  
        legend boxoff 
        %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/h_nav1p8.svg'],'svg');    elseif i == 30
    elseif i == 31
        ylabel('u - NaV1.8 (mA/cm^2)');
        ylim([0 1]);
        legend boxoff 
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/h_nav1p8.svg'],'svg');       
    elseif i == 32     
        ylabel('m_nav1p9 (mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/m_nav1p9.svg'],'svg');
     elseif i == 33
        ylabel('h_nav1p9(mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/h_nav1p9.svg'],'svg');
    elseif i == 34
        ylabel('s_nav1p9(mA/cm^2)');
        ylim([0 1]);
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/h_nav1p9.svg'],'svg'); 
    elseif i == 35
        ylabel('m - NaV1.7(mA/cm^2)');
        ylim([0 1]);
        legend boxoff 
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/O_nav17.svg'],'svg');
     elseif i == 36
        ylabel('h - NaV1.7 (mA/cm^2)');
        ylim([0 1]);
        legend boxoff 
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/C1_nav17.svg'],'svg');
    elseif i == 37
        ylabel('s - NaV1.7 (mA/cm^2)');
        ylim([0 1]);
        legend boxoff 
           %saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/C1_nav17.svg'],'svg');       
    hold off 
    end
    end
 %% Params over time at TEST PULSE 
    tstart = 95;
    tstop = 125;
    start_test = find(t==tstart);
    end_test = find(t==tstop);
    
    c = jet(9);
    h = figure 
    hold on 
    plot(t([start_test:end_test],1),axon{1,1}([start_test:end_test],1),'color',c(1,:));
    plot(t([start_test:end_test],1),axon{1,5}([start_test:end_test],1),'color',c(2,:));
    plot(t([start_test:end_test],1),axon{1,14}([start_test:end_test],1),'color',c(3,:));
    plot(t([start_test:end_test],1),axon{1,16}([start_test:end_test],1),'color',c(4,:));
    plot(t([start_test:end_test],1),axon{1,18}([start_test:end_test],1),'color',c(5,:));
    plot(t([start_test:end_test],1),axon{1,20}([start_test:end_test],1),'color',c(6,:));
    plot(t([start_test:end_test],1),axon{1,22}([start_test:end_test],1),'color',c(7,:));
    plot(t([start_test:end_test],1),axon{1,32}([start_test:end_test],1),'color',c(8,:));
    plot(t([start_test:end_test],1),axon{1,36}([start_test:end_test],1),'color',c(9,:));
    axis([tstart tstop -120 40])
    lgd = legend({'node 1','node 5','node 14','node 16','node 18','node 20','node 22','node 32','node 36'}, 'Location','northeast');
    lgd.FontSize = 7;
    xlabel('Time, ms');
    ylabel('Membrane potential, mV');
    set(gcf, 'color','white')
    hold off 
    
    saveas(h,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOMvm.svg'],'svg');
    
    for i = 2:1:31
    k = figure 
    hold on 
    plot(t([start_test:end_test],1),axon{1,1}([start_test:end_test],i),'color',c(1,:));
    plot(t([start_test:end_test],1),axon{1,5}([start_test:end_test],i),'color',c(2,:));
    plot(t([start_test:end_test],1),axon{1,14}([start_test:end_test],i),'color',c(3,:));
    plot(t([start_test:end_test],1),axon{1,16}([start_test:end_test],i),'color',c(4,:));
    plot(t([start_test:end_test],1),axon{1,18}([start_test:end_test],i),'color',c(5,:));
    plot(t([start_test:end_test],1),axon{1,20}([start_test:end_test],i),'color',c(6,:));
    plot(t([start_test:end_test],1),axon{1,22}([start_test:end_test],i),'color',c(7,:));
    plot(t([start_test:end_test],1),axon{1,32}([start_test:end_test],i),'color',c(8,:));
    plot(t([start_test:end_test],1),axon{1,36}([start_test:end_test],i),'color',c(9,:));
    lgd2 = legend({'node 1','node 5','node 14','node 16','node 18','node 20','node 22','node 32','node 36'}, 'Location','northwest');
    lgd2.FontSize = 7;
    xlim([tstart tstop]);
    xlabel('Time, ms');
    set(gcf, 'color','white')
    if i == 2 
        ylabel('e_extracellular (mV)');
        saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_e_extracellular.svg'],'svg');
%         ylim([0 1]);
    elseif i == 3
        ylabel('ik (mA/cm^2)');
        legend('Location','northeast');
        saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_ik.svg'],'svg');
%         ylim([0 1]);
    elseif i == 4
        ylabel('ina (mA/cm^2) ');
        legend('Location','southeast');
        saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_ina.svg'],'svg');
%         ylim([0 1]);
    elseif i == 5
         ylabel('ki (mM)');
         saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_ki.svg'],'svg');
         legend('Location','southwest');
%         ylim([0 1]);
    elseif i == 6
        ylabel('nai (mM)');
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_nai.svg'],'svg');
           %lgd2.Location('northeast');
%         ylim([0 1]);
    elseif i == 7
        ylabel('ko (mM)');
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_ko.svg'],'svg');
%         ylim([0 1]);
    elseif i == 8
        ylabel('nao (mM)');
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_nao.svg'],'svg');
%         ylim([0 1]);
    elseif i == 9
        ylabel('ih_hcn (mA/cm^2)');
        legend('Location','northeast');
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_ih_hcn.svg'],'svg');
%         ylim([0 1]);
    elseif i == 10
        ylabel('ina_nav17 (mA/cm^2)');
        legend('Location','southeast');
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_ina_nav17.svg'],'svg');
%         ylim([0 1]);
    elseif i == 11
        ylabel('ina_nav18 (mA/cm^2)');
        legend('Location','southwest');
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_ina_nav18.svg'],'svg');
%         ylim([0 1]);
    elseif i == 12
        ylabel('ina_nav19 (mA/cm^2)');
        legend('Location','northeast');
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_ina_nav19.svg'],'svg');;
%         ylim([0 1]);
    elseif i == 13
        ylabel('ink_NaKpump (mA/cm^2)');
        legend('Location','northeast');
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_ink_NaKpump.svg'],'svg');
%         ylim([0 1]);
    elseif i == 14
        ylabel('ik_knatype (mA/cm^2)');
        legend('Location','northeast');
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_ik_knatype.svg'],'svg');
%         ylim([0 1]);
    elseif i == 15
        ylabel('ik_kmtype (mA/cm^2)');
        legend('Location','northeast');
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_ik_kmtype.svg'],'svg');
%         ylim([0 1]);
    elseif i == 16
        ylabel('ik_Kdr (mA/cm^2)');
        legend('Location','northeast');
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_ik_Kdr.svg'],'svg');
%         ylim([0 1]);
    elseif i == 17
        ylabel('jika_ka (mA/cm^2)');
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_jika_ka.svg'],'svg');
%         ylim([0 1]);
    elseif i == 18
        ylabel('mf_hcn (mA/cm^2)');
        ylim([0 1]);
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_mf_hcn.svg'],'svg');
    elseif i == 19
        ylabel('msl_hcn (mA/cm^2)');
        ylim([0 1]);
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_msl_hcn.svg'],'svg');
    elseif i == 20
        ylabel('a_ka (mA/cm^2)');
        ylim([0 1]);
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_a_ka.svg'],'svg');
    elseif i == 21
        ylabel('b_ka (mA/cm^2)');
        ylim([0 1]);
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_b_ka.svg'],'svg');
    elseif i == 22
        ylabel('m_Kdr (mA/cm^2)');
        ylim([0 1]);
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_m_Kdr.svg'],'svg');
    elseif i == 23
        ylabel('h_Kdr (mA/cm^2)');
        ylim([0 1]);
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_h_Kdr.svg'],'svg');
    elseif i == 24
        ylabel('n_kmtype (mA/cm^2)');
        ylim([0 1]);
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_n_kmtype.svg'],'svg');
    elseif i == 25
        ylabel('w_knatype (mA/cm^2)');
        ylim([0 1]);
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_w_knatype.svg'],'svg');
    elseif i == 26
        ylabel('m_nav1p8 (mA/cm^2)');
        ylim([0 1]);
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_m_nav1p8.svg'],'svg');
    elseif i == 27
        ylabel('h_nav1p8 (mA/cm^2)');
        ylim([0 1]);
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_h_nav1p8.svg'],'svg');
    elseif i == 28
        ylabel('m_nav1p9 (mA/cm^2)');
        ylim([0 1]);
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_m_nav1p9.svg'],'svg');
     elseif i == 29
        ylabel('h_nav1p9(mA/cm^2)');
        ylim([0 1]);
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_h_nav1p9.svg'],'svg');
     elseif i == 30
        ylabel('O_nav17 (mA/cm^2)');
        ylim([0 1]);
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_O_nav17.svg'],'svg');
     elseif i == 31
        ylabel('C1_nav17 (mA/cm^2)');
        ylim([0 1]);
           saveas(k,[pwd,'/Figs/NOTEST_PARAMS_SINE_2um_1mmIED_11000kHz_-57mA/ZOOM_C1_nav17.svg'],'svg');
    hold off 
    end
    end
    
%% Params over time specifically at TEST PULSE 
    tstart = 95;
    tstop = 125;
    start_test = find(t==tstart);
    end_test = find(t==tstop);
    
    c = jet(9);
    h = figure 
    hold on 
    plot(t([start_test:end_test],1),axon{1,1}([start_test:end_test],1),'color',c(1,:));
    plot(t([start_test:end_test],1),axon{1,5}([start_test:end_test],1),'color',c(2,:));
    plot(t([start_test:end_test],1),axon{1,14}([start_test:end_test],1),'color',c(3,:));
    plot(t([start_test:end_test],1),axon{1,16}([start_test:end_test],1),'color',c(4,:));
    plot(t([start_test:end_test],1),axon{1,18}([start_test:end_test],1),'color',c(5,:));
    plot(t([start_test:end_test],1),axon{1,20}([start_test:end_test],1),'color',c(6,:));
    plot(t([start_test:end_test],1),axon{1,22}([start_test:end_test],1),'color',c(7,:));
    plot(t([start_test:end_test],1),axon{1,32}([start_test:end_test],1),'color',c(8,:));
    plot(t([start_test:end_test],1),axon{1,36}([start_test:end_test],1),'color',c(9,:));
    axis([tstart tstop -120 40])
    lgd = legend({'node 1','node 5','node 14','node 16','node 18','node 20','node 22','node 32','node 36'}, 'Location','northeast');
    lgd.FontSize = 7;
    xlabel('Time, ms');
    ylabel('NaKpump ina');
    set(gcf, 'color','white')
    hold off 
    
    saveas(h,[pwd,'/Figs/B_PARAMS_SINE_2um_1mmIED_11000kHz_-58mA/ZOOMnakpumpina.svg'],'svg');
    
    
    h2 = figure 
    hold on 
    plot(t([start_test:end_test],1),axon{1,1}([start_test:end_test],2),'color',c(1,:));
    plot(t([start_test:end_test],1),axon{1,5}([start_test:end_test],2),'color',c(2,:));
    plot(t([start_test:end_test],1),axon{1,14}([start_test:end_test],2),'color',c(3,:));
    plot(t([start_test:end_test],1),axon{1,16}([start_test:end_test],2),'color',c(4,:));
    plot(t([start_test:end_test],1),axon{1,18}([start_test:end_test],2),'color',c(5,:));
    plot(t([start_test:end_test],1),axon{1,20}([start_test:end_test],2),'color',c(6,:));
    plot(t([start_test:end_test],1),axon{1,22}([start_test:end_test],2),'color',c(7,:));
    plot(t([start_test:end_test],1),axon{1,32}([start_test:end_test],2),'color',c(8,:));
    plot(t([start_test:end_test],1),axon{1,36}([start_test:end_test],2),'color',c(9,:));
    axis([tstart tstop -120 40])
    lgd = legend({'node 1','node 5','node 14','node 16','node 18','node 20','node 22','node 32','node 36'}, 'Location','northeast');
    lgd.FontSize = 7;
    xlabel('Time, ms');
    ylabel('NaKpump ik');
    set(gcf, 'color','white')
    hold off 
    
    saveas(h2,[pwd,'/Figs/B_PARAMS_SINE_2um_1mmIED_11000kHz_-58mA/ZOOMnakpumpik.svg'],'svg');
    
  
   


    