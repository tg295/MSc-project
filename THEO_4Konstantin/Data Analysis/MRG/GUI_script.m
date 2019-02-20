close all 
clear all 
clc
%% READ FILE PARAMETERS AND DATA MATRIX
nnodes = 36;
var = 15; %number of variables being looked at
% Open the file.
fullFileName = fullfile(pwd,'Sine','5.7um','1mm IED','10000_Hz_-1.6_mA.txt');
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
	elseif contains(textLine(1:6), 'fiberD', 'IgnoreCase', true)
		nseg = str2double(textLine(8:end))
	elseif contains(textLine(1:2), 'L ', 'IgnoreCase', true)
		L = str2double(textLine(3:end))
	elseif contains(textLine(1:8), 'stim_amp ', 'IgnoreCase', true)
		stim_amp = str2double(textLine(10:end))
	elseif contains(textLine, 'batch_run from', 'IgnoreCase', true)
		% Now, after this batch_run line we're expecting 5 lines with 5 numbers per line.
		m = zeros(tfinal/dt+1, nnodes*var); % Preallocate space for speed.
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
axon = mat2cell(m, [tfinal/dt+1], repmat(var, 1, nnodes));
t = [0:dt:tfinal]';

%Calculating im 
for i = 1:nnodes    
i_m(:,i) = axon{1,i}(:,11) + axon{1,i}(:,12) + axon{1,i}(:,13) + axon{1,i}(:,14) + axon{1,i}(:,15);
end 
%Calculating iNa + iK
for i = 1:nnodes    
i_nak(:,i) = axon{1,i}(:,11) + axon{1,i}(:,12);
end 
%% CALC CONDUCTION VELOCITY
%first read spike time file
L = 35000;
nnodes = 36;
ffname = fullfile(pwd,'Sine','8.7um','1mm IED','10000_Hz_-1.0_mA_ST.txt');
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

C = zeros(1,length(B));
D = zeros(1,length(B));
for i=1:length(B)-1
    C(1,i) = find(B(:,i)~=0, 1, 'first');
end 
for i=1:length(B)-1
D(1,i) = B(C(1,i),i);
end 
% Remove first 3 segments because test stim occurs in the 4th
E = D(:,4:length(B));
%% calculate conduction velocity 
deltax = L/nnodes;
deltat = zeros(1,length(E));
cond_vel = zeros(1,length(E));

for i = 1:length(E)-1
    deltat(1,i) = E(1,i+1)-E(1,i);
end 

for i = 1:length(E)-1
    cond_vel(1,i) = deltax/deltat(1,i)/1000
end 

   %% Reorganise for plotting 
  start_onset = find(t==5);
  end_onset = find(t==15);
  start_test = find(t==35);
  end_test = find(t==45);
  rows_r = (end_test-start_test)+(end_onset-start_onset);
  
    for i = 1:18
    axondata(:,i) = axon{1,2*i}(:,1);
    end 
    
    specaxondata = zeros(rows_r,18);
    specaxondata = axondata([start_onset:end_onset,start_test:end_test],:);
    %specaxondata = axondata([1000:2000],[6000:8000],:);
    
    %% 3D plot of Vm over axon length 
    figure
    hold on 
    l = L/1000;
    y = 0:l;
    x = 1:4002;
    z = specaxondata(:,[1:18]);
    waterfall(x', y', z'); 
    axis([0 rows_r 0 18 -80 50]);
     set(gca,...
    'XTickLabel',{'0','5','10','15','20'},...
    'XTick',[0  1000  2000 3000 4000],...
    'YTick',[0:2:18],...
    'ZTick',[],...
    'Zcolor','none');
%     xlabel('Time, ms')
%     ylabel('Length along axon, mm')
    set(gcf, 'color','white')
    view(-40,65)
    box off
    hold off 
    
    %% Vm over all nodes 
    mycolours = parula(nnodes);
    figure 
    hold on;
    for i = 1:nnodes
    plot(t(:,1),axon{1,i}(:,1),'color',mycolours(i,:));
    end 
    colorbar('Ticks',[0,1],'TickLabels',{'Node 1','Node 36'})
    hold off;
    
    %% Plotting sodium current and potassium current over all time
    c = jet(5);
    for i = 11:12 
    figure 
    hold on 
    plot(t(:,1),axon{1,14}(:,i),'color',c(1,:));
    plot(t(:,1),axon{1,16}(:,i),'color',c(2,:));
    plot(t(:,1),axon{1,18}(:,i),'color',c(3,:));
    plot(t(:,1),axon{1,20}(:,i),'color',c(4,:));
    plot(t(:,1),axon{1,22}(:,i),'color',c(5,:));
    xlim([0 50]);
     
    xlabel('Time, ms');
    if i == 11
    ylabel('Sodium current (mA/cm^2)');
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','southeast');
    elseif i == 12 
    ylabel('Potassium current (mA/cm^2)');
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northwest');
    end
    set(gcf, 'color','white')
    hold off 
    end 
    %% Params over nodes close to AC stim
    st = 34.5;
    et = 37;
    start_test = find(t==st);
    end_test = find(t==et);
    
    %Plotting Vm
    c = jet(5);
    figure 
    hold on 
    plot(t([start_test:end_test],1),axon{1,14}([start_test:end_test],1),'color',c(1,:),'linewidth',1);
    plot(t([start_test:end_test],1),axon{1,16}([start_test:end_test],1),'color',c(2,:),'linewidth',1);
    plot(t([start_test:end_test],1),axon{1,18}([start_test:end_test],1),'color',c(3,:),'linewidth',1);
    plot(t([start_test:end_test],1),axon{1,20}([start_test:end_test],1),'color',c(4,:),'linewidth',1);
    plot(t([start_test:end_test],1),axon{1,22}([start_test:end_test],1),'color',c(5,:),'linewidth',1);
    axis([st et -120 40])
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northwest');
    xlabel('Time (ms)','fontsize',14);
    ylabel('V_m (mV)','fontsize',20);
    set(gcf, 'color','white')
    hold off 

    % i = 3 : m_axnode
    % i = 4 : s_axnode
    % i = 5 : h_axnode 
    % i = 6 : mp_axnode
    % i = 7 : tau_m_axnode
    % i = 8 : tau_mp_axnode
    % i = 9 : tau_h_axnode
    % i = 10 : tau_s_axnode
    % i = 11 : ina_axnode
    % i = 12 : ik_axnode
    % i = 13 : inap_axnode 
    % i = 14 : il_axnode 
    % i = 15 : i_cap
    
    %Plotting channel activations
    for i = [3:6]
    figure 
    hold on 
    plot(t([start_test:end_test],1),axon{1,14}([start_test:end_test],i),'color',c(1,:),'linewidth',1);
    plot(t([start_test:end_test],1),axon{1,16}([start_test:end_test],i),'color',c(2,:),'linewidth',1);
    plot(t([start_test:end_test],1),axon{1,18}([start_test:end_test],i),'color',c(3,:),'linewidth',1);
    plot(t([start_test:end_test],1),axon{1,20}([start_test:end_test],i),'color',c(4,:),'linewidth',1);
    plot(t([start_test:end_test],1),axon{1,22}([start_test:end_test],i),'color',c(5,:),'linewidth',1);
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northwest');
    xlim([st et]);
    xlabel('Time (ms)','fontsize',14);
    set(gcf, 'color','white')
    if i == 3 
        ylabel('m','fontsize',20);
        ylim([0 1]);
        legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northeast');
    elseif i == 4
        ylabel('s','fontsize',20);
        ylim([0 1]);
        legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','southeast');
    elseif i == 5
        ylabel('h','fontsize',20);
        ylim([0 1]);
    elseif i == 6
        ylabel('m_p','fontsize',20);
        ylim([0 1]);
        legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','southeast');
         hold off 
    end
    end
  
    %Plotting currents 
    
    for i = [11:15] 
    figure 
    hold on 
    plot(t([start_test:end_test],1),axon{1,14}([start_test:end_test],i),'color',c(1,:),'linewidth',1);
    plot(t([start_test:end_test],1),axon{1,16}([start_test:end_test],i),'color',c(2,:),'linewidth',1);
    plot(t([start_test:end_test],1),axon{1,18}([start_test:end_test],i),'color',c(3,:),'linewidth',1);
    plot(t([start_test:end_test],1),axon{1,20}([start_test:end_test],i),'color',c(4,:),'linewidth',1);
    plot(t([start_test:end_test],1),axon{1,22}([start_test:end_test],i),'color',c(5,:),'linewidth',1);
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northwest');
    xlim([st et]);
    xlabel('Time (ms)','fontsize',14);
    set(gcf, 'color','white')
    if i == 11 
        ylabel('i_N_a_f  (mA/cm^2)','fontsize',20);
        %ylim([0 1]);
         legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','southeast');
    elseif i == 12
        ylabel('i_K_s  (mA/cm^2)','fontsize',20);
        %ylim([0 1]);
         legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northeast');
    elseif i == 13
        ylabel('i_N_a_p  (mA/cm^2)','fontsize',20);
        %ylim([0 1]);
    elseif i == 14
        ylabel('Leakage current (mA/cm^2)','fontsize',20);
        %ylim([0 1]);
    elseif i == 15
        ylabel('Capacitative current (mA/cm^2)');
        %ylim([0 1]);
         legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northeast');
         hold off 
    end
    end
     
    %Total Im 
    figure 
    hold on 
    plot(t([start_test:end_test],1),i_m([start_test:end_test],14),'color',c(1,:),'linewidth',1);
    plot(t([start_test:end_test],1),i_m([start_test:end_test],16),'color',c(2,:),'linewidth',1);
    plot(t([start_test:end_test],1),i_m([start_test:end_test],18),'color',c(3,:),'linewidth',1);
    plot(t([start_test:end_test],1),i_m([start_test:end_test],20),'color',c(4,:),'linewidth',1);
    plot(t([start_test:end_test],1),i_m([start_test:end_test],22),'color',c(5,:),'linewidth',1);
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','southeast');
    xlim([st et]);
    xlabel('Time (ms)','fontsize',14);
    ylabel('Total membrane current (mA/cm^2)'); 
    set(gcf, 'color','white')
    hold off
 %{   
    % iNa + iK 
     figure 
    hold on 
    plot(t([start_test:end_test],1),i_nak([start_test:end_test],14),'color',c(1,:));
    plot(t([start_test:end_test],1),i_nak([start_test:end_test],16),'color',c(2,:));
    plot(t([start_test:end_test],1),i_nak([start_test:end_test],18),'color',c(3,:));
    plot(t([start_test:end_test],1),i_nak([start_test:end_test],20),'color',c(4,:));
    plot(t([start_test:end_test],1),i_nak([start_test:end_test],22),'color',c(5,:));
     legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','southeast');
    xlim([st et]);
    xlabel('Time, ms');
    ylabel('i_{Na} + i_K (mA/cm^2)'); 
    set(gcf, 'color','white')
    hold off
%}
     %% Params over nodes close to AC stim OVER WHOLE TIME
    %Plotting Vm
    c = jet(5);
    figure 
    hold on 
    plot(t(:,1),axon{1,14}(:,1),'color',c(1,:));
    plot(t(:,1),axon{1,16}(:,1),'color',c(2,:));
    plot(t(:,1),axon{1,18}(:,1),'color',c(3,:));
    plot(t(:,1),axon{1,20}(:,1),'color',c(4,:));
    plot(t(:,1),axon{1,22}(:,1),'color',c(5,:));
    axis([0 tfinal -120 40])
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northwest');
    xlabel('Time, ms');
    ylabel('Membrane potential, mV');
    set(gcf, 'color','white')
    hold off 

    % i = 3 : m_axnode
    % i = 4 : s_axnode
    % i = 5 : h_axnode 
    % i = 6 : mp_axnode
    % i = 7 : tau_m_axnode
    % i = 8 : tau_mp_axnode
    % i = 9 : tau_h_axnode
    % i = 10 : tau_s_axnode
    % i = 11 : ina_axnode
    % i = 12 : ik_axnode
    % i = 13 : inap_axnode 
    % i = 14 : il_axnode 
    % i = 15 : i_cap
    
    %Plotting channel activations
    for i = [3:6]
    figure 
    hold on 
    plot(t(:,1),axon{1,14}(:,i),'color',c(1,:));
    plot(t(:,1),axon{1,16}(:,i),'color',c(2,:));
    plot(t(:,1),axon{1,18}(:,i),'color',c(3,:));
    plot(t(:,1),axon{1,20}(:,i),'color',c(4,:));
    plot(t(:,1),axon{1,22}(:,i),'color',c(5,:));
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northwest');
    xlim([0 tfinal]);
    xlabel('Time, ms');
    set(gcf, 'color','white')
    if i == 3 
        ylabel('Fast Na+ Channel Activation (m)');
        ylim([0 1]);
        legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northeast');
    elseif i == 4
        ylabel('Slow K+ Channel Activation (s)');
        ylim([0 1]);
        legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','southeast');
    elseif i == 5
        ylabel('Fast Na+ Channel Inactivation (h)');
        ylim([0 1]);
    elseif i == 6
        ylabel('Persisent Na+ Channel Activation (mp)');
        ylim([0 1]);
        legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','southeast');
         hold off 
    end
    end
  
    %Plotting currents 
    
    for i = [11:15] 
    figure 
    hold on 
    plot(t(:,1),axon{1,14}(:,i),'color',c(1,:));
    plot(t(:,1),axon{1,16}(:,i),'color',c(2,:));
    plot(t(:,1),axon{1,18}(:,i),'color',c(3,:));
    plot(t(:,1),axon{1,20}(:,i),'color',c(4,:));
    plot(t(:,1),axon{1,22}(:,i),'color',c(5,:));
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northwest');
    xlim([0 tfinal]);
    xlabel('Time, ms');
    set(gcf, 'color','white')
    if i == 11 
        ylabel('Fast Na+ current (mA/cm^2)');
        %ylim([0 1]);
         legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','southeast');
    elseif i == 12
        ylabel('Slow K+ current (mA/cm^2)');
        %ylim([0 1]);
         legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northeast');
    elseif i == 13
        ylabel('Persistent Na+ current (mA/cm^2)');
        %ylim([0 1]);
    elseif i == 14
        ylabel('Leakage current (mA/cm^2)');
        %ylim([0 1]);
    elseif i == 15
        ylabel('Capacitative current (mA/cm^2)');
        %ylim([0 1]);
         legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northeast');
         hold off 
    end
    end
     
    %Total Im 
    figure 
    hold on 
    plot(t(:,1),i_m(:,14),'color',c(1,:));
    plot(t(:,1),i_m(:,16),'color',c(2,:));
    plot(t(:,1),i_m(:,18),'color',c(3,:));
    plot(t(:,1),i_m(:,20),'color',c(4,:));
    plot(t(:,1),i_m(:,22),'color',c(5,:));
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','southeast');
    xlim([0 tfinal]);
    xlabel('Time, ms');
    ylabel('Total membrane current (mA/cm^2)'); 
    set(gcf, 'color','white')
    hold off
    
    % iNa + iK 
     figure 
    hold on 
    plot(t(:,1),i_nak(:,14),'color',c(1,:));
    plot(t(:,1),i_nak(:,16),'color',c(2,:));
    plot(t(:,1),i_nak(:,18),'color',c(3,:));
    plot(t(:,1),i_nak(:,20),'color',c(4,:));
    plot(t(:,1),i_nak(:,22),'color',c(5,:));
     legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','southeast');
    xlim([0 tfinal]);
    xlabel('Time, ms');
    ylabel('i_{Na} + i_K (mA/cm^2)'); 
    set(gcf, 'color','white')
    hold off

   %% Params over nodes when test stim occurs 
    c = jet(5);
    figure 
    hold on 
    plot(t([10000:20000],1),axon{1,14}([10000:20000],1),'color',c(1,:));
    plot(t([10000:20000],1),axon{1,16}([10000:20000],1),'color',c(2,:));
    plot(t([10000:20000],1),axon{1,18}([10000:20000],1),'color',c(3,:));
    plot(t([10000:20000],1),axon{1,20}([10000:20000],1),'color',c(4,:));
    plot(t([10000:20000],1),axon{1,22}([10000:20000],1),'color',c(5,:));
    axis([11 16 -120 40])
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northwest');
    xlabel('Time, ms');
    ylabel('Membrane potential, mV');
    set(gcf, 'color','white')
    hold off 
    
    for i = 3:1:7
    figure 
    hold on 
    plot(t([10000:20000],1),axon{1,14}([10000:20000],i),'color',c(1,:));
    plot(t([10000:20000],1),axon{1,16}([10000:20000],i),'color',c(2,:));
    plot(t([10000:20000],1),axon{1,18}([10000:20000],i),'color',c(3,:));
    plot(t([10000:20000],1),axon{1,20}([10000:20000],i),'color',c(4,:));
    plot(t([10000:20000],1),axon{1,22}([10000:20000],i),'color',c(5,:));
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northwest');
    xlim([11 16]);
    xlabel('Time, ms');
    set(gcf, 'color','white')
    if i == 3 
        ylabel('Na+ Channel Activation (m)');
        ylim([0 1]);
    elseif i == 4
        ylabel('Na+ Channel Inactivation (h)');
        ylim([0 1]);
    elseif i == 5
        ylabel('K+ Channel Activation (n)');
        ylim([0 1]);
    elseif i == 6
        ylabel('Na+ Current (mA/cm^2)');
        legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','southeast');
%         ylim([0 1]);
    elseif i == 7
        ylabel('K+ Current (mA/cm^2)');
%         ylim([0 1]);
    hold off 
    end
    end 
 


    