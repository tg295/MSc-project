close all 
clear all 
clc
%% READ FILE PARAMETERS AND DATA MATRIX
var = 7; %number of variables being looked at
% Open the file.
fullFileName = fullfile(pwd,'Sine','2um','1mm IED','negative','10000_Hz_-88.0_mA.txt');
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
ffname = fullfile(pwd,'Sine','1.5um','1mm IED','10000_Hz_0.0_mA_ST.txt');
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
    
    for i = 1:18
    axondata(:,i) = axon{1,2*i}(:,1);
    end 
    
    specaxondata = zeros(20002,9);
    specaxondata = axondata([50000:60000,70000:80000],:);
    
    %% 3D plot of Vm over axon length 
    figure
    hold on 
    length = L/1000;
    y = 0:0.5:length-1;
    x = 1:20002;
    z = specaxondata(:,[1:17]);
    waterfall(x', y', z'); 
    axis([0 20000 0 8 -80 50]);
     set(gca,...
    'XTickLabel',{'0','5','10','15','20'},...
    'XTick',[0 5000 10000 15000 20000],...
    'YTick',[1:8],...
    'ZTick',[],...
    'Zcolor','none');
%     xlabel('Time, ms')
%     ylabel('Length along axon, mm')
    set(gcf, 'color','white')
    view(-47,68)
    box off
    hold off 

    %% Vm over all nodes 
    mycolours = parula(nseg);
    figure 
    hold on;
    for i = 1:nseg
    plot(t([69500:75000],1),axon{1,i}([69500:75000],1),'color',mycolours(i,:));
    end 
    colorbar('Ticks',[0,1],'TickLabels',{'Node 1','Node 36'})
    hold off;
    
    %% Params over nodes close to AC stim
    c = jet(5);
    figure 
    hold on 
    plot(t([71000:76000],1),axon{1,14}([71000:76000],1),'color',c(1,:));
    plot(t([71000:76000],1),axon{1,16}([71000:76000],1),'color',c(2,:));
    plot(t([71000:76000],1),axon{1,18}([71000:76000],1),'color',c(3,:));
    plot(t([71000:76000],1),axon{1,20}([71000:76000],1),'color',c(4,:));
    plot(t([71000:76000],1),axon{1,22}([71000:76000],1),'color',c(5,:));
    axis([71 76 -120 40])
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northwest');
    xlabel('Time, ms');
    ylabel('Membrane potential, mV');
    set(gcf, 'color','white')
    hold off 
    
    for i = 3:1:7
    figure 
    hold on 
    plot(t([71000:76000],1),axon{1,14}([71000:76000],i),'color',c(1,:));
    plot(t([71000:76000],1),axon{1,16}([71000:76000],i),'color',c(2,:));
    plot(t([71000:76000],1),axon{1,18}([71000:76000],i),'color',c(3,:));
    plot(t([71000:76000],1),axon{1,20}([71000:76000],i),'color',c(4,:));
    plot(t([71000:76000],1),axon{1,22}([71000:76000],i),'color',c(5,:));
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northwest');
    xlim([71 76]);
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
    
     %% Params over nodes close to AC stim OVER WHOLE TIME 
    c = jet(5);
    figure 
    hold on 
    plot(t(:,1),axon{1,14}(:,1),'color',c(1,:));
    plot(t(:,1),axon{1,16}(:,1),'color',c(2,:));
    plot(t(:,1),axon{1,18}(:,1),'color',c(3,:));
    plot(t(:,1),axon{1,20}(:,1),'color',c(4,:));
    plot(t(:,1),axon{1,22}(:,1),'color',c(5,:));
    axis([50 tfinal -120 40])
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northwest');
    xlabel('Time, ms');
    ylabel('Membrane potential, mV');
    set(gcf, 'color','white')
    hold off 
    
    for i = 3:1:7
    figure 
    hold on 
    plot(t(:,1),axon{1,14}(:,i),'color',c(1,:));
    plot(t(:,1),axon{1,16}(:,i),'color',c(2,:));
    plot(t(:,1),axon{1,18}(:,i),'color',c(3,:));
    plot(t(:,1),axon{1,20}(:,i),'color',c(4,:));
    plot(t(:,1),axon{1,22}(:,i),'color',c(5,:));
    legend({'node 14','node 16','node 18','node 20','node 22'}, 'Location','northwest');
    xlim([50 tfinal]);
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
   %% Params over nodes when test stim occurs 
    c = jet(9);
    figure 
    hold on 
    plot(t(:,1),axon{1,1}(:,1),'color',c(1,:));
    plot(t(:,1),axon{1,5}(:,1),'color',c(2,:));
    plot(t(:,1),axon{1,14}(:,1),'color',c(3,:));
    plot(t(:,1),axon{1,16}(:,1),'color',c(4,:));
    plot(t(:,1),axon{1,18}(:,1),'color',c(5,:));
    plot(t(:,1),axon{1,20}(:,1),'color',c(6,:));
    plot(t(:,1),axon{1,22}(:,1),'color',c(7,:));
    plot(t(:,1),axon{1,32}(:,1),'color',c(8,:));
    plot(t(:,1),axon{1,36}(:,1),'color',c(9,:));
    axis([0 tfinal -120 40])
    legend({'node 1','node 5','node 14','node 16','node 18','node 20','node 22','node 32','node 36'}, 'Location','northwest');
    xlabel('Time, ms');
    ylabel('Membrane potential, mV');
    set(gcf, 'color','white')
    hold off 
    
    for i = 2:1:17
    figure 
    hold on 
    plot(t(:,1),axon{1,1}(:,i),'color',c(1,:));
    plot(t(:,1),axon{1,5}(:,i),'color',c(2,:));
    plot(t(:,1),axon{1,14}(:,i),'color',c(3,:));
    plot(t(:,1),axon{1,16}(:,i),'color',c(4,:));
    plot(t(:,1),axon{1,18}(:,i),'color',c(5,:));
    plot(t(:,1),axon{1,20}(:,i),'color',c(6,:));
    plot(t(:,1),axon{1,22}(:,i),'color',c(7,:));
    plot(t(:,1),axon{1,32}(:,i),'color',c(8,:));
    plot(t(:,1),axon{1,36}(:,i),'color',c(9,:));
    legend({'node 11','node 5','node 14','node 16','node 18','node 20','node 22','node 32','node 36'}, 'Location','northwest');
    xlim([0 tfinal]);
    xlabel('Time, ms');
    set(gcf, 'color','white')
    if i == 3 
        ylabel('Na+ Channel Activation (m)');
%         ylim([0 1]);
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
  

    