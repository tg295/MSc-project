clc
close all
clear all
%% Loading data and saving into struct elements

Amplitudes = [ ];
Frequencies = [ ];
nseg = [1:36];
Data = struct ();

files_list=dir('*.txt');

for i=1:length(files_list)
    
    f=fopen(files_list(i).name,'r');
    headers = textscan(f, '%s %f');
    fclose(f);
    
    params=cell2mat(headers(2));
    
    Data(i).dt = params(1);      %expressed in milliseconds
    Data(i).tfinal = params(2);
    Data(i).acfreq = params(3);
    Data(i).amplitude = -params(4);
    Data(i).nseg = params(5);
    Data(i).L = params(6);
    Data(i).diam = params(7);
    Data(i).stim_amp = params(8);
    
    tmp=importdata(files_list(i).name,' ',length(params)+1);
    
    Data(i).T= 1/Data(i).acfreq * 1000 ; %ms
    Data(i).nT = round(Data(i).T / Data(i).dt); %samples per period
    Data(i).total_num_of_T = round(Data(i).tfinal / Data(i).T);
    Data(i).t = [0 : length(tmp.data)-1] / ( length(tmp.data) - 1 ) * Data(i).tfinal;
    
%     Data(i).vm05 = tmp.data(:,1);
%     Data(i).e_ext = tmp.data(:,2);
%     Data(i).m = tmp.data(:,3);
%     Data(i).h = tmp.data(:,4);
%     Data(i).n = tmp.data(:,5);
%     Data(i).ina = tmp.data(:,6);
%     Data(i).ik = tmp.data(:,7);

    Data(i).vm05 = tmp.data(:,1);
    Data(i).m = tmp.data(:,2);
    Data(i).s = tmp.data(:,3);
    Data(i).h = tmp.data(:,4);
    Data(i).mp = tmp.data(:,5);
    Data(i).inap = tmp.data(:,6);

%    % loop over the nodes
%    for j = 1 : length(nseg)
%        % State variables
%        Data(i).nodes(j).node = j;
%        
%        %Voltages
%        
%        Data(i).nodes(j).vm = tmp.data(:,j+1);
%         Data(i).nodes(j).ve = tmp.data(:,j+36);
%         Data(i).nodes(j).vi = Data(i).nodes(j).vm + Data(i).nodes(j).ve;
%        
%    end
%     
%     Data(i).nodes(j+1).vm = tmp.data(:,end);
%     Data(i).nodes(j+1).ve = tmp.data(:,end);
%     Data(i).nodes(j+1).vi = Data(i).nodes(j+1).vm + Data(i).nodes(j+1).ve;
    
    
    if ~ any( Amplitudes == Data(i).amplitude)
        Amplitudes = [ Amplitudes Data(i).amplitude ];
    end
    
    if ~ any( Frequencies == Data(i).acfreq)
        Frequencies = [ Frequencies Data(i).acfreq ];
    end
    
end

Frequencies = sort(Frequencies);
Amplitudes = sort(Amplitudes);

%% Sort structure 

Afields = fieldnames(Data);
Acell = struct2cell(Data);
sz = size(Acell)
% Convert to a matrix
Acell = reshape(Acell, sz(1), []);      % Px(MxN)

% Make each field a column
Acell = Acell';                         % (MxN)xP

% Sort by first field "name"
Acell = sortrows(Acell, 3)

% Put back into original cell array format
Acell = reshape(Acell', sz);

% Convert to Struct
Asorted = cell2struct(Acell, Afields, 1);

for id = 1:length(Asorted)
    fprintf('%d\n',id)
    disp(Asorted(id))
end

%% Plot activation variables for all amps

linS = {'-','--',':','-.','-'};
%c=gray(5);
%{
c=[0,0,0;0.1,0.1,0.1;0.2,0.2,0.2;0.3,0.3,0.3;0.4,0.4,0.4]
figure  
hold on
for i = 1:length(Asorted)
    plot(Asorted(i).t([33000:35000]), Asorted(i).m([33000:35000]),'color',c(i,:),'linestyle',linS{i},'linewidth',1);
end 
legend({'60mA','70mA','80mA','90mA','100mA'})
ylabel('Na+ Channel Activation (m)');
xlim([33 35])
ylim([0 1])
hold off 
%}
  st = 34.5;
  et = 40;
  start_test = find(Asorted(1).t==st);
  end_test = find(Asorted(1).t==et);
    
d=parula(5);
figure  
hold on
for i = 1:length(Asorted)
    plot(Asorted(i).t(:), Asorted(i).m(:),'color',d(i,:),'linewidth',1);
end 
legend({'1.3mA','1.4mA','1.5mA','1.6mA','1.7mA'})
% xlim([33 35])
% xticks([33 34 35])
xlabel('Time (ms)','fontsize',14);
ylim([0 1])
ylabel('m','fontsize',20);
legend boxoff
hold off 

figure  
hold on
for i = 1:length(Asorted)
    plot(Asorted(i).t(:), Asorted(i).s(:),'color',d(i,:),'linewidth',1);
end 
legend({'1.3mA','1.4mA','1.5mA','1.6mA','1.7mA'})
xlabel('Time (ms)','fontsize',14);
ylabel('s','fontsize',20);
% xlim([32 40])
ylim([0 1])
legend boxoff
hold off

figure
hold on
for i = 1:length(Asorted)
    plot(Asorted(i).t(:), Asorted(i).h(:),'color',d(i,:),'linewidth',1);
end 
legend({'1.3mA','1.4mA','1.5mA','1.6mA','1.7mA'})
xlabel('Time (ms)','fontsize',14);
ylabel('h','fontsize',20);
%xlim([32 40])
ylim([0 1])
legend boxoff
hold off

figure
hold on
for i = 1:length(Asorted)
    plot(Asorted(i).t(:), Asorted(i).mp(:),'color',d(i,:),'linewidth',1);
end 
legend({'1.3mA','1.4mA','1.5mA','1.6mA','1.7mA'})
xlabel('Time (ms)','fontsize',14);
ylabel('mp','fontsize',20);
%xlim([32 40])
ylim([0 1])
legend boxoff
hold off


figure
hold on
for i = 1:length(Asorted)
    plot(Asorted(i).t(:), Asorted(i).inap(:),'color',d(i,:),'linewidth',1);
end 
legend({'1.3mA','1.4mA','1.5mA','1.6mA','1.7mA'})
xlabel('Time (ms)','fontsize',14);
ylabel('inap','fontsize',20);
ylim([-0.00000001 0.00000001])
legend boxoff
hold off

%% Plot activation variables for all amps

linS = {'-','--',':','-.','-'};
%c=gray(5);
%{
c=[0,0,0;0.1,0.1,0.1;0.2,0.2,0.2;0.3,0.3,0.3;0.4,0.4,0.4]
figure  
hold on
for i = 1:length(Asorted)
    plot(Asorted(i).t([33000:35000]), Asorted(i).m([33000:35000]),'color',c(i,:),'linestyle',linS{i},'linewidth',1);
end 
legend({'60mA','70mA','80mA','90mA','100mA'})
ylabel('Na+ Channel Activation (m)');
xlim([33 35])
ylim([0 1])
hold off 
%}
  st = 35;
  et = 37;
  start_test = find(Asorted(1).t==st);
  end_test = find(Asorted(1).t==et);
    
d=parula(5);
figure  
hold on
for i = 1:length(Asorted)
    plot(Asorted(i).t([start_test:end_test]), Asorted(i).m([start_test:end_test]),'color',d(i,:),'linewidth',1);
end 
legend({'1.3mA','1.4mA','1.5mA','1.6mA','1.7mA'},'location','southeast')
% xlim([33 35])
% xticks([33 34 35])
xlabel('Time (ms)','fontsize',14);
ylim([0 1])
ylabel('m','fontsize',20);
legend boxoff
hold off 

figure  
hold on
for i = 1:length(Asorted)
    plot(Asorted(i).t([start_test:end_test]), Asorted(i).s([start_test:end_test]),'color',d(i,:),'linewidth',1);
end 
legend({'1.3mA','1.4mA','1.5mA','1.6mA','1.7mA'})
xlabel('Time (ms)','fontsize',14);
ylabel('s','fontsize',20);
% xlim([32 40])
ylim([0.8 0.9])
legend boxoff
hold off

figure
hold on
for i = 1:length(Asorted)
    plot(Asorted(i).t([start_test:end_test]), Asorted(i).h([start_test:end_test]),'color',d(i,:),'linewidth',1);
end 
legend({'1.3mA','1.4mA','1.5mA','1.6mA','1.7mA'},'location','southeast')
xlabel('Time (ms)','fontsize',14);
ylabel('h','fontsize',20);
%xlim([32 40])
ylim([0 0.1])
legend boxoff
hold off

figure
hold on
for i = 1:length(Asorted)
    plot(Asorted(i).t([start_test:end_test]), Asorted(i).mp([start_test:end_test]),'color',d(i,:),'linewidth',1);
end 
legend({'1.3mA','1.4mA','1.5mA','1.6mA','1.7mA'})
xlabel('Time (ms)','fontsize',14);
ylabel('mp','fontsize',20);
%xlim([32 40])
ylim([0 1])
legend boxoff
hold off

figure
hold on
for i = 1:length(Asorted)
    plot(Asorted(i).t([start_test:end_test]), Asorted(i).inap([start_test:end_test]),'color',d(i,:),'linewidth',1);
end 
legend({'1.3mA','1.4mA','1.5mA','1.6mA','1.7mA'})
xlabel('Time (ms)','fontsize',14);
ylabel('inap','fontsize',20);
%xlim([32 40])
ylim([0 1])
legend boxoff
hold off
%% Plot block and transmission regions

steady_state = 40; %ms
stim_delay = 50; %ms

block_region = zeros(length(Amplitudes), length(Frequencies));
trans_region = zeros(length(Amplitudes), length(Frequencies));
fire_region = zeros(length(Amplitudes), length(Frequencies));

for i=1:length(Frequencies)
    
    for j=1:length(Amplitudes)
        
        ind = find([Data.acfreq] == Frequencies(i) & [Data.amplitude] == Amplitudes(j));
        
        nT = Data(ind).nT;
        total_Ts=  Data(ind).total_num_of_T;
        t = Data(ind).t;
        k = round(steady_state /  Data(ind).tfinal *  total_Ts);
        
        vm1 = Data(ind).nodes(5).vm;
%         vm1 = tsmovavg(tsmovavg(vm1, 's',nT,1), 's',nT,1);
        
        vm2 = Data(ind).nodes(36).vm;
%         vm2 = tsmovavg(tsmovavg(vm2, 's',nT,1), 's',nT,1);
        
        [pks1, locs1] = findpeaks(vm1, 'MinPeakHeight',0, 'MinPeakDistance',100);
        [pks2, locs2] = findpeaks(vm2, 'MinPeakHeight',-20, 'MinPeakDistance',100);

        if sum(t(locs2)>stim_delay) == 0 & (sum(t(locs1)<stim_delay & t(locs1)>stim_delay-10) == 0) ==1
            block_region(j, i) = 1; % JUST SAYING IF THEYRE ARE NO SPIKES ARTER STIM AND NO SPIKES BETWEEN 20 AND 30 (ONSET OF TEST STIM) THEN IT HAS BLOCKED 
        end
        
        if (sum(t(locs2)>stim_delay) )> 0 | (sum(t(locs1)<stim_delay & t(locs1)>stim_delay-10)> 0) == 1
            trans_region(j, i) = 1; %OPPOSITE
        end
        
       % if sum(t(locs2)<stim_delay) > 1 & sum(t(locs2)>stim_delay) > 1
       %     fire_region(j, i) = 1; %FIRE REGION IF THERE ARE SPIKES BEFORE AND AFTER ONSET OF STIM
        %end
        
    end
    
end

%Check block or transmission
% block_region = avg_dep > -70;
block_region = bwboundaries(block_region);
trans_region = bwboundaries(trans_region);
fire_region = bwboundaries(fire_region);

f=figure;
hold on
for k=1:length(block_region)
    b = fill(Frequencies(block_region{k,1}(:,2)) * 0.001 ,Amplitudes(block_region{k,1}(:,1)),'b','EdgeColor','blue','FaceColor','none','LineWidth',2);
end

for k=1:length(trans_region)
    trans = fill(Frequencies(trans_region{k,1}(:,2)) * 0.001 ,Amplitudes(trans_region{k,1}(:,1)),'r','EdgeColor','r','FaceColor','none','LineWidth',2);
end

%for k=1:length(fire_region)
%    fire = fill(Frequencies(fire_region{k,1}(:,2)) * 0.001 ,Amplitudes(fire_region{k,1}(:,1)),'r','EdgeColor','m','FaceColor','none','LineWidth',1);
%end

axis([min(Frequencies)/1000 max(Frequencies)/1000 min(Amplitudes) max(Amplitudes) ] )
set(gca, 'ydir','reverse');
% title(['Fiber ' num '\mum, \delta 1mm'])
xlabel('Frequency of the external eletric field, kHz')
ylabel('Amplitude of AC stimulation, mA')
hold off

%% Plotting map of average depolarization under the electrode at steady state

avg_dep = [];
steady_state = 40; %ms
stim_delay = 50; %ms
% node = find(Sections == 304);
% fiberD = 5.7;%um

block_region = zeros(length(Amplitudes), length(Frequencies));
trans_region = zeros(length(Amplitudes), length(Frequencies));
fire_region = zeros(length(Amplitudes), length(Frequencies));

for i=1:length(Frequencies)
    
    for j=1:length(Amplitudes)
        
        ind = find([Data.acfreq] == Frequencies(i) & [Data.amplitude] == Amplitudes(j));
        
        nT = Data(ind).nT;
        total_Ts=  Data(ind).total_num_of_T;
        t = Data(ind).t;
        k = round(steady_state /  Data(ind).tfinal *  total_Ts);
        
        
%         node1 = find(Sections == 0);
%         node2 = find(Sections == 594);
        
        vm1 = Data(ind).nodes(5).vm;
%         vm1 = tsmovavg(tsmovavg(vm1, 's',nT,1), 's',nT,1);
        
        vm2 = Data(ind).nodes(30).vm;
%         vm2 = tsmovavg(tsmovavg(vm2, 's',nT,1), 's',nT,1);
        
        [pks1, locs1] = findpeaks(vm1, 'MinPeakHeight',0, 'MinPeakDistance',100);
        [pks2, locs2] = findpeaks(vm2, 'MinPeakHeight',-20, 'MinPeakDistance',100);

        if sum(t(locs2)>stim_delay) == 0 %& (sum(t(locs1)<stim_delay & t(locs1)>stim_delay-10) == 0) ==1
            block_region(j, i) = 1; % JUST SAYING IF THEYRE ARE NO SPIKES ARTER STIM AT END NODE AND NO SPIKES BETWEEN 20 AND 30 (ONSET OF TEST STIM) THEN IT HAS BLOCKED 
        end
        %if the sum of spikes between 10 secs before onset of stim and stim
        %is equal to zero in the first node 
        if (sum(t(locs2)>stim_delay) )> 0 %| (sum(t(locs1)<stim_delay & t(locs1)>stim_delay-10)> 0) == 1
            trans_region(j, i) = 1; %OPPOSITE
        end
        
        if sum(t(locs2)<stim_delay & t(locs2)>stim_delay) > 0
            fire_region(j, i) = 1; %FIRE REGION IF THERE ARE SPIKES BEFORE AND AFTER ONSET OF STIM
        end
        
        avg_dep(j, i) = mean(Data(ind).nodes(18).vm((k-1)*nT+1:k*nT)) ;
        
    end
    
end

%Check block or transmission
% block_region = avg_dep > -70;
block_region = bwboundaries(block_region);
trans_region = bwboundaries(trans_region);
fire_region = bwboundaries(fire_region);

%Plotting
f=figure;
avg_dep = avg_dep +80;
surf(Frequencies * 0.001, Amplitudes, avg_dep);
view(2)
shading interp
set(gca, 'ydir','reverse');
axis([min(Frequencies)/1000 max(Frequencies)/1000 min(Amplitudes) max(Amplitudes) ] )
% title(['Fiber ' num '\mum, \delta 1mm'])
xlabel('Frequency of the external eletric field, kHz')
ylabel('Amplitude of AC stimulation, mA')
m = colorbar;

hold on

for k=1:length(block_region)
    b = fill3(Frequencies(block_region{k,1}(:,2)) * 0.001 ,Amplitudes(block_region{k,1}(:,1)),ones(size(block_region{k,1}(:,2))).*max(max(avg_dep)),'b','EdgeColor','blue','FaceColor','none','LineWidth',2);
end

for k=1:length(trans_region)
    trans = fill3(Frequencies(trans_region{k,1}(:,2)) * 0.001 ,Amplitudes(trans_region{k,1}(:,1)),ones(size(trans_region{k,1}(:,2))).*max(max(avg_dep)),'r','EdgeColor','r','FaceColor','none','LineWidth',2);
end

%{
for k=1:length(block_region)
    b = fill2(Frequencies(block_region{k,1}(:,2)) * 0.001 ,Amplitudes(block_region{k,1}(:,1)),'b','EdgeColor','blue','FaceColor','none','LineWidth',2);
end

for k=1:length(trans_region)
    trans = fill2(Frequencies(trans_region{k,1}(:,2)) * 0.001 ,Amplitudes(trans_region{k,1}(:,1)),'EdgeColor','r','FaceColor','none','LineWidth',2);
end
%}
% for k=1:length(fire_region)
%     fire = fill3(Frequencies(fire_region{k,1}(:,2)) * 0.001 ,Amplitudes(fire_region{k,1}(:,1)),ones(size(fire_region{k,1}(:,2))).*max(max(avg_dep)),'r','EdgeColor','b','FaceColor','none');
% end

% text(mean(b.XData), mean(b.YData), max(max(avg_dep)),'Block','Color','black','FontSize',14)
% text(mean(trans.XData), mean(trans.YData), max(max(avg_dep)),'Transmission','Color','black','FontSize',14)

ylabel(m, 'V_m - Resting potential, mV')
