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
        
        %Voltages
        Data(i).params.vm = tmp.data(:,1);
        Data(i).params.e_ext = tmp.data(:,2);
        Data(i).params.vmEND = tmp.data(:,38);
%         Data(i).nodes(j).vi = Data(i).nodes(j).vm + Data(i).nodes(j).ve;

        %Currents 
        Data(i).params.ik = tmp.data(:,3);
        Data(i).params.ina = tmp.data(:,4);
        Data(i).params.ina_hcn = tmp.data(:,9);
        Data(i).params.ik_hcn = tmp.data(:,10);
        Data(i).params.ina_nav1p7 = tmp.data(:,11);
        Data(i).params.ina_nav1p8 = tmp.data(:,12);
        Data(i).params.ina_nav1p9 = tmp.data(:,13);
        Data(i).params.ina_NaKpump = tmp.data(:,14);
        Data(i).params.ik_NaKpump = tmp.data(:,15);
        Data(i).params.ik_knatype = tmp.data(:,16);
        Data(i).params.ik_km = tmp.data(:,17);
        Data(i).params.ik_kdr = tmp.data(:,18);
        Data(i).params.ik_ka = tmp.data(:,19);
        
        %Concentrations 
        Data(i).params.ki = tmp.data(:,5);
        Data(i).params.nai = tmp.data(:,6);
        Data(i).params.ko = tmp.data(:,7);
        Data(i).params.nao = tmp.data(:,8);
        
        %Gating variables 
        Data(i).params.nf_hcn = tmp.data(:,20);
        Data(i).params.hs_hcn = tmp.data(:,21);
        Data(i).params.n_ka = tmp.data(:,22);
        Data(i).params.h_ka = tmp.data(:,23);
        Data(i).params.n_kdr = tmp.data(:,24);
        Data(i).params.ns_km = tmp.data(:,25);
        Data(i).params.nf_km = tmp.data(:,26);
        Data(i).params.w_knatype = tmp.data(:,27);
        Data(i).params.m_nav1p8 = tmp.data(:,28);
        Data(i).params.h_nav1p8 = tmp.data(:,29);
        Data(i).params.s_nav1p8 = tmp.data(:,30);
        Data(i).params.u_nav1p8 = tmp.data(:,31);
        Data(i).params.m_nav1p9 = tmp.data(:,32);
        Data(i).params.h_nav1p9 = tmp.data(:,33);
        Data(i).params.s_nav1p9 = tmp.data(:,34);
        Data(i).params.m_nav1p7 = tmp.data(:,35);
        Data(i).params.h_nav1p7 = tmp.data(:,36);
        Data(i).params.s_nav1p7 = tmp.data(:,37);
         
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

t = Data(1).t;
%% Compute block threshold
steady_state = 80; %ms
stim_delay = 100; %ms

block_region = zeros(length(Amplitudes), 1);
trans_region = zeros(length(Amplitudes), 1);

for i = 1:length(Amplitudes)
    
     ind = find([Data.amplitude] == Amplitudes(i));
        
        nT = Data(ind).nT;
        total_Ts=  Data(ind).total_num_of_T;
        t = Data(ind).t;
        k = round(steady_state /  Data(ind).tfinal *  total_Ts);
        
        vmEND = Data(ind).params.vmEND;
        
    [pks1, locs1] = findpeaks(vmEND, 'MinPeakHeight',0, 'MinPeakDistance',100);
        
        if sum(t(locs1)>stim_delay) == 0 
            block_region(i) = 1; 
        end
       
        if (sum(t(locs1)>stim_delay) )> 0 
            trans_region(i) = 1; 
        end
        
end 

BLOCK = [Amplitudes; block_region']

%% Currents 
%ik
c = jet(5)
C1 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ik(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ik, (mA/cm^2)');
hold off
%fname = 'D:\Thesis\data\Unmyelin_LDiff2\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ik.svg';

%ina
C2 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ina(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ina, (mA/cm^2)');
hold off
%ina_hcn
C3 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ina_hcn(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ina hcn, (mA/cm^2)');
hold off
%ik_hcn
C4 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ik_hcn(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ik hcn, (mA/cm^2)');
hold off
%ina_nav1.7
C5 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ina_nav1p7(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ina nav1.7, (mA/cm^2)');
hold off
%ina_nav1.8
C6 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ina_nav1p8(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ina nav1.8, (mA/cm^2)');
hold off
%ina_nav1.9
C7 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ina_nav1p9(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ina nav1.9, (mA/cm^2)');
hold off
%ina_NaKpump
C8 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ina_NaKpump(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ina pump, (mA/cm^2)');
hold off
%ik_NaKpump 
C9 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ik_NaKpump(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ik pump, (mA/cm^2)');
hold off
%ik_knatype
C10 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ik_knatype(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ik kna, (mA/cm^2)');
hold off
%ik_km
C11 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ik_km(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ik km, (mA/cm^2)');
hold off
%ik_kdr 
C12 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ik_kdr(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ik kdr, (mA/cm^2)');
hold off
%ik_ka
C13 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ik_ka(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ik ka, (mA/cm^2)');
hold off
%SAVEALL
% saveas(C1, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ik.svg'],'svg');
% saveas(C2, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ina.svg'],'svg');
% saveas(C3, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ina_hcn.svg'],'svg');
% saveas(C4, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ik_hcn.svg'],'svg');
% saveas(C5, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ina_nav1p7.svg'],'svg');
% saveas(C6, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ina_nav1p8.svg'],'svg');
% saveas(C7, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ina_nav1p9.svg'],'svg');
% saveas(C8, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ina_NaKpump.svg'],'svg');
% saveas(C9, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ik_NaKpump.svg'],'svg');
% saveas(C10, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ik_knatype.svg'],'svg');
% saveas(C11, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ik_km.svg'],'svg');
% saveas(C12, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ik_kdr.svg'],'svg');
% saveas(C13, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ik_ka.svg'],'svg');

%% Concentrations 
% ki
Co1 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ki(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ki, (mM)');
hold off
% ko
Co2 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ko(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('ko, (mM)');
hold off
% nai
Co3 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.nai(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('nai, (mM)');
hold off
% nao 
Co4 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.nao(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
xlim([110 150])
xlabel('Time, ms');
ylabel('nao, (mM)');
hold off
%SAVEALL
% saveas(Co1, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ki.svg'],'svg');
% saveas(Co2, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ko.svg'],'svg');
% saveas(Co3, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\nai.svg'],'svg');
% saveas(Co4, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\nao.svg'],'svg');
%% K+ gating variables
%n_kdr 
c = jet(5)
G1 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.n_kdr(find(t==100):find(t==150)),'color',c(i,:),'linewidth',1);
end 
legend('80mA','90mA','100mA','110mA','120mA')
legend boxoff
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('n - K_d_r');
hold off
%n_ka
G2 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.n_ka(find(t==100):find(t==150)),'color',c(i,:),'linewidth',1);
end 
legend('80mA','90mA','100mA','110mA','120mA')
legend boxoff
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('n - K_a');
hold off
%h_ka
G3 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.h_ka(find(t==100):find(t==150)),'color',c(i,:),'linewidth',1);
end 
legend('80mA','90mA','100mA','110mA','120mA')
legend boxoff
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('h - K_a');
hold off
%ns_km
G4 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.ns_km(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('ns km');
hold off
%nf_km 
G5 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.nf_km(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('nf km');
hold off
%ns_hcn
G6 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.hs_hcn(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('ns hcn');
hold off
%nf_hcn
G7 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.nf_hcn(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('nf hcn');
hold off
%w_knatype
G8 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.w_knatype(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('w kna');
hold off
%SAVEALL
% saveas(G1, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\n_kdr.svg'],'svg');
% saveas(G2, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\n_ka.svg'],'svg');
% saveas(G3, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\h_ka.svg'],'svg');
% saveas(G4, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ns_km.svg'],'svg');
% saveas(G5, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\nf_km.svg'],'svg');
% saveas(G6, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\ns_hcn.svg'],'svg');
% saveas(G7, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\nf_hcn.svg'],'svg');
% saveas(G8, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\w_knatype.svg'],'svg');

%% Na+ gating variables 
%m_nav1.7
G9 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==130)), Asorted(i).params.m_nav1p7(find(t==100):find(t==130)),'color',c(i,:),'linewidth',0.7);
end 
legend('80mA','90mA','100mA','110mA','120mA')
legend boxoff
axis([110.01 130 0 1])
xlabel('Time, ms');
ylabel('m - NaV1.7');
hold off

%h_nav1.7
G10 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.h_nav1p7(find(t==100):find(t==150)),'color',c(i,:),'linewidth',0.7);
end 
legend('80mA','90mA','100mA','110mA','120mA')
legend boxoff
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('h - NaV1.7');
hold off

%s_nav1.7
G11 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.s_nav1p7(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('s nav1p7');
hold off

%m_nav1.8
G12 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==130)), Asorted(i).params.m_nav1p8(find(t==100):find(t==130)),'color',c(i,:),'linewidth',0.7);
end 
legend('80mA','90mA','100mA','110mA','120mA')
legend boxoff
axis([110.01 130 0 1])
xlabel('Time, ms');
ylabel('m - NaV1.8');
hold off

%h_nav1.8
G13 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==140)), Asorted(i).params.h_nav1p8(find(t==100):find(t==140)),'color',c(i,:),'linewidth',0.7);
end 
legend('80mA','90mA','100mA','110mA','120mA')
legend boxoff
legend('location','southeast')
axis([110.01 140 0 1])
xlabel('Time, ms');
ylabel('h- NaV1.8');
hold off

%s_nav1.8
G14 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.s_nav1p8(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('s nav1p8');
hold off

%u_nav1.8
G15 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.u_nav1p8(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('u nav1p8');
hold off

%m_nav1.9
G16 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.m_nav1p9(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('m nav1p9');
hold off

%h_nav1.9
G17 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.h_nav1p9(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('h nav1p9');
hold off

%s_nav1.9
G18 = figure
hold on 
for i = find(Frequencies == 7000):find(Frequencies == 11000)
    plot(t(find(t==100):find(t==150)), Asorted(i).params.s_nav1p9(find(t==100):find(t==150)),'color',c(i,:));
end 
legend('80mA','90mA','100mA','110mA','120mA')
axis([110.01 150 0 1])
xlabel('Time, ms');
ylabel('s nav1p9');
hold off

%SAVEALL
% saveas(G9, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\m_nav1p7.svg'],'svg');
% saveas(G10, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\h_nav1p7.svg'],'svg');
% saveas(G11, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\s_nav1p7.svg'],'svg');
% saveas(G12, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\m_nav1p8.svg'],'svg');
% saveas(G13, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\h_nav1p8.svg'],'svg');
% saveas(G14, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\s_nav1p8.svg'],'svg');
% saveas(G15, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\u_nav1p8.svg'],'svg');
% saveas(G16, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\m_nav1p9.svg'],'svg');
% saveas(G17, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\h_nav1p9.svg'],'svg');
% saveas(G18, [pwd,'\Figs\PARAMS_SINE_2um_1mmIED_CONST_11kHz\s_nav1p9.svg'],'svg');

