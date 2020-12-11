%Electrophysiology manipulation
%   This script is to be used for the loading, cleaning, filtering, and
%processing of elecrophysiological data. It is suggested that this code be
%run section by section as each section take a lot of time to complete and
%may crash Matlab if run all at once. Run sections IN ORDER as they are
%dependent on the sections before it.
%
%   This script is built for cases where stimulus presentations include CS+
%and CS- trials. Data gained from Blackrock must be placed into a .mat file
% including the event time points and only the .NS6 data from utilized electrodes.
%
%Assumptions: SinoFilt_joe, mtcsg, BlackRock NPMK package
%Version:1
%Created on:4/28/2018
%Created by: Morgan Harnois (morgan.harnois@gmail.com)

%% Formatting files
%This is used if more channels were being recorded from than the ones being
%used. Alter the channels listed in open NSx if necessary to be the
%channels utilized (add 2 to the number due to blackrock numbering). Check
%that they are the correct channels before running the save function.

% openNSx('L:\Morgan\IL\Ephys\Mrgn11\FC_Recall\FC_Recall.ns4','c:1,2,3')
% save('FC_Recall.mat','NS4','-append')

%% 0 Initialization

clear all
close all
clc

Preprocessdir = 'C:\Users\Elikhtik\Desktop\Morgan\'; 
Processeddir = 'E:\Morgan\Processed\';
dirname=Processeddir; 
Animaldir=dir([dirname,'Mrgn*']);
%load('IL_ephys_workspace.mat') %This will automatically load the data post-processing in the final structures for manipulation and figure creation. If you load this then skip to section 9!

%% 1 Seperating ephys data by trial and filtering
%Setting the directory and finding all animal files within the directory 

%This for loop will cut each electrode recording for each mouse into the trials and filter it
for animals=1:length(Animaldir)
    clear LFPs
    try
    A=dir([dirname,'Itnik*']);
    %filename=[dirname,Animaldir(animals).name,'\',A.name];
    filename=[dirname,A(animals).name];
    openNEV([filename,'\FC_Recall\FC_Recall.nev'],'noread','nosave','overwrite')
    load([filename,'\FC_Recall\FC_Recall.mat'])
    load([filename,'\FC_Recall\CSorder.mat'])
    openNSx([filename,'\FC_Recall\FC_Recall.ns6'])
    nCSCs=length(NS4.ElectrodesInfo);
    switch nCSCs
        case 3
            names={'BLA','BF','PFC'};
        case 6
            names={'rBLA','BLA','rBF','BF','rPFC','PFC'};
        case 4
            names={'BLA','rBF','BF','PFC'};
    end
    ds_rate=NS4.MetaTags.SamplingFreq/2000;
    %Your ds_rate is dependent on the type of recording... your goal is to
    %downsample so that you are showing 2000 samples per second, which is
    %an NS3. Samples/second determined through your blackrock output.
    %NS1=500, NS2=1000, NS3=2000, NS4=10000, NS5 = 30000, NS6=30000 (raw)
    for a=1:length(names)
        LFPs.(names{a})=downsample(double(NS4.Data(a,:)),ds_rate);
    end
    LFPs.ts=linspace(1,NS4.MetaTags.DataPoints/ds_rate,NS4.MetaTags.DataPoints/ds_rate);
    LFPs.filename=Animaldir(animals).name;
    
    trialstart=NEV.Data.SerialDigitalIO.TimeStampSec(NEV.Data.SerialDigitalIO.UnparsedData==65503);%65487
    trialstop=NEV.Data.SerialDigitalIO.TimeStampSec(NEV.Data.SerialDigitalIO.UnparsedData==65535);
    
    LFPs.OnTimes=trialstart;
    LFPs.OffTimes=trialstart + 30;  %LFPs.OffTimes=trialstop;
    
    for a=1:length(LFPs.OnTimes)
        start = LFPs.OnTimes(a)*(NS4.MetaTags.SamplingFreq/ds_rate);
        precsStart = LFPs.OnTimes(a)*(NS4.MetaTags.SamplingFreq/ds_rate)-(30*(NS4.MetaTags.SamplingFreq/ds_rate));
        stop = LFPs.OffTimes(a)*(NS4.MetaTags.SamplingFreq/ds_rate);
        precsStop = LFPs.OnTimes(a)*(NS4.MetaTags.SamplingFreq/ds_rate);
        LFPs.CSts{a}= LFPs.ts(1,start:stop);
        LFPs.pCSts{a}=LFPs.ts(1,precsStart:precsStop);
    end
  
    for b=1:length(names)
        for c=1:length(LFPs.CSts)
            ind=cell2mat(LFPs.CSts(c));
            ind(1);
            ind(length(ind));
            [LFPs.Filt.(names{b}).CSP{c},foo,LFPs.Pwr.(names{b}).CSP{c}]=SinoFilt_Joe((LFPs.(names{b})(ind(1):ind(length(ind)))),2000,9,4,12);
            [LFPs.Filt.(names{b}).pCSP{c},foo,LFPs.Pwr.(names{b}).preCSP{c}]=SinoFilt_Joe((LFPs.(names{b})(ind(1):ind(length(ind)))),2000,9,4,12);
        end
    end
  
    save([Postprocessdir,LFPs.filename],'LFPs','-v7.3');
    clear LFPs
end

%% 2 Sorting animals into experimental and control
%Here we change the directory into the newly created processed folder
%containing all the cut and filtered data
dirname=Processeddir;

LFPsall.eyfp=struct('Mrgn1',[],'Mrgn2',[],'Mrgn3',[],'Mrgn8',[],'Mrgn9',[],'Mrgn12',[]);
LFPsall.arch=struct('Mrgn4',[],'Mrgn5',[],'Mrgn10',[],'Mrgn11',[]);
LFPsall.Internal=struct('Mrgn6',[],'Mrgn7',[]);
numCS= 10; %zeros(length(Morgandir),1);

%This for loop splits the animals into experimental or control groups and 
%places them into the structure we will be pulling from for the rest of the code.
for i=1:length(Animaldir)
    A=dir([dirname,Animaldir(i).name]);
    filename=[dirname,Animaldir(i).name];
    load(filename)
    anim = i;
    switch anim
        case 1
        LFPsall.eyfp.(['Mrgn',num2str(i)])=LFPs;
        case 2
        LFPsall.eyfp.(['Mrgn',num2str(i)])=LFPs;
        case 3
        LFPsall.eyfp.(['Mrgn',num2str(i)])=LFPs;
        case 4
        LFPsall.arch.(['Mrgn',num2str(i)])=LFPs;
        case 5
        LFPsall.arch.(['Mrgn',num2str(i)])=LFPs;
        case 6
        LFPsall.Internal.(['Mrgn',num2str(i)])=LFPs;
        case 7
        LFPsall.Internal.(['Mrgn',num2str(i)])=LFPs;
        case 8
        LFPsall.eyfp.(['Mrgn',num2str(i)])=LFPs;
        case 9
        LFPsall.eyfp.(['Mrgn',num2str(i)])=LFPs;
        case 10
        LFPsall.arch.(['Mrgn',num2str(i)])=LFPs;
        case 11
        LFPsall.arch.(['Mrgn',num2str(i)])=LFPs;
        case 12
        LFPsall.eyfp.(['Mrgn',num2str(i)])=LFPs;
    end
        clear LFPs
end


%% 3 Basic visualization of raw data, spectrogram, and power
%Within this section we can do a basic visualization of the data to check
%for any errors before running analysis. 

names={'BLA','BF','PFC'};
anat = 1; %Region selected
TR = 1; %CS trial selected
fs=2000; %Sampling rate

l1=length(LFPsall.eyfp.Mrgn1.Filt.(names{anat}).CSP{TR});
l2=length(LFPsall.eyfp.Mrgn2.Filt.(names{anat}).CSP{TR});
l3=length(LFPsall.eyfp.Mrgn3.Filt.(names{anat}).CSP{TR});
l4=length(LFPsall.arch.Mrgn4.Filt.(names{anat}).CSP{TR});
l5=length(LFPsall.arch.Mrgn5.Filt.(names{anat}).CSP{TR});
l8=length(LFPsall.eyfp.Mrgn8.Filt.(names{anat}).CSP{TR});

y1= fft(LFPsall.eyfp.Mrgn1.Filt.(names{anat}).CSP{TR});
y2= fft(LFPsall.eyfp.Mrgn2.Filt.(names{anat}).CSP{TR});
y3= fft(LFPsall.eyfp.Mrgn3.Filt.(names{anat}).CSP{TR});
y4= fft(LFPsall.arch.Mrgn4.Filt.(names{anat}).CSP{TR});
y5= fft(LFPsall.arch.Mrgn5.Filt.(names{anat}).CSP{TR});
y8= fft(LFPsall.eyfp.Mrgn8.Filt.(names{anat}).CSP{TR});

p1=abs(y1/l1);
p12=p1(1:l1/2+1);
p12(2:end-1)=2*p12(2:end-1);
p2=abs(y2/l2);
p22=p2(1:l2/2+1);
p22(2:end-1)=2*p22(2:end-1);
p3=abs(y3/l3);
p32=p3(1:l3/2+1);
p32(2:end-1)=2*p32(2:end-1);
p4=abs(y4/l4);
p42=p4(1:l4/2+1);
p42(2:end-1)=2*p42(2:end-1);
p5=abs(y5/l5);
p52=p5(1:l5/2+1);
p52(2:end-1)=2*p52(2:end-1);
p8=abs(y8/l8);
p82=p8(1:l8/2+1);
p82(2:end-1)=2*p82(2:end-1);

f=fs*(0:(l4/2))/l4;

ySubDivisions = (0:1/fs:20); %This creates the y axis frequency subdivisions
wind = kaiser(2048); %This sets our spetrogram window as kaiser with a length of 2048
overl = length(wind)/2; %This sets our overlap as the maximum possible overlap
figure('color','w')
subplot(6,3,1)
plot(LFPsall.eyfp.Mrgn1.Filt.(names{anat}).CSP{TR})
xlabel('Time (sec)')
title('EYFP filtered signal (Mrgn1)')
subplot(6,3,2)
spectrogram(LFPsall.eyfp.Mrgn1.Filt.(names{anat}).CSP{TR},wind,overl,ySubDivisions,fs,'yaxis')
title('EYFP spectrogram (Mrgn1)')
subplot(6,3,3)
plot(f,p12,'k')
xlim([0 10])

subplot(6,3,4)
plot(LFPsall.eyfp.Mrgn2.Filt.(names{anat}).CSP{TR})
xlabel('Time (sec)')
title('EYFP filtered signal (Mrgn2)')
subplot(6,3,5)
spectrogram(LFPsall.eyfp.Mrgn2.Filt.(names{anat}).CSP{TR},wind,overl,ySubDivisions,fs,'yaxis')
title('EYFP spectrogram (Mrgn2)')
subplot(6,3,6)
plot(f,p22,'k')
xlim([0 10])

subplot(6,3,7)
plot(LFPsall.eyfp.Mrgn3.Filt.(names{anat}).CSP{TR})
xlabel('Time (sec)')
title('EYFP filtered signal (Mrgn3)')
subplot(6,3,8)
spectrogram(LFPsall.eyfp.Mrgn3.Filt.(names{anat}).CSP{TR},wind,overl,ySubDivisions,fs,'yaxis')
title('EYFP spectrogram (Mrgn3)')
subplot(6,3,9)
plot(f,p32,'k')
xlim([0 10])

subplot(6,3,10)
plot(LFPsall.eyfp.Mrgn8.Filt.(names{anat}).CSP{TR})
xlabel('Time (sec)')
title('EYFP filtered signal (Mrgn8)')
subplot(6,3,11)
spectrogram(LFPsall.eyfp.Mrgn8.Filt.(names{anat}).CSP{TR},wind,overl,ySubDivisions,fs,'yaxis')
title('EYFP spectrogram (Mrgn8)')
subplot(6,3,12)
plot(f,p82,'k')
xlim([0 10])

subplot(6,3,13)
plot(LFPsall.arch.Mrgn4.Filt.(names{anat}).CSP{TR})
title('ARCH filtered signal (Mrgn4)')
xlabel('Time (sec)')
subplot(6,3,14)
spectrogram(LFPsall.arch.Mrgn4.Filt.(names{anat}).CSP{TR},wind,overl,ySubDivisions,fs,'yaxis')
title('ARCH spectrogram (Mrgn4)')
subplot(6,3,15)
plot(f,p42,'k')
xlim([0 10])

subplot(6,3,16)
plot(LFPsall.arch.Mrgn5.Filt.(names{anat}).CSP{TR})
title('ARCH filtered signal (Mrgn5)')
xlabel('Time (sec)')
subplot(6,3,17)
spectrogram(LFPsall.arch.Mrgn5.Filt.(names{anat}).CSP{TR},wind,overl,ySubDivisions,fs,'yaxis')
title('ARCH spectrogram (Mrgn5)')
subplot(6,3,18)
plot(f,p52,'k')
xlim([0 10])


%% 4 EYFP Power and Coherence Calculation
%This section will calculate the power and coherence for each trial and
%pre-trial section for each EYFP animal for each electrode. This
%information is added into our LFP structure within each animals level.
names={'BLA','BF','PFC'};
animalse={'Mrgn1','Mrgn2','Mrgn3','Mrgn9','Mrgn12'};
% animalse={'Mrgn1','Mrgn2','Mrgn3','Mrgn12'};

for a=1:length(animalse)
    for b=1:length(names)
        for c=1:numCS
            [LFPsall.eyfp.(animalse{a}).TonePwr.(names{b}).CSp_freqtime{c},Freq,Time]=mtcsg(LFPsall.eyfp.(animalse{a}).(names{b})(floor(LFPsall.eyfp.(animalse{a}).CSts{c})),4096,2000,2048,1024,2);%2048,2000,500,250,2 
            [LFPsall.eyfp.(animalse{a}).TonePwr.(names{b}).pCSp_freqtime{c},foo1,foo2]=mtcsg(LFPsall.eyfp.(animalse{a}).(names{b})(floor(LFPsall.eyfp.(animalse{a}).pCSts{c})),4096,2000,2048,1024,2);
        end
    end
end

animalse=intersect(animalse_bf,animalse_pfc);
for a=1:length(animalse)
    for b=1:length(names)
        for c=1:numCS
%             [LFPsall.eyfp.(animalse{a}).ToneCoh.BLApfc.CSp_freqtime{c},Freq_c,Time_c]=mscohere(LFPsall.eyfp.(animalse{a}).(names{1})(floor(LFPsall.eyfp.(animalse{a}).CSts{c})),LFPsall.eyfp.(animalse{a}).(names{2})(floor(LFPsall.eyfp.(animalse{a}).CSts{c})),2048,1024,4096,2000);
%             [LFPsall.eyfp.(animalse{a}).ToneCoh.BLApfc.pCSp_freqtime{c},foo1,foo2]=mscohere(LFPsall.eyfp.(animalse{a}).(names{1})(floor(LFPsall.eyfp.(animalse{a}).pCSts{c})),LFPsall.eyfp.(animalse{a}).(names{2})(floor(LFPsall.eyfp.(animalse{a}).pCSts{c})),2048,1024,4096,2000);
%             
%             [LFPsall.eyfp.(animalse{a}).ToneCoh.BLAbf.CSp_freqtime{c},Freq_c,Time_c]=mscohere(LFPsall.eyfp.(animalse{a}).(names{1})(floor(LFPsall.eyfp.(animalse{a}).CSts{c})),LFPsall.eyfp.(animalse{a}).(names{3})(floor(LFPsall.eyfp.(animalse{a}).CSts{c})),2048,1024,4096,2000);
%             [LFPsall.eyfp.(animalse{a}).ToneCoh.BLAbf.pCSp_freqtime{c},foo1,foo2]=mscohere(LFPsall.eyfp.(animalse{a}).(names{1})(floor(LFPsall.eyfp.(animalse{a}).pCSts{c})),LFPsall.eyfp.(animalse{a}).(names{3})(floor(LFPsall.eyfp.(animalse{a}).pCSts{c})),2048,1024,4096,2000);
%             
            [LFPsall.eyfp.(animalse{a}).ToneCoh.BFpfc.CSp_freqtime{c},Freq_c,Time_c]=mscohere(LFPsall.eyfp.(animalse{a}).(names{3})(floor(LFPsall.eyfp.(animalse{a}).CSts{c})),LFPsall.eyfp.(animalse{a}).(names{2})(floor(LFPsall.eyfp.(animalse{a}).CSts{c})),4096,2000,2048,1024); %4096,2000,2048,1024,2
            [LFPsall.eyfp.(animalse{a}).ToneCoh.BFpfc.pCSp_freqtime{c},foo1,foo2]=mscohere(LFPsall.eyfp.(animalse{a}).(names{3})(floor(LFPsall.eyfp.(animalse{a}).pCSts{c})),LFPsall.eyfp.(animalse{a}).(names{2})(floor(LFPsall.eyfp.(animalse{a}).pCSts{c})),4096,2000,2048,1024);
        end
    end
end

%% 5 ARCH Power and Coherence Calculation
%This section will calculate the power and coherence for each trial and
%pre-trial section for each ARCH animal for each electrode. This
%information is added into our LFP structure within each animals level.

animalsa={'Mrgn5','Mrgn10','Mrgn11'};
%animalsa={'Mrgn4','Mrgn5','Mrgn10','Mrgn11'};

for a=1:length(animalsa)
    for b=1:length(names)
        for c=1:numCS
            [LFPsall.arch.(animalsa{a}).TonePwr.(names{b}).CSp_freqtime{c},Freq,Time]=mtcsg(LFPsall.arch.(animalsa{a}).(names{b})(floor(LFPsall.arch.(animalsa{a}).CSts{c})),4096,2000,2048,1024,2);%2048,2000,500,250,2
            [LFPsall.arch.(animalsa{a}).TonePwr.(names{b}).pCSp_freqtime{c},foo1,foo2]=mtcsg(LFPsall.arch.(animalsa{a}).(names{b})(floor(LFPsall.arch.(animalsa{a}).pCSts{c})),4096,2000,2048,1024,2);
        end
    end
end

animalsa=intersect(animalsa_bf,animalsa_pfc);
for a=1:length(animalsa)
    for b=1:length(names)
        for c=1:numCS
%             [LFPsall.arch.(animalsa{a}).ToneCoh.BLApfc.CSp_freqtime{c},Freq_c,Time_c]=mscohere(LFPsall.arch.(animalsa{a}).(names{1})(floor(LFPsall.arch.(animalsa{a}).CSts{c})),LFPsall.arch.(animalsa{a}).(names{2})(floor(LFPsall.arch.(animalsa{a}).CSts{c})),2048,1024,4096,2000);
%             [LFPsall.arch.(animalsa{a}).ToneCoh.BLApfc.pCSp_freqtime{c},foo1,foo2]=mscohere(LFPsall.arch.(animalsa{a}).(names{1})(floor(LFPsall.arch.(animalsa{a}).pCSts{c})),LFPsall.arch.(animalsa{a}).(names{2})(floor(LFPsall.arch.(animalsa{a}).pCSts{c})),2048,1024,4096,2000);
%             
%             [LFPsall.arch.(animalsa{a}).ToneCoh.BLAbf.CSp_freqtime{c},Freq_c,Time_c]=mscohere(LFPsall.arch.(animalsa{a}).(names{1})(floor(LFPsall.arch.(animalsa{a}).CSts{c})),LFPsall.arch.(animalsa{a}).(names{3})(floor(LFPsall.arch.(animalsa{a}).CSts{c})),2048,1024,4096,2000);
%             [LFPsall.arch.(animalsa{a}).ToneCoh.BLAbf.pCSp_freqtime{c},foo1,foo2]=mscohere(LFPsall.arch.(animalsa{a}).(names{1})(floor(LFPsall.arch.(animalsa{a}).pCSts{c})),LFPsall.arch.(animalsa{a}).(names{3})(floor(LFPsall.arch.(animalsa{a}).pCSts{c})),2048,1024,4096,2000);
            
            [LFPsall.arch.(animalsa{a}).ToneCoh.BFpfc.CSp_freqtime{c},Freq_c,Time_c]=mscohere(LFPsall.arch.(animalsa{a}).(names{3})(floor(LFPsall.arch.(animalsa{a}).CSts{c})),LFPsall.arch.(animalsa{a}).(names{2})(floor(LFPsall.arch.(animalsa{a}).CSts{c})),2048,1024,4096,2000);
            [LFPsall.arch.(animalsa{a}).ToneCoh.BFpfc.pCSp_freqtime{c},foo1,foo2]=mscohere(LFPsall.arch.(animalsa{a}).(names{3})(floor(LFPsall.arch.(animalsa{a}).pCSts{c})),LFPsall.arch.(animalsa{a}).(names{2})(floor(LFPsall.arch.(animalsa{a}).pCSts{c})),2048,1024,4096,2000);
         end
    end
end

%% 6 Restructuring Power for Arch vs Eyfp
%Within this section we will seperate out the power analysis for each
%region, for each experimental group, and for each animal. After removing
%it into a seperate structure, we will norm the CS to the preCS. 

%Setting up strcutures to hold information about CS,pCS pwr for each of the
%animals with verified placements
animalse_bla=[animalse(1),animalse(2),animalse(3),animalse(5)];
animalse_bf=[animalse(1),animalse(2),animalse(3),animalse(4),animalse(5)];
animalse_pfc=[animalse(1),animalse(2),animalse(3),animalse(4)];

animalsa_bla=[animalsa(1),animalsa(2),animalsa(3)];
animalsa_bf=[animalsa(1)];
animalsa_pfc=[animalsa(2)];

%Including Pwr calc only for verified Eyfp animals
for n=1:length(animalse_bla)
    Pwr.BLA.eyfp.(animalse_bla{n})=struct('pwr_CS',[],'pwr_pCS',[],'norm',[]);
    for p=1:numCS
            Pwr.BLA.eyfp.(animalse_bla{n}).pwr_CS(:,p)=mean(LFPsall.eyfp.(animalse_bla{n}).TonePwr.BLA.CSp_freqtime{p},2);
            Pwr.BLA.eyfp.(animalse_bla{n}).pwr_pCS(:,p)=mean(LFPsall.eyfp.(animalse_bla{n}).TonePwr.BLA.pCSp_freqtime{p},2);
            Pwr.BLA.eyfp.(animalse_bla{n}).norm(:,p)=(Pwr.BLA.eyfp.(animalse_bla{n}).pwr_CS(:,p))./(Pwr.BLA.eyfp.(animalse_bla{n}).pwr_pCS(:,p));
    end
    foo=struct2cell(Pwr.BLA.eyfp.(animalse_bla{n}));
    Pwr.BLA.eyfp.(animalse_bla{n}).avg_CS=mean(foo{1,1},2);
    Pwr.BLA.eyfp.(animalse_bla{n}).avg_pCS=mean(foo{2,1},2);
    Pwr.BLA.eyfp.(animalse_bla{n}).avg_normCS=mean(foo{3,1},2);
    Pwr.BLA.eyfp.CS(:,n)=mean(foo{1,1},2);
    Pwr.BLA.eyfp.pCS(:,n)=mean(foo{2,1},2);
    Pwr.BLA.eyfp.norm(:,n)=mean(foo{3,1},2);
    Pwr.BLA.eyfp.CS1_3(:,n)= mean(foo{1,1}(:,1:3),2);
    Pwr.BLA.eyfp.pCS1_3(:,n)= mean(foo{2,1}(:,1:3),2);
    Pwr.BLA.eyfp.norm1_3(:,n)= mean(foo{3,1}(:,1:3),2);
    Pwr.BLA.eyfp.CS8_10(:,n)= mean(foo{1,1}(:,8:10),2);
    Pwr.BLA.eyfp.pCS8_10(:,n)= mean(foo{2,1}(:,8:10),2);
    Pwr.BLA.eyfp.norm8_10(:,n)= mean(foo{3,1}(:,8:10),2);
end

for n=1:length(animalse_bf)
    Pwr.BF.eyfp.(animalse_bf{n})=struct('pwr_CS',[],'pwr_pCS',[],'norm',[]);
    for p=1:numCS
            Pwr.BF.eyfp.(animalse_bf{n}).pwr_CS(:,p)=mean(LFPsall.eyfp.(animalse_bf{n}).TonePwr.BF.CSp_freqtime{p},2);
            Pwr.BF.eyfp.(animalse_bf{n}).pwr_pCS(:,p)=mean(LFPsall.eyfp.(animalse_bf{n}).TonePwr.BF.pCSp_freqtime{p},2);
            Pwr.BF.eyfp.(animalse_bf{n}).norm(:,p)=(Pwr.BF.eyfp.(animalse_bf{n}).pwr_CS(:,p))./(Pwr.BF.eyfp.(animalse_bf{n}).pwr_pCS(:,p));
    end
    foo=struct2cell(Pwr.BF.eyfp.(animalse_bf{n}));
    Pwr.BF.eyfp.(animalse_bf{n}).avg_CS=mean(foo{1,1},2);
    Pwr.BF.eyfp.(animalse_bf{n}).avg_pCS=mean(foo{2,1},2);
    Pwr.BF.eyfp.(animalse_bf{n}).avg_normCS=mean(foo{3,1},2);
    Pwr.BF.eyfp.CS(:,n)=mean(foo{1,1},2);
    Pwr.BF.eyfp.pCS(:,n)=mean(foo{2,1},2);
    Pwr.BF.eyfp.norm(:,n)=mean(foo{3,1},2);
    Pwr.BF.eyfp.CS1_3(:,n)= mean(foo{1,1}(:,1:3),2);
    Pwr.BF.eyfp.pCS1_3(:,n)= mean(foo{2,1}(:,1:3),2);
    Pwr.BF.eyfp.norm1_3(:,n)= mean(foo{3,1}(:,1:3),2);
    Pwr.BF.eyfp.CS8_10(:,n)= mean(foo{1,1}(:,8:10),2);
    Pwr.BF.eyfp.pCS8_10(:,n)= mean(foo{2,1}(:,8:10),2);
    Pwr.BF.eyfp.norm8_10(:,n)= mean(foo{3,1}(:,8:10),2);
end

for n=1:length(animalse_pfc)
    Pwr.PFC.eyfp.(animalse_pfc{n})=struct('pwr_CS',[],'pwr_pCS',[],'norm',[]);
    for p=1:numCS
            Pwr.PFC.eyfp.(animalse_pfc{n}).pwr_CS(:,p)=mean(LFPsall.eyfp.(animalse_pfc{n}).TonePwr.PFC.CSp_freqtime{p},2);
            Pwr.PFC.eyfp.(animalse_pfc{n}).pwr_pCS(:,p)=mean(LFPsall.eyfp.(animalse_pfc{n}).TonePwr.PFC.pCSp_freqtime{p},2);
            Pwr.PFC.eyfp.(animalse_pfc{n}).norm(:,p)=(Pwr.PFC.eyfp.(animalse_pfc{n}).pwr_CS(:,p))./(Pwr.PFC.eyfp.(animalse_pfc{n}).pwr_pCS(:,p));
    end
    foo=struct2cell(Pwr.PFC.eyfp.(animalse_pfc{n}));
    Pwr.PFC.eyfp.(animalse_pfc{n}).avg_CS=mean(foo{1,1},2);
    Pwr.PFC.eyfp.(animalse_pfc{n}).avg_pCS=mean(foo{2,1},2);
    Pwr.PFC.eyfp.(animalse_pfc{n}).avg_normCS=mean(foo{3,1},2);
    Pwr.PFC.eyfp.CS(:,n)=mean(foo{1,1},2);
    Pwr.PFC.eyfp.pCS(:,n)=mean(foo{2,1},2);
    Pwr.PFC.eyfp.norm(:,n)=mean(foo{3,1},2);
    Pwr.PFC.eyfp.CS1_3(:,n)= mean(foo{1,1}(:,1:3),2);
    Pwr.PFC.eyfp.pCS1_3(:,n)= mean(foo{2,1}(:,1:3),2);
    Pwr.PFC.eyfp.norm1_3(:,n)= mean(foo{3,1}(:,1:3),2);
    Pwr.PFC.eyfp.CS8_10(:,n)= mean(foo{1,1}(:,8:10),2);
    Pwr.PFC.eyfp.pCS8_10(:,n)= mean(foo{2,1}(:,8:10),2);
    Pwr.PFC.eyfp.norm8_10(:,n)= mean(foo{3,1}(:,8:10),2);
end


%Including Pwr calc only for verified Ach animals
for n=1:length(animalsa_bla)
    Pwr.BLA.arch.(animalsa_bla{n})=struct('pwr_CS',[],'pwr_pCS',[],'norm',[]);
    for p=1:numCS
            Pwr.BLA.arch.(animalsa_bla{n}).pwr_CS(:,p)=mean(LFPsall.arch.(animalsa_bla{n}).TonePwr.BLA.CSp_freqtime{p},2);
            Pwr.BLA.arch.(animalsa_bla{n}).pwr_pCS(:,p)=mean(LFPsall.arch.(animalsa_bla{n}).TonePwr.BLA.pCSp_freqtime{p},2);
            Pwr.BLA.arch.(animalsa_bla{n}).norm(:,p)=(Pwr.BLA.arch.(animalsa_bla{n}).pwr_CS(:,p))./(Pwr.BLA.arch.(animalsa_bla{n}).pwr_pCS(:,p));
    end
    foo=struct2cell(Pwr.BLA.arch.(animalsa_bla{n}));
    Pwr.BLA.arch.(animalsa_bla{n}).avg_CS=mean(foo{1,1},2);
    Pwr.BLA.arch.(animalsa_bla{n}).avg_pCS=mean(foo{2,1},2);
    Pwr.BLA.arch.(animalsa_bla{n}).avg_normCS=mean(foo{3,1},2);
    Pwr.BLA.arch.CS(:,n)=mean(foo{1,1},2);
    Pwr.BLA.arch.pCS(:,n)=mean(foo{2,1},2);
    Pwr.BLA.arch.norm(:,n)=mean(foo{3,1},2);
    Pwr.BLA.arch.CS1_3(:,n)= mean(foo{1,1}(:,1:3),2);
    Pwr.BLA.arch.pCS1_3(:,n)= mean(foo{2,1}(:,1:3),2);
    Pwr.BLA.arch.norm1_3(:,n)= mean(foo{3,1}(:,1:3),2);
    Pwr.BLA.arch.CS8_10(:,n)= mean(foo{1,1}(:,8:10),2);
    Pwr.BLA.arch.pCS8_10(:,n)= mean(foo{2,1}(:,8:10),2);
    Pwr.BLA.arch.norm8_10(:,n)= mean(foo{3,1}(:,8:10),2);
end

for n=1:length(animalsa_bf)
    Pwr.BF.arch.(animalsa_bf{n})=struct('pwr_CS',[],'pwr_pCS',[],'norm',[]);
    for p=1:numCS
            Pwr.BF.arch.(animalsa_bf{n}).pwr_CS(:,p)=mean(LFPsall.arch.(animalsa_bf{n}).TonePwr.BF.CSp_freqtime{p},2);
            Pwr.BF.arch.(animalsa_bf{n}).pwr_pCS(:,p)=mean(LFPsall.arch.(animalsa_bf{n}).TonePwr.BF.pCSp_freqtime{p},2);
            Pwr.BF.arch.(animalsa_bf{n}).norm(:,p)=(Pwr.BF.arch.(animalsa_bf{n}).pwr_CS(:,p))./(Pwr.BF.arch.(animalsa_bf{n}).pwr_pCS(:,p));
    end
    foo=struct2cell(Pwr.BF.arch.(animalsa_bf{n}));
    Pwr.BF.arch.(animalsa_bf{n}).avg_CS=mean(foo{1,1},2);
    Pwr.BF.arch.(animalsa_bf{n}).avg_pCS=mean(foo{2,1},2);
    Pwr.BF.arch.(animalsa_bf{n}).avg_normCS=mean(foo{3,1},2);
    Pwr.BF.arch.CS(:,n)=mean(foo{1,1},2);
    Pwr.BF.arch.pCS(:,n)=mean(foo{2,1},2);
    Pwr.BF.arch.norm(:,n)=mean(foo{3,1},2);
    Pwr.BF.arch.CS1_3(:,n)= mean(foo{1,1}(:,1:3),2);
    Pwr.BF.arch.pCS1_3(:,n)= mean(foo{2,1}(:,1:3),2);
    Pwr.BF.arch.norm1_3(:,n)= mean(foo{3,1}(:,1:3),2);
    Pwr.BF.arch.CS8_10(:,n)= mean(foo{1,1}(:,8:10),2);
    Pwr.BF.arch.pCS8_10(:,n)= mean(foo{2,1}(:,8:10),2);
    Pwr.BF.arch.norm8_10(:,n)= mean(foo{3,1}(:,8:10),2);
end

for n=1:length(animalsa_pfc)
    Pwr.PFC.arch.(animalsa_pfc{n})=struct('pwr_CS',[],'pwr_pCS',[],'norm',[]);
    for p=1:numCS
            Pwr.PFC.arch.(animalsa_pfc{n}).pwr_CS(:,p)=mean(LFPsall.arch.(animalsa_pfc{n}).TonePwr.PFC.CSp_freqtime{p},2);
            Pwr.PFC.arch.(animalsa_pfc{n}).pwr_pCS(:,p)=mean(LFPsall.arch.(animalsa_pfc{n}).TonePwr.PFC.pCSp_freqtime{p},2);
            Pwr.PFC.arch.(animalsa_pfc{n}).norm(:,p)=(Pwr.PFC.arch.(animalsa_pfc{n}).pwr_CS(:,p))./(Pwr.PFC.arch.(animalsa_pfc{n}).pwr_pCS(:,p));
    end
    foo=struct2cell(Pwr.PFC.arch.(animalsa_pfc{n}));
    Pwr.PFC.arch.(animalsa_pfc{n}).avg_CS=mean(foo{1,1},2);
    Pwr.PFC.arch.(animalsa_pfc{n}).avg_pCS=mean(foo{2,1},2);
    Pwr.PFC.arch.(animalsa_pfc{n}).avg_normCS=mean(foo{3,1},2);
    Pwr.PFC.arch.(animalsa_pfc{n}).avg_CS1_3=mean(foo{1,1}(:,1:3),2);
    Pwr.PFC.arch.(animalsa_pfc{n}).avg_pCS1_3=mean(foo{2,1}(:,1:3),2);
    Pwr.PFC.arch.(animalsa_pfc{n}).avg_norm1_3=mean(foo{3,1}(:,1:3),2);
    Pwr.PFC.arch.CS(:,n)=mean(foo{1,1},2);
    Pwr.PFC.arch.pCS(:,n)=mean(foo{2,1},2);
    Pwr.PFC.arch.norm(:,n)=mean(foo{3,1},2);
    Pwr.PFC.arch.CS1_3(:,n)= mean(foo{1,1}(:,1:3),2);
    Pwr.PFC.arch.pCS1_3(:,n)= mean(foo{2,1}(:,1:3),2);
    Pwr.PFC.arch.norm1_3(:,n)= mean(foo{3,1}(:,1:3),2);
    Pwr.PFC.arch.CS8_10(:,n)= mean(foo{1,1}(:,8:10),2);
    Pwr.PFC.arch.pCS8_10(:,n)= mean(foo{2,1}(:,8:10),2);
    Pwr.PFC.arch.norm8_10(:,n)= mean(foo{3,1}(:,8:10),2);
end


%% Plotting Power
% BF and BLA power preCS/CS
figure('color',[1 1 1],'name','Absolute Power'); 
subplot(2,1,1)
errorbarplot_joe(Freq,Pwr.BLA.arch.CS',[0 1 0],[0.2 0.8 0.2])
hold on
errorbarplot_joe(Freq,Pwr.BLA.eyfp.CS',[0 0 0],[0.5 0.5 0.5])
xlim([0 55])
title('BLA Arch CS, eYFP CS')

subplot(2,1,2)
errorbarplot_joe(Freq,Pwr.BF.arch.CS',[0 1 0],[0.2 0.8 0.2])
hold on
errorbarplot_joe(Freq,Pwr.BF.eyfp.CS',[0 0 0],[0.5 0.5 0.5])
xlim([0 55])
title('BF Arch CS, eYFP CS')

% BF/BLA power Normalized 
figure('color',[1 1 1],'name','Normalized Power'); 
subplot(3,1,1)
errorbarplot_joe(Freq,Pwr.BLA.arch.norm1_3',[1 0 0],[0.9 0.2 0.2])
hold on
errorbarplot_joe(Freq,Pwr.BLA.arch.norm8_10',[0 0 0],[0.5 0.5 0.5])
xlim([0 55])
title('BLA trials 1-3, trials 8_10, n=3')

subplot(3,1,2)
errorbarplot_joe(Freq,Pwr.BF.eyfp.norm1_3',[1 0 0],[0.8 0.2 0.2])
hold on
errorbarplot_joe(Freq,Pwr.BF.eyfp.norm8_10',[0 0 0],[0.5 0.5 0.5])
xlim([0 55])
title('BF trials 1-3, trials 8_10, n=1')

subplot(3,1,3)
errorbarplot_joe(Freq,Pwr.PFC.arch.norm1_3',[1 0 0],[0.8 0.2 0.2])
hold on
errorbarplot_joe(Freq,Pwr.PFC.arch.norm8_10',[0 0 0],[0.5 0.5 0.5])
xlim([0 55])
title('PFC trials 1-3, trials 8_10, n=1')

figure('color',[1 1 1],'name','BF Power'); 
errorbarplot_joe(Freq,Pwr.BF.arch.norm1_3',[1 0 0],[0.8 0.2 0.2])
hold on
errorbarplot_joe(Freq,Pwr.BF.arch.norm8_10',[0 0 1],[0.2 0.2 0.8])
xlim([0 50])
title('BF arch Norm Power Tr 1-3, 8-10')

figure('color',[1 1 1],'name','BF arch M2 (1:3) eyfp M2 Trials 8:10'); 
errorbarplot_joe(Freq,[Pwr.BF.eyfp.Mrgn2.pwr_CS(:,1:3)./1000]',[1 0 0],[0.8 0.2 0.2])
hold on
errorbarplot_joe(Freq,[Pwr.BF.eyfp.Mrgn2.pwr_CS(:,8:10)./1000]',[0 0 1],[0.2 0.2 0.8])
xlim([10 50])
ylim([0 10^3])
ylabel('Power (mV/Hz^2)')
xlabel('Frequency (Hz)')

figure('color',[1 1 1],'name','BF arch M5 (1:3) arch M5 Trials 8:10'); 
errorbarplot_joe(Freq,[Pwr.BF.arch.Mrgn5.pwr_CS(:,1:3)./10000]',[1 0 0],[0.8 0.2 0.2])
hold on
%errorbarplot_joe(Freq,[Pwr.BF.arch.Mrgn5.pwr_CS(:,8:10)./10000]',[0 0 1],[0.2 0.2 0.8])
errorbarplot_joe(Freq,[Pwr.BF.arch.Mrgn5.pwr_CS(:,8:10)./10000]',[0 0 1],[0.2 0.2 0.8])
xlim([10 50])
ylim([0 10^3])
ylabel('Power (mV/Hz^2)')
xlabel('Frequency (Hz)')

delta_freq_idx=find(Freq>=0.1 & Freq<4);
theta_freq_idx=find(Freq>=4 & Freq<12);
gamma_freq_idx=find(Freq>=30 & Freq<=50);
hgamma_freq_idx=find(Freq>=70 & Freq<=120);

arch_CS_SEM = std(Pwr.BLA.arch.norm(theta_freq_idx))/sqrt(size(Pwr.BLA.arch.norm,2));
eyfp_CS_SEM = std(Pwr.BLA.eyfp.norm(theta_freq_idx))/sqrt(size(Pwr.BLA.eyfp.norm,2));

arch_CS_thetadelta_SEM = std(Pwr.BLA.arch.norm(theta_freq_idx)./Pwr.BLA.arch.norm(delta_freq_idx))/sqrt(size(Pwr.BLA.arch.norm,2));
eyfp_CS_thetadelta_SEM = std(Pwr.BLA.eyfp.norm(theta_freq_idx)./Pwr.BLA.eyfp.norm(delta_freq_idx))/sqrt(size(Pwr.BLA.arch.norm,2));


figure5=figure('color','w');
axes1 = axes('Parent',figure5);
hold(axes1,'on');
bar(1,mean(Pwr.BLA.arch.norm(theta_freq_idx)),'g')
hold on
errorbar(1,mean(Pwr.BLA.arch.norm(theta_freq_idx)),arch_CS_SEM,'g.')
hold on
bar(2,mean(Pwr.BLA.eyfp.norm(theta_freq_idx)),'w')
hold on
errorbar(2,mean(Pwr.BLA.eyfp.norm(theta_freq_idx)),eyfp_CS_SEM,'k.')
xlim([0 3]); %ylim([0 8E6])
set(axes1,'XTick',[1 2],'XTickLabel',{'Arch','eYFP'});
title('BLA norm theta (4-12 Hz) All trials')

figure6=figure('color','w');
axes1 = axes('Parent',figure6);
hold(axes1,'on');
bar(1,mean(Pwr.BLA.arch.norm1_3(theta_freq_idx)./mean(Pwr.BLA.arch.norm1_3(delta_freq_idx))),'g')
hold on
errorbar(1,mean(Pwr.BLA.arch.norm1_3(theta_freq_idx))./mean(Pwr.BLA.arch.norm1_3(delta_freq_idx)),arch_CS_thetadelta_SEM,'g.')
hold on
bar(2,mean(Pwr.BLA.eyfp.norm1_3(theta_freq_idx)./mean(Pwr.BLA.eyfp.norm1_3(delta_freq_idx))),'w')
hold on
errorbar(2,mean(Pwr.BLA.eyfp.norm1_3(theta_freq_idx)./mean(Pwr.BLA.arch.norm1_3(delta_freq_idx),arch_CS_thetadelta_SEM,'k.')))
xlim([0 3]); %ylim([0 8E6])
set(axes1,'XTick',[1 2],'XTickLabel',{'Arch','eYFP'});
title('BLA theta:delta (4-12 Hz) trials 1-3')

figure('color',[1 1 1],'name','PFC Power eYFP arch'); 
errorbarplot_joe(Freq,Pwr.PFC.eyfp.CS1_3',[0 0 0],[0.5 0.5 0.5])
hold on
errorbarplot_joe(Freq,Pwr.PFC.arch.CS1_3',[0 1 0],[0.5 0.8 0.5])
xlim([0 50])
title('PFC yfp arch tr1-3')


figure('color',[1 1 1],'name','BLA Power eYFP & eArch'); 
errorbarplot_joe(Freq,Pwr.BLA.eyfp.norm',[0 0 0],[0.5 0.5 0.5])
hold on
errorbarplot_joe(Freq,Pwr.BLA.arch.norm',[0 1 0],[0.2 0.8 0.2])
xlim([0 50])

figure('color',[1 1 1],'name','eYFP vs Arch Power All Trials');
subplot(3,1,1)
errorbarplot_joe(Freq,Pwr.BLA.eyfp.norm',[0 0 0],[0.5 0.5 0.5])
hold on
errorbarplot_joe(Freq,Pwr.BLA.arch.norm',[0 1 0],[0.2 0.8 0.2])
xlim([0 80])
legend('Eyfp','','Arch','')
title('BLA')

subplot(3,1,2)
errorbarplot_joe(Freq,Pwr.BF.eyfp.norm',[0 0 0],[0.5 0.5 0.5])
hold on
errorbarplot_joe(Freq,Pwr.BF.arch.norm',[0 1 0],[0.2 0.8 0.2])
xlim([0 80])
legend('Eyfp','','Arch','')
title('BF')

subplot(3,1,3)
errorbarplot_joe(Freq,Pwr.PFC.eyfp.norm',[0 0 0],[0.5 0.5 0.5])
hold on
errorbarplot_joe(Freq,Pwr.PFC.arch.norm',[0 1 0],[0.2 0.8 0.2])
xlim([0 80])
legend('Eyfp','','Arch','')
title('PFC')

%% Plotting power
figure('color',[1 1 1],'name','BF')
plot(Freq, mean(Pwr.BF.arch.Mrgn5.pwr_CS(:,1:3),2),'-r')
xlim([0 55])
hold on
plot(Freq, mean(Pwr.BF.arch.Mrgn5.pwr_CS(:,8:10),2),'-k')
title('mrgn 11')
legend ('BF Trials 1-3','Trials 8-10'), box off

%% 7 Restructuring Coherence for EYFP vs ARCH
%This section reorganizes the coherence so that it is in a seperate
%structure for manipulation or visualization.

animalse=intersect(animalse_bf,animalse_pfc);

for n=1:length(animalse)  
%     Coh.BLApfc.eyfp.CSp(n).CSp=cell2mat(LFPsall.eyfp.(animalse{n}).ToneCoh.BLApfc.CSp_freqtime);
%     Coh.BLApfc.eyfp.pCSp(n).preCSp=cell2mat(LFPsall.eyfp.(animalse{n}).ToneCoh.BLApfc.pCSp_freqtime);
%     
%     Coh.BLAbf.eyfp.CSp(n).CSp=cell2mat(LFPsall.eyfp.(animalse{n}).ToneCoh.BLAbf.CSp_freqtime);
%     Coh.BLAbf.eyfp.pCSp(n).preCSp=cell2mat(LFPsall.eyfp.(animalse{n}).ToneCoh.BLAbf.pCSp_freqtime);
%     
    Coh.BFpfc.eyfp.CSp(n).CSp=cell2mat(LFPsall.eyfp.(animalse{n}).ToneCoh.BFpfc.CSp_freqtime);
    Coh.BFpfc.eyfp.pCSp(n).preCSp=cell2mat(LFPsall.eyfp.(animalse{n}).ToneCoh.BFpfc.pCSp_freqtime);
end

animalsa=intersect(animalsa_bla,animalsa_bf);

for n=1:length(animalsa)
%     Coh.BLApfc.arch.CSp(n).CSp=cell2mat(LFPsall.arch.(animalsa{n}).ToneCoh.BLApfc.CSp_freqtime);
%     Coh.BLApfc.arch.pCSp(n).preCSp=cell2mat(LFPsall.arch.(animalsa{n}).ToneCoh.BLApfc.pCSp_freqtime);
%     
    Coh.BLAbf.arch.CSp(n).CSp=cell2mat(LFPsall.arch.(animalsa{n}).ToneCoh.BLAbf.CSp_freqtime);
    Coh.BLAbf.arch.pCSp(n).preCSp=cell2mat(LFPsall.arch.(animalsa{n}).ToneCoh.BLAbf.pCSp_freqtime);
%     
%     Coh.BFpfc.arch.CSp(n).CSp=cell2mat(LFPsall.arch.(animalsa{n}).ToneCoh.BFpfc.CSp_freqtime);
%     Coh.BFpfc.arch.pCSp(n).preCSp=cell2mat(LFPsall.arch.(animalsa{n}).ToneCoh.BFpfc.pCSp_freqtime);
end

%% 8 Calculating Correlation Between Regions
animalsa=intersect(animalsa_bf,animalsa_pfc);

for a=1:length(animalsa)
    for c=1:numCS
        
%         [lags, Corr.BLApfc.arch.(animalsa{a}).corr_CS(:,c),Corr.BLApfc.arch.maxLagCS(a,c)]=amp_crosscorr(LFPsall.arch.(animalsa{a}).(names{1})(floor(LFPsall.arch.(animalsa{a}).CSts{c})),LFPsall.arch.(animalsa{a}).(names{3})(floor(LFPsall.arch.(animalsa{a}).CSts{c})),2000,20,50);
%         [lags, Corr.BLApfc.arch.(animalsa{a}).corr_pCS(:,c),Corr.BLApfc.arch.maxLagpCS(a,c)]=amp_crosscorr(LFPsall.arch.(animalsa{a}).(names{1})(floor(LFPsall.arch.(animalsa{a}).pCSts{c})),LFPsall.arch.(animalsa{a}).(names{3})(floor(LFPsall.arch.(animalsa{a}).pCSts{c})),2000,20,50);
%          
%         [lags, Corr.BLAbf.arch.(animalsa{a}).corr_CS(:,c),Corr.BLAbf.arch.maxLagCS(a,c)]=amp_crosscorr(LFPsall.arch.(animalsa{a}).(names{1})(floor(LFPsall.arch.(animalsa{a}).CSts{c})),LFPsall.arch.(animalsa{a}).(names{2})(floor(LFPsall.arch.(animalsa{a}).CSts{c})),2000,20,50);
%         [lags, Corr.BLAbf.arch.(animalsa{a}).corr_pCS(:,c),Corr.BLAbf.arch.maxLagpCS(a,c)]=amp_crosscorr(LFPsall.arch.(animalsa{a}).(names{1})(floor(LFPsall.arch.(animalsa{a}).pCSts{c})),LFPsall.arch.(animalsa{a}).(names{2})(floor(LFPsall.arch.(animalsa{a}).pCSts{c})),2000,20,50);
%          
        [lags, Corr.BFpfc.arch.(animalsa{a}).corr_CS(:,c),Corr.BFpfc.arch.maxLagCS(a,c)]=amp_crosscorr(LFPsall.arch.(animalsa{a}).(names{2})(floor(LFPsall.arch.(animalsa{a}).CSts{c})),LFPsall.arch.(animalsa{a}).(names{3})(floor(LFPsall.arch.(animalsa{a}).CSts{c})),2000,20,50);
        [lags, Corr.BFpfc.arch.(animalsa{a}).corr_pCS(:,c),Corr.BFpfc.arch.maxLagpCS(a,c)]=amp_crosscorr(LFPsall.arch.(animalsa{a}).(names{2})(floor(LFPsall.arch.(animalsa{a}).pCSts{c})),LFPsall.arch.(animalsa{a}).(names{3})(floor(LFPsall.arch.(animalsa{a}).pCSts{c})),2000,4,12); 
     close all
    end
end

animalse=intersect(animalse_bf,animalse_pfc);

for a=1:length(animalse)
    for c=1:numCS
%         [lags, Corr.BLApfc.eyfp.(animalse{a}).corr_CS(:,c),Corr.BLApfc.eyfp.maxLagCS(a,c)]=amp_crosscorr(LFPsall.eyfp.(animalse{a}).(names{1})(floor(LFPsall.eyfp.(animalse{a}).CSts{c})),LFPsall.eyfp.(animalse{a}).(names{3})(floor(LFPsall.eyfp.(animalse{a}).CSts{c})),2000,20,50);
%         [lags, Corr.BLApfc.eyfp.(animalse{a}).corr_pCS(:,c),Corr.BLApfc.eyfp.maxLagpCS(a,c)]=amp_crosscorr(LFPsall.eyfp.(animalse{a}).(names{1})(floor(LFPsall.eyfp.(animalse{a}).pCSts{c})),LFPsall.eyfp.(animalse{a}).(names{3})(floor(LFPsall.eyfp.(animalse{a}).pCSts{c})),2000,20,50);
        
%         [lags, Corr.BLAbf.eyfp.(animalse{a}).corr_CS(:,c),Corr.BLAbf.eyfp.maxLagCS(a,c)]=amp_crosscorr(LFPsall.eyfp.(animalse{a}).(names{1})(floor(LFPsall.eyfp.(animalse{a}).CSts{c})),LFPsall.eyfp.(animalse{a}).(names{2})(floor(LFPsall.eyfp.(animalse{a}).CSts{c})),2000,20,50);
%         [lags, Corr.BLAbf.eyfp.(animalse{a}).corr_pCS(:,c),Corr.BLAbf.eyfp.maxLagpCS(a,c)]=amp_crosscorr(LFPsall.eyfp.(animalse{a}).(names{1})(floor(LFPsall.eyfp.(animalse{a}).pCSts{c})),LFPsall.eyfp.(animalse{a}).(names{2})(floor(LFPsall.eyfp.(animalse{a}).pCSts{c})),2000,20,50);
        
        [lags, Corr.BFpfc.eyfp.(animalse{a}).corr_CS(:,c),Corr.BFpfc.eyfp.maxLagCS(a,c)]=amp_crosscorr(LFPsall.eyfp.(animalse{a}).(names{2})(floor(LFPsall.eyfp.(animalse{a}).CSts{c})),LFPsall.eyfp.(animalse{a}).(names{3})(floor(LFPsall.eyfp.(animalse{a}).CSts{c})),2000,20,50);
        [lags, Corr.BFpfc.eyfp.(animalse{a}).corr_pCS(:,c),Corr.BFpfc.eyfp.maxLagpCS(a,c)]=amp_crosscorr(LFPsall.eyfp.(animalse{a}).(names{2})(floor(LFPsall.eyfp.(animalse{a}).pCSts{c})),LFPsall.eyfp.(animalse{a}).(names{3})(floor(LFPsall.eyfp.(animalse{a}).pCSts{c})),2000,20,50); 
    close all
    end
end           

Corr.Lags = lags;

%% Graphing Power Cross-Correlations

%BLA-pfc CS and pCS
MaxLags.eyfp.CS=struct('BLApfc',reshape(Corr.BLApfc.eyfp.maxLagCS,[],1),'BFpfc',reshape(Corr.BFpfc.eyfp.maxLagCS,[],1),'BLAbf',reshape(Corr.BLAbf.eyfp.maxLagCS,[],1));
MaxLags.eyfp.pCS=struct('BLApfc',reshape(Corr.BLApfc.eyfp.maxLagpCS,[],1),'BFpfc',reshape(Corr.BFpfc.eyfp.maxLagpCS,[],1),'BLAbf',reshape(Corr.BLAbf.eyfp.maxLagpCS,[],1));

MaxLags.arch.CS=struct('BLApfc',reshape(Corr.BLApfc.arch.maxLagCS,[],1),'BFpfc',reshape(Corr.BFpfc.arch.maxLagCS,[],1),'BLAbf',reshape(Corr.BLAbf.arch.maxLagCS,[],1));
MaxLags.arch.pCS=struct('BLApfc',reshape(Corr.BLApfc.arch.maxLagpCS,[],1),'BFpfc',reshape(Corr.BFpfc.arch.maxLagpCS,[],1),'BLAbf',reshape(Corr.BLAbf.arch.maxLagpCS,[],1));

MaxLags.arch.CS=struct('BLApfc',reshape(Corr.BLApfc.arch.maxLagCS,[],1),'BLAbf',reshape(Corr.BLAbf.arch.maxLagCS,[],1));
MaxLags.arch.pCS=struct('BLApfc',reshape(Corr.BLApfc.arch.maxLagpCS,[],1),'BLAbf',reshape(Corr.BLAbf.arch.maxLagpCS,[],1));

figplots=[{'pCS'},{'CS'}];

figure('color',[1 1 1],'name','eyfp BF-PFC cross Corrs 20-50Hz')
subplot(1,2,1)
histogram(MaxLags.eyfp.pCS.BFpfc,30)
hold on; plot([0 0],[0 9],'k:','linewidth',2)
xlim([-80 80])
title(figplots(1))
    
subplot(1,2,2)
histogram(MaxLags.eyfp.CS.BFpfc,30)
hold on; plot([0 0],[0 10],'k:','linewidth',2)
xlim([-80 80])
title(figplots(2))
    
signrank(MaxLags.eyfp.CS.BFpfc)
mean(MaxLags.eyfp.CS.BFpfc)

%4-8hz BLA-PFC Arch pre-CS is diff from 0, median 10ms mPFC lead, 
%4-8hz BLA-PFC Arch CS is not diff from 0

%4-8hz BLA-PFC eYFP pCS is diff from 0, median 6.75ms, 
%4-12 hz pCS is diff from 0 (p=.000000841), median -14.25ms, mean -15.96ms

%4-8hz BLA-PFC eYPF CS is diff from 0, median 7.5ms
%4-12 CS is diff from 0 p=.0073, mean 8.6ms, median 7.5ms

%4-8hz BF-PFC Arch pCS is diff from 0, median 7ms mPFC lead, mean is 10.68ms
%4-8hz BF-PFC Arch CS is not diff from 0

%4-8hz BF-PFC eYFP pCS is diff from 0, median 3.25ms, 
%4-12Hz pCS diff from 0, median -1.92ms

%4-8hz BF-PFC eYFP CS is diff from 0, median 3.5, 
%4-12 Hz CS is diff from 0 p=.00082, median 3.5ms

%Is PFC starting to lead BF and BLA in 4-8Hz theta during CS in eYFP but
%not in eArch mice? I may need to remove 2 mice with weird data from the
%eyfp group (mouse #4 and 5 in the matrices, which are Mrgn8 and Mrgn9).
test=horzcat(MaxLags.arch.pCS.BLApfc,MaxLags.arch.CS.BLApfc);
p=kruskalwallis(test)

%% Under Construction... Ignore for now.

for c=1 %#ok<COLND>
BF1(:,c) = (LFPsall.eyfp.(animals1{1}).BF(floor(LFPsall.eyfp.(animals1{1}).CSts{c})));
BF2(:,c) = (LFPsall.eyfp.(animals1{2}).BF(floor(LFPsall.eyfp.(animals1{2}).CSts{c})));
BF3(:,c) = (LFPsall.eyfp.(animals1{3}).BF(floor(LFPsall.eyfp.(animals1{3}).CSts{c})));
BF4(:,c) = LFPsall.arch.(animals2{1}).BF(floor(LFPsall.arch.(animals2{1}).CSts{c}));
BF5(:,c) = LFPsall.arch.(animals2{2}).BF(floor(LFPsall.arch.(animals2{2}).CSts{c}));
BF8(:,c) = (LFPsall.eyfp.(animals1{4}).BF(floor(LFPsall.eyfp.(animals1{4}).CSts{c})));
BF9(:,c) = (LFPsall.eyfp.(animals1{5}).BF(floor(LFPsall.eyfp.(animals1{4}).CSts{c})));
BF10(:,c) = LFPsall.arch.(animals2{3}).BF(floor(LFPsall.arch.(animals2{2}).CSts{c}));
BF11(:,c) = LFPsall.arch.(animals2{4}).BF(floor(LFPsall.arch.(animals2{2}).CSts{c}));
BF12(:,c) = (LFPsall.eyfp.(animals1{6}).BF(floor(LFPsall.eyfp.(animals1{4}).CSts{c})));
end

    
%% Under Construction... Ignore for now.
k=[];
j=[];
s=[];
t=[];
%cohpair=1; %1=blapfc, 2=bfpfc, 3=blabf



figure('Color',[1 1 1],'Name','BLAbf')
for j=1:size(Coh_BLAbf(j).CSp,2)
    subplot(1,size(Coh_BLAbf(j).CSp,2),j)
    plot(Freq_c,Coh_BLAbf(j).CSp(:,1),'r')
    hold on
    plot(Freq_c,Coh_BLAbf(j).CSm(:,1),'b')
    xlim([0 15])
    ylim([0 1])
end
