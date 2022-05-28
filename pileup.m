%pulse pile up
clc; clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           Reading spectrum                 %%
%%---------------------------------------------%%
 global SPwoBg Pt
%opening spectrum file
file=fopen('Cs137.txt','r');
%opening backgroud file
file2=fopen('background.txt','r');

formatSpec='%f';
%reading spectrum data
Spectrum=fscanf(file,formatSpec);
%reading background data
BG=fscanf(file2,formatSpec);
%excluding background from spectrum
SPwoBg=Spectrum-BG;
%--------------------------------------%
% load('SpectrumCs137.mat');


%sum will be required for probabilities when sampling from spectrum
Pt=sum(SPwoBg);

%---------------------------------%
%             Input               %
%---------------------------------%
TCR=[8000 40000 80000 400000 600000];  %true count rate
%TCR=600000;     %true count rate
tau=1E-6;       %resolving time
n_pu=2;         %number of pileup pulses to simulate
%-----------------------------------------------------%
for cr=1:length(TCR)
%probability of pile up for True count rate 
%Using Paralyzable model
P_pu=zeros(1,n_pu);
for j=1:n_pu
P_pu(j)=exp(-TCR(cr)*tau)*(1-exp(-TCR(cr)*tau))^(j-1);
end
%normalize probability
P_pu=P_pu./(sum(P_pu));
%-----------------------------------------------------%

%Number of samples for simulation
Samples=2E4;
%vector to save count for both pileup and without pileup counts
channel_counts=zeros(2,length(SPwoBg));
%----------------------------------------------------%
%loop over each sample
for i=1:Samples
    %determine the number of piled up pulses
    %impiles wether the pulse will be piled up/at what degree
    D_p=degree_pileup(P_pu,n_pu);
    %Sampling random pulse from orignal spectrum
     [~,CH1]=Pulse_height();
    %vector to save channel number in case of pile up pulses
     CH=zeros(1,D_p);
        %loop over each pilled up degree
        for t=1:D_p
            %Sampling random pulses from orignal spectrum
            [~,CH(t)]=Pulse_height();
        end
        %saving pilled up count in the corresponding channel
        channel_counts(1,sum(CH))=channel_counts(1,sum(CH))+1;
        %saving non pilled or free count in channel
        channel_counts(2,CH1)=channel_counts(2,CH1)+1;
end

%--------------------------------------------------------------%
%                       Plots                                  %
%--------------------------------------------------------------%
figure(cr)
%plot of piled up spectrum
plot(1:length(SPwoBg),channel_counts(1,:))
hold on
%plot of spectrum with out pile up
plot(1:length(SPwoBg),channel_counts(2,:))
hold off
xlabel('Channel number'); ylabel('dN/dH');
title(['Differential Pulse Height Spectrum at Count rate ' num2str(TCR(cr))]);
legend('Pile up','without Pile up')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            Functions                 %%%
%%----------------------------------------%%
%Number of pile up pulses/degree of pile up
function Np=degree_pileup(P,np)
    r=rand;  %random number btw 0,1
    for i=1:np
        x_l=sum(P(1:i-1)); 
        x_u=sum(P(1:i));  
        %if random number lies with in range of proabability
        % of free pulse or higher order pile up
        if (r>x_l)&&(r<x_u)
            Np=i;
            break
        end
    end
end

%Sampling from spectrum
function [x,i]=Pulse_height()
    global SPwoBg Pt
    K=length(SPwoBg);
    r=rand; %random number
    %pulse heights are scaled btw 0,1 
    %for which random number lies within range, that pulse/Channel is
    %sampled
    for i=1:K
        P_min=sum(SPwoBg(1:i-1))/Pt;
        P_max=sum(SPwoBg(1:i))/Pt;
        if (r>P_min)&&(r<P_max)
            x=SPwoBg(i);
            break
        end
    end
end
