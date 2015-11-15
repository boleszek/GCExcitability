% Code for generating Figure 4

% This script simulates the full model with 45 MCs and 720 dGCs over the
% full range of excitability and across a range of inhibitory weights.
% Parameters are in OB_params_GCE_Fig4.txt

% Boleslaw Osinski (2015)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
input_file = 'OB_params_GCE_Fig4.txt';
[dt,tsim,numtp,nmit,ngradist,ngraprox,sampf,timevec] ...
    = InitNetwork_GCE(input_file);

%% Part Ai,ii
% runtime: ~1.5hr

% P1 controls inhibitory weight from GCs onto MCs, wGABAGR
P1.line = 44;
P1.name = 'wGABAGR ';
P1.val = 0.013:0.001:0.017;

% P2 is Vrest, the GC excitability
P2.line = 56;
P2.name = 'Vrest ';
P2.val = 1e-3 * [-75:1:-55];

% Simulate the model over chosen parameter range
[MCMAT dGCMAT ICMAT param MClfpMAT dGClfpMAT] = ...
    ParamSweep_GCE(P1,P2,numtp,input_file);

%%% Calculate LFP frequency and power %%%

FmaxMAT = zeros(length(P1.val),length(P2.val));
maxpwrMAT = zeros(length(P1.val),length(P2.val));

trim = 500; % Ignore first 500 timepoints

L = length(timevec(trim:end-100));  % Length of data that will be analyzed
NFFT = 2^nextpow2(L); % Next power of 2 from L
f = sampf/2*linspace(0,1,NFFT/2+1);
ROI = ceil(8/(f(2)-f(1))):ceil(140/(f(2)-f(1)));
% ROI gives us a window from ~7Hz - 115Hz from which to select peak power.
% This avoids selecting low freuqency components that may dominate the
% spectrum.

% Run FFT on each simulated LFP
for n1 = 1:length(P1.val)
    for n2 = 1:length(P2.val)
        
        MitLFP = MClfpMAT{n1,n2}.GradistMitGlobal;
        mitFFT = fft(detrend(MitLFP(trim:end-100,1),'constant'),NFFT)/L;

        absmitFFT = 2*abs(mitFFT(1:NFFT/2+1));
        maxpwrMAT(n1,n2) = max(absmitFFT(ROI));

        maxind = find(absmitFFT == maxpwrMAT(n1,n2));
        FmaxMAT(n1,n2) = f(maxind);
    end
end


%%% Plot Ai,ii %%%

Vrest = 1e-3 * [-75:1:-55];
Weights = 0.013:0.001:0.017;
winds = 1:length(Weights);

scrsz = get(0,'ScreenSize');
figH=figure;
set(figH,'position',[0,400,scrsz(3)-0.2*scrsz(3),scrsz(4)-0.5*scrsz(4)]);
subplot(1,2,1)
hold on
for ii = 1:length(winds)
    plot(Vrest*1e3,FmaxMAT(winds(ii),:),'o-','color',[(1-ii/length(winds)),(1-ii/length(winds)),(1-ii/length(winds))])
end
hold off
    set(gca,'fontsize',16)
    xlabel('V_{rest,GC} (mV)');ylabel('LFP Fq (Hz)')
    xlim([-75 -55])
    ylim([10 90])
    
subplot(1,2,2)
hold on
for ii = 1:length(winds)
    plot(Vrest*1e3,maxpwrMAT(winds(ii),:),'o-','color',[(1-ii/length(winds)),(1-ii/length(winds)),(1-ii/length(winds))])
end
plot(Vrest*1e3,1.3e-4*ones(1,length(Vrest)),'k--')
hold off
    set(gca,'fontsize',16)
    xlabel('V_{rest,GC} (mV)');ylabel('Power')
    xlim([-75 -55])
    ylim([0 2.3e-3])
    
    
%% Part B
% runtime: immediate

% plot N-Type and NMDA activation curves
V = 1e3*((E_L-0.1):0.001:(E_L+0.1)); % in mV
v = V*1e-3; % in V

E_L = 0;

Mg_conc = 1;
E_Mg = 0;
gamma = 0.016;
eta = 0.28;
Mg_block = 1./(1+eta*Mg_conc*exp(-(v-E_Mg)/gamma));

scrsz = get(0,'ScreenSize');
figH=figure;
set(figH,'position',[0,400,scrsz(3)-0.6*scrsz(3),scrsz(4)-0.7*scrsz(4)]);
plot(V,mbarN,'-k',V,Mg_block,'.k');xlim([-80 0])
set(gca,'fontsize',17)
legend('N-Type','NMDA','location','southeast')
legend boxoff
xlabel('V_{rest,GC} (mV)');ylabel('Activation')


%% Parts C and D
% runtime: ~2hrs

% P1 controls ing=hibitory weight from GCs onto MCs, wGABAGR
P1.line = 44;
P1.name = 'wGABAGR ';
P1.val = 0.002:0.0005:0.026;

% P2 is Vrest, the GC excitability
P2.line = 56;
P2.name = 'Vrest ';
P2.val = 1e-3 * [-72 -60];


% Simulate the model over chosen parameter range
[MCMAT dGCMAT ICMAT param MClfpMAT dGClfpMAT] = ...
    ParamSweep_GCE(P1,P2,numtp,input_file);


%%% Calculate LFP frequency and power %%%

FmaxMAT = zeros(length(P1.val),length(P2.val));
maxpwrMAT = zeros(length(P1.val),length(P2.val));

trim = 500; % Ignore first 500 timepoints

L = length(timevec(trim:end-100));  % Length of data that will be analyzed
NFFT = 2^nextpow2(L); % Next power of 2 from L
f = sampf/2*linspace(0,1,NFFT/2+1);
ROI = ceil(8/(f(2)-f(1))):ceil(140/(f(2)-f(1)));
% ROI gives us a window from ~7Hz - 115Hz from which to select peak power.
% This avoids selecting low freuqency components that may dominate the
% spectrum.

% Run FFT on each simulated LFP
for n1 = 1:length(P1.val)
    for n2 = 1:length(P2.val)
        
        MitLFP = MClfpMAT{n1,n2}.GradistMitGlobal;
        mitFFT = fft(detrend(MitLFP(trim:end-100,1),'constant'),NFFT)/L;

        absmitFFT = 2*abs(mitFFT(1:NFFT/2+1));
        maxpwrMAT(n1,n2) = max(absmitFFT(ROI));

        maxind = find(absmitFFT == maxpwrMAT(n1,n2));
        FmaxMAT(n1,n2) = f(maxind);
    end
end


%%% PLot C and D %%%
scrsz = get(0,'ScreenSize');
figH=figure;
set(figH,'position',[0,400,scrsz(3)-0.65*scrsz(3),scrsz(4)-0.5*scrsz(4)]);
subplot(2,1,1)
plot(Weights,FmaxMAT(:,1),'b*-',Weights,FmaxMAT(:,2),'ro-')
set(gca,'fontsize',17,'YTickLabel',['20 ';'   ';'60 ';'   ';'100';'   ';'140'])
xlim([Weights(1) Weights(end)])
ylim([10 140])
ylabel('LFP fq (Hz)')
legend('V_{rest, GC} = -72 mV','V_{rest, GC} = -60 mV','location','northeast')
legend boxoff
subplot(2,1,2)
hold on
plot(Weights,maxpwrMAT(:,1),'b*-',Weights,maxpwrMAT(:,2),'ro-')
plot(Weights,1.3e-4*ones(1,length(Weights)),'k--')
hold off
xlim([Weights(1) Weights(end)])
ylim([0 2.3e-3])
set(gca,'fontsize',17)
xlabel('W_{GABA,MC}');ylabel('Power')

