% Code for generating Figure 5

clear all
input_file = 'OB_params_GCE_Fig5.txt';
[dt,tsim,numtp,nmit,ngradist,ngraprox,sampf,timevec] ...
    = InitNetwork_GCE(input_file);

%% Simulate LFP from MC IPSCs

% P1 controls inhibitory weight from GCs onto MCs, wGABAGR
P1.line = 44;
P1.name = 'wGABAGR ';
P1.val = 0.0165;
% Note: wGABAR is not varied here

% P2 is Vrest, the GC excitability
P2.line = 56;
P2.name = 'Vrest ';
P2.val = 1e-3 * [-73 -69 -60];

% Simulate the model over chosen parameter range
[MCMAT dGCMAT ICMAT param MClfpMAT dGClfpMAT] = ...
    ParamSweep_GCE(P1,P2,numtp,input_file);


%% Raster plots
% call to custom function PlotRasterplot.m

colorcell{1} = [0,0.1,1];
colorcell{2} = [0.5,0.1,0.5];
colorcell{3} = [1,0.1,0];

for ii = 1:3
    Rasterplot(MCMAT{ii},param,colorcell{ii},18);
    xlim([0 500])
end


%% LFP plots

colorcell{1} = [0,0.1,1];
colorcell{2} = [0.5,0.1,0.5];
colorcell{3} = [1,0.1,0];

scrsz = get(0,'ScreenSize');
for ii = 1:3
    figH = figure;
    set(figH,'position',[0,400,scrsz(3)-0.6*scrsz(3),scrsz(4)-0.8*scrsz(4)]);
    mitH = plot(timevec,MClfpMAT{ii}.GradistMitGlobal,'color',colorcell{ii});
    set(gca,'fontsize',18)
    xlim([0 500]);ylim([-.012 0])
    xlabel('time (ms)');ylabel('LFP')
end

%% Prelease plots

% GC graded activation (Prlelease)
scrsz = get(0,'ScreenSize');
for ii = 1:3
    figH=figure;
    set(figH,'position',[0,400,scrsz(3)-0.6*scrsz(3),scrsz(4)-0.8*scrsz(4)]);
    hold on
    for jj = 1:ngradist
        colorcell{1} = [0,0.1,jj/ngradist];
        colorcell{2} = [jj/(ngradist*2),0.1,jj/(ngradist*2)];
        colorcell{3} = [jj/ngradist,0.1,0];
        plot(timevec,ICMAT{ii}.Prelease(jj,:),'color',colorcell{ii})
    end
    hold off
    set(gca,'fontsize',18)
    xlim([0 500]);ylim([0 1])
    xlabel('time (ms)');ylabel('P_{release}')
end


%% Power spectra

trim = 1000; % trim beginning
L = length(timevec(trim:end-100));  % Length of simulation
NFFT = 2^nextpow2(L); % Next power of 2 from length of simulation
f = sampf/2*linspace(0,1,NFFT/2+1);

scrsz = get(0,'ScreenSize');
for ii = 1:3
    mitFFTG = fft(detrend(MClfpMAT{ii}.GradistMitGlobal(trim:end-100,1),'constant'),NFFT)/L;
    mitFFT1 = fft(detrend(MClfpMAT{ii}.GradistMit1(trim:end-100,1),'constant'),NFFT)/L;
    mitFFT2 = fft(detrend(MClfpMAT{ii}.GradistMit2(trim:end-100,1),'constant'),NFFT)/L;
    mitFFT3 = fft(detrend(MClfpMAT{ii}.GradistMit3(trim:end-100,1),'constant'),NFFT)/L;

    figH=figure;
    set(figH,'position',[0,400,scrsz(3)-0.7*scrsz(3),scrsz(4)-0.7*scrsz(4)]);
    hold on
    mitH = plot(f,2*abs(mitFFTG(1:NFFT/2+1)),f,2*abs(mitFFT1(1:NFFT/2+1)),f,2*abs(mitFFT2(1:NFFT/2+1)),f,2*abs(mitFFT3(1:NFFT/2+1)));
    plot(f,2.6e-4*ones(1,length(f)),'k--')
    hold off
    xlim([1 100]);ylim([0 4.5e-3])
    set(gca,'fontsize',15)
    legend('G','1','2','3')
    legend boxoff
    xlabel('Frequency (Hz) ');ylabel('Power')
end




