fs=97656.25; %TDT sampling rate

tvec=1/fs:1/fs:1; %time in s
tdat=tvec;
tvec2=1/fs:1/fs:.5; %duration of each piece
rftime=.001; %rise fall time in s
rflen=round(rftime*fs); %length of rise or fall
%flist=[2000 8000 19000 32000]; %frequency list in Hz
flist=8000;%tone carrier in Hz
rise(1:rflen)=1-(cos(2*pi*(1/(1.33*pi*rftime))*tvec(1:rflen))).^2;
fall(1:rflen)=rise(rflen:-1:1);
noisecar=randn(1,length(tvec)); %noise carrier
tonecar=sin(2*pi*flist*tvec); %tone carrier
noisecar(find(noisecar>3.5))=3.5; %remove large outliers in noise
noisecar(find(noisecar<-3.5))=-3.5; %remove large outliers in noise
% AMmodfq1=linspace(3.5,10.5,length(tvec2));
AMmodfq1(1:length(tvec2))=linspace(4,6.5,length(tvec2));
AMmodfq1(length(tvec2)+1:length(tvec))=linspace(6.5,10.25,length(tvec2));
%Made AMmodfq1 a single continuous vector of tvec length
 
AMmoddepth1=1;
AMfreqvec=2.^AMmodfq1;
% AMmod1=.5*(1-AMmoddepth1*cos(2*pi*AMfreqvec.*tvec+.001));
AMmod1=.5*(1-AMmoddepth1*cos(AMfreqvec.*tvec+.001));
%Made AMmod1 a single continuous vector of tvec length

% amfmnoisevec(1:rflen)=amfmnoisevec(1:rflen).*rise;
% amfmnoisevec((length(tvec2)-rflen+1):length(tvec2))=amfmnoisevec((length(tvec2)-rflen+1):length(tvec2)).*fall;
amfmnoisevec(1:length(tvec))=noisecar.*AMmod1;
amfmtonevec=tonecar(1:length(tvec)).*AMmod1;
%%Made amfmtonevec a single continuous vector of tvec length

% amfmtonevec(1:rflen)=amfmtonevec(1:rflen).*rise;
% amfmtonevec(length(tvec2)+1:length(tvec))=amfmtonevec(length(tvec2):-1:1);

subplot(2,1,1)
spectrogram(amfmtonevec,1200,1000,[],fs,'yaxis')
ylim([6.5 9.5])
subplot(2,1,2)
plot(tvec,amfmtonevec,'b-')
