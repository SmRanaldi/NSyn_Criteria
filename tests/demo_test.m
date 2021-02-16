clearvars -except sw par;
close all;
clc;


sw=1;
par = input('Par: ');

[a,b]=butter(4,10/500,'low');

switch sw
    
    case 1
        
        [~,~,~,EMG,idx]=genEMGSyn(12,200,par,20,[20,30]);
        for i=1:size(EMG,1)
            env(i,:)=filtfilt(a,b,abs(EMG(i,:)));
        end
        env(env<0.001)=0.001;
        M=env(:,idx(2):idx(end-1));
        idx=idx-idx(2)+1;
        idx=idx(2:end-1);
        [W,H,meanH,VAF,AIC,R2,VAFGlobal,R2Global,DoFTot,VAFMuscles,DoF,L,nPower,E]=synergiesAICWavelet(env(:,1:end),'rand','min',idx);
      [tmptmp,nSyn] = min(AIC);
      figure;plot(AIC);hold on;
%         plot(L,'k');plot(DoF,'r');
        plot(nSyn,AIC(nSyn),'g.','markersize',20);
%         plot(L,'k');plot(DoF,'r');
        %         figure;
        %         for i=1:12
        %             subplot(4,3,i)
        %             plot(EMG(i,1:lCycle));
        %         end
        
    case 2
        
        nSamples=128;
        fe=par;
        [emgs,~,~,idx,envTemp]=pedallingEMG('SOGGETTO1.mat',nSamples,20,fe,0);
        
        
      [W,H,meanH,VAF,AIC,R2,VAFGlobal,R2Global,DoFTot,VAFMuscles,DoF,L,nPower,E]=sAICREvents(envTemp(:,1:end),'rand','min',idx);
      [tmptmp,nSyn] = min(AIC);
      figure;plot(AIC);hold on;
%         plot(L,'k');plot(DoF,'r');
        plot(nSyn,AIC(nSyn),'g.','markersize',20);
        
        
    case 3
        
        f0=10;
        f1=20;
        c=1;
        x=zeros(1,500);
        x(f0:f1)=hanning(f1-f0+1)';
        %         x(6+a+c:6+a+c+b-1)=0.9*hanning(b)';
        clc
        x=x+0.01*rand(size(x))+0.0001;
        plot(x,'k','linewidth',2)
        out=maxFCorr(x);hold on;plot(out,x(out),'r.','markersize',20)
        
    case 4
        
        load('Dati_Dublino.mat');
        L=200;s=4;fe=par;
        [envsOut,idx,envs]=emgDublin(EMG{s},fs_emg,BIOMECH{s},fs_biomech,fe,L);
        [W,~,~,AIC,nSyn,~,VAF,~,~,~,DoF,L]=sAICREvents(envs,'rand','min',idx);
        lAIC=AIC(1:end-1);L=L(1:end-1);DoF=DoF(1:end-1);
        plot(AIC);hold on;
        plot(nSyn,AIC(nSyn),'g.','markersize',20);
%         plot(L,'k');plot(DoF,'r')
        
end