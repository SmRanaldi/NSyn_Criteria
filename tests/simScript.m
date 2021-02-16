clear;close all; clc;

init='rand';
nReal=50;
nMuscles=15;
nResample=1000;
SNRtot={[2,5],[7 10],[15 18]};
fEnv=[2,5,10,20]./500;
% nSynergies=2:ceil(nMuscles/2)+1;
nSynergies=2:nMuscles;
nCycles=5; %Minimo 3
ada=0; %Ada=1 tests the adaptive algorithm

tTot=tic;

for nSynTrue=1:length(nSynergies)
    
    for realizationW=1:nReal
        
        for ss = 1:length(SNRtot)
            
            SNR = SNRtot{ss};
            
            disp(' ');
            
            disp(['SNR mean: ', num2str(mean(SNR))]);
            
            disp(' ');
            
            tstart=tic;
            
            [~,emgs,~,~,fs]=generateEMGSynergies(nMuscles,nSynergies(nSynTrue),SNR,200,1000*ceil(nSynergies(nSynTrue)/5),nCycles+2,500);
            newS = 1:nResample:size(emgs,2);
            newS = newS - 1;
            
            for fEnvelope=1:length(fEnv)+ada
                
                if fEnvelope < length(fEnv)+1
                    
                    [a,b]=butter(3,fEnv(fEnvelope),'low');
                    
                    clear M MProv;
                    
                    for i=1:nMuscles
                        
                        MProv(i,:)=filtfilt(a,b,abs(emgs(i,:)));
                        
                    end
                    
                else
                    
                    clear M MProv;
                    
                    for i=1:nMuscles
                        
                        MProv(i,:)=adaptiveEnvelope(conditionEMG(emgs(i,:)),'mincontrol',true);
                        
                    end
                    
                end
                M = MProv;
                
                for i=1:size(M,1)
                    M(i,M(i,:)<=0)=0.0001;
                end
                
                if max(isnan(M(:)))>0
                    disp(' ');
                end
                
                
                [~,~,~,~,AIC,~,VAFGlobal,R2Global,d,VAFMuscles]=synergiesAICWavelet(M,init,'min',newS);
                
                VAFTot(nSynTrue,fEnvelope,:,realizationW,ss)=VAFGlobal;
                VAFMTot(nSynTrue,fEnvelope,:,:,realizationW,ss)=VAFMuscles;
                R2Tot(nSynTrue,fEnvelope,:,realizationW,ss)=R2Global;
                AICTot(nSynTrue,fEnvelope,:,realizationW,ss)=AIC;
                [mm,am] = min(AIC);
                disp(['True value: ', num2str(nSynergies(nSynTrue)),'. Estimated value: ', num2str(am), '.']);
                
            end
            
        end
        
        ttt=toc(tstart);
        
        disp(' ');
        disp(['Simulated synergies: ', num2str(nSynergies(nSynTrue)), ' of ', num2str(max(nSynergies)), '. Realization ', num2str(realizationW), ' of ', num2str(nReal),'. Elapsed time: ', num2str(ttt),' s.']);
        disp([num2str( (((nSynTrue-1)*nReal+realizationW)/(nReal*max(nSynergies-1)))*100), '%']);
        
    end
    
    disp('-----------------------------------------------------------------------------------------------------');
    
end

timeEnd=now;
filename=['Simulation_',datestr(timeEnd,'HHMM-mmddyy')];

disp(['Saving ', filename, '...']);

save(filename);

tocTot=toc(tTot);

disp(['Done in ' num2str(tocTot), ' s.']);

clc;

% showResultsAIC(filename);

