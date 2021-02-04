clear;close all; clc;

init='sparse';
nReal=2;
nMuscles=8;
nResample=1000;
SNR=[10 30];
fEnv=[5, 10, 15, 30]./500;
% nSynergies=2:ceil(nMuscles/2)+1;
nSynergies=2:nMuscles-3;
nCycles=3; %Minimo 3
lengthCycle=1000;
ada=0; %Ada=1 tests the adaptive algorithm

tTot=tic;

for nSynTrue=1:length(nSynergies)
    
    for realizationW=1:nReal
        
        tstart=tic;
        
        [~,emgs,~,~,fs]=generateEMGSynergies(nMuscles,nSynergies(nSynTrue),SNR,400,1000*ceil(nSynergies(nSynTrue)/5),nCycles+2,nResample);
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
            
            VAFTot(nSynTrue,fEnvelope,:,realizationW)=VAFGlobal;
            VAFMTot(nSynTrue,fEnvelope,:,:,realizationW)=VAFMuscles;
            R2Tot(nSynTrue,fEnvelope,:,realizationW)=R2Global;
            AICTot(nSynTrue,fEnvelope,:,realizationW)=AIC;
           
            
        end
        
        ttt=toc(tstart);
        
        disp(['Simulated synergies: ', num2str(nSynTrue+1), ' of ', num2str(max(nSynergies)), '. Realization ', num2str(realizationW), ' of ', num2str(nReal),'. Elapsed time: ', num2str(ttt),' s.']);
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

showResultsAIC(filename);

