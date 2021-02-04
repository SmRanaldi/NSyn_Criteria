function showResultsAIC(filename)

clc;

load(filename);

disp(['Simulation ended ', datestr(timeEnd,'HH:MM dd/mm/yyyy'),'.']);
disp(['Number of realizations: ', num2str(nReal),'.']);
disp(['Number of muscles: ', num2str(nMuscles),'.']);
disp(['Samples per cycle: ', num2str(newS),'.']);
disp(['SNR limits: ', num2str(min(SNR)),'-',num2str(max(SNR)),' dB.']);
disp(['Number of synergies from ', num2str(min(nSynergies)),' to ',num2str(max(nSynergies)),'.']);
disp(['Frequencies from ', num2str(min(fEnv)*(500)),' to ',num2str(max(fEnv)*500),'.']);

AICTotMean=nanmedian(AICTot,4);
[~,~,~,~,~,~,VAFRand]=synergiesAICWavelet(rand(nMuscles,lengthCycle),'rand','min',[0,lengthCycle]);
for i=1:length(nSynergies)
    for j=1:size(VAFTot,2)
        for k=1:nReal
            [~,idxm]=min(AICTot(i,j,:,k),[],3);
            nSynAIC(i,j,k)=idxm;
            nSynVAFThs(i,j,k)=nSynThs(squeeze(VAFTot(i,j,:,k)),0.9);
            nSynVAFMus(i,j,k)=nSynThs(squeeze(VAFMTot(i,j,:,:,k)),0.9);
            nSynVAFThs1(i,j,k)=nSynThs(squeeze(VAFTot(i,j,:,k)),0.95);
            nSynVAF5Perc(i,j,k)=nSyn5Perc(squeeze(VAFTot(i,j,:,k)));
            nSynVAFRand(i,j,k)=nSynRand(squeeze(VAFTot(i,j,:,k)),VAFRand);
        end
    end
end
VAFTotMean=nanmedian(VAFTot,4);

for i=1:length(nSynergies)
    
    for j=1:size(VAFTot,2)
        
        eAIC(i,j,:)=nSynAIC(i,j,:)-nSynergies(i);
        eVAFThs(i,j,:)=nSynVAFThs(i,j,:)-nSynergies(i);
        eVAFMus(i,j,:)=nSynVAFMus(i,j,:)-nSynergies(i);
        eVAFThs1(i,j,:)=nSynVAFThs1(i,j,:)-nSynergies(i);
        eVAF5Perc(i,j,:)=nSynVAF5Perc(i,j,:)-nSynergies(i);
        eVAFRand(i,j,:)=nSynVAFRand(i,j,:)-nSynergies(i);
        
    end
    
end

figure;

B=0.5:1:nMuscles+.5;

for i=1:length(nSynergies)
    for j=1:size(nSynAIC,2)
        subplot(size(nSynAIC,1),size(nSynAIC,2),(i-1)*size(nSynAIC,2)+j);
        histogram(squeeze(nSynAIC(i,j,:)),B);
        hold on;
        plot([mode(squeeze(nSynAIC(i,j,:))), mode(squeeze(nSynAIC(i,j,:)))],[0,length(squeeze(nSynAIC(i,j,:)))],'r','linewidth',2);
        plot([mean(squeeze(nSynAIC(i,j,:))), mean(squeeze(nSynAIC(i,j,:)))],[0,length(squeeze(nSynAIC(i,j,:)))],'g','linewidth',2);
        title(['nSyn = ', num2str(i+1)]);
    end
end

if ~exist('nSynVAFThs')
    
    for i=1:length(nSynergies)
        for j=1:size(nSynAIC,2)
            meanNSyn(i,j)=median(nSynAIC(i,j,:));
            stdNSyn(i,j)=sqrt(sum(((squeeze(nSynAIC(i,j,:)))-meanNSyn(i,j)).^2)/(length(nSynAIC(i,j,:)-1)));
        end
    end
    
    figure;
    
    for i=1:length(nSynergies)
        sPlot(i)=subplot(1,length(nSynergies),i);
        barwitherr(stdNSyn(i,:),meanNSyn(i,:));
        hold on;
        plot([0,4],[nSynergies(i),nSynergies(i)],'k','linewidth',2);
    end
    linkaxes(sPlot,'xy');
    
else
    
    for i=1:length(nSynergies)
        
        for j=1:size(nSynAIC,2)
            
            meanNSyn(i,j)=median(nSynAIC(i,j,:));
            stdNSyn(i,j)=sqrt(sum(((squeeze(nSynAIC(i,j,:)))-meanNSyn(i,j)).^2)/(length(nSynAIC(i,j,:)-1)));
            pAIC(i,j)=(sum(eAIC(i,j,:)==0)/length(eAIC(i,j,:)))*100;
            
            meanNSyn5Perc(i,j)=median(nSynVAF5Perc(i,j,:));
            stdNSyn5Perc(i,j)=sqrt(sum(((squeeze(nSynVAF5Perc(i,j,:)))-meanNSyn5Perc(i,j)).^2)/(length(nSynVAF5Perc(i,j,:)-1)));
            pVAF5Perc(i,j)=(sum(eVAF5Perc(i,j,:)==0)/length(eVAF5Perc(i,j,:)))*100;
            
            meanNSynThs(i,j)=median(nSynVAFThs(i,j,:));
            stdNSynThs(i,j)=sqrt(sum(((squeeze(nSynVAFThs(i,j,:)))-meanNSynThs(i,j)).^2)/(length(nSynVAFThs(i,j,:)-1)));
            pVAFThs(i,j)=(sum(eVAFThs(i,j,:)==0)/length(eVAFThs(i,j,:)))*100;
            
            meanNSynThs1(i,j)=median(nSynVAFThs1(i,j,:));
            stdNSynThs1(i,j)=sqrt(sum(((squeeze(nSynVAFThs1(i,j,:)))-meanNSynThs1(i,j)).^2)/(length(nSynVAFThs1(i,j,:)-1)));
            pVAFThs1(i,j)=(sum(eVAFThs1(i,j,:)==0)/length(eVAFThs1(i,j,:)))*100;
            
            meanNSynRand(i,j)=median(nSynVAFRand(i,j,:));
            stdNSynRand(i,j)=sqrt(sum(((squeeze(nSynVAFRand(i,j,:)))-meanNSynRand(i,j)).^2)/(length(nSynVAFRand(i,j,:)-1)));
            pVAFRand(i,j)=(sum(eVAFRand(i,j,:)==0)/length(eVAFRand(i,j,:)))*100;
            
            meanNSynMus(i,j)=median(nSynVAFMus(i,j,:));
            stdNSynMus(i,j)=sqrt(sum(((squeeze(nSynVAFMus(i,j,:)))-meanNSynMus(i,j)).^2)/(length(nSynVAFMus(i,j,:)-1)));
            pVAFMus(i,j)=(sum(eVAFMus(i,j,:)==0)/length(eVAFMus(i,j,:)))*100;
            
        end
        
    end
    
    figure;
    
    labels={'AIC','90%','95%','5%','Rand'};
    
    for i=1:length(nSynergies)
        sPlot(i)=subplot(2,ceil(length(nSynergies)/2),i);
        barwitherr([stdNSyn(i,:),0,stdNSynThs(i,:),0,stdNSynThs1(i,:),0,stdNSyn5Perc(i,:),0,stdNSynRand(i,:),0,stdNSynMus(i,:)],[meanNSyn(i,:),0,meanNSynThs(i,:),0,meanNSynThs1(i,:),0,meanNSyn5Perc(i,:),0,meanNSynRand(i,:),0,meanNSynMus(i,:)]);
        hold on;
        plot([0,24],[nSynergies(i),nSynergies(i)],'r','linewidth',2);
        axis([0 24 0 max(nSynergies)+1]);
%         xticks([2,6,10,14,18]);
%         xticklabels(labels);
    end
    linkaxes(sPlot,'xy');
    
    figure;
    
    for i=1:length(nSynergies)
        sPlot(i)=subplot(2,ceil(length(nSynergies)/2),i);
        bar([pAIC(i,:),0,pVAFThs(i,:),0,pVAFThs(i,:),0,pVAF5Perc(i,:),0,pVAFRand(i,:),0,pVAFMus(i,:)])
        axis([0 25 0 100]);
%         xticks([2,6,10,14,18]);
%         xticklabels(labels);
    end
    linkaxes(sPlot,'xy');
    
    figure;
    
    for i=1:length(fEnv)
        aaa(i)=subplot(length(fEnv),1,i);
        plot(pAIC(:,i),'k','linewidth',2);
        hold on;
        plot(pVAFThs(:,i),'c','linewidth',2);
        plot(pVAFThs1(:,i),'r','linewidth',2);
        plot(pVAF5Perc(:,i),'g','linewidth',2);
        plot(pVAFRand(:,i),'b','linewidth',2);
        plot(pVAFMus(:,i),'y','linewidth',2);
    end
    linkaxes(aaa,'xy');
    
end