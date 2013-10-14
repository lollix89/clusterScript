
function plotRMSE()
close all
fileList=dir('./results/simulationResultJob_*.mat');
shortRange=[];
mediumRange=[];
longRange=[];

NRobots=0;

for i=1:length(fileList)
    L=load(strcat('./results/', fileList(i).name));
    NRobots= size(L.RMSEValuesAllScenarios, 2);
        
    for j=1:size(L.RMSEValuesAllScenarios, 2)        %#robots in the simulation
        
        if mod(L.jobID, 3)== 1
            shortRange(j,:,end+1)= L.RMSEValuesAllScenarios(:,j);
        elseif mod(L.jobID, 3)== 2
            mediumRange(j,:,end+1)= L.RMSEValuesAllScenarios(:,j);
        else
            longRange(j,:,end+1)= L.RMSEValuesAllScenarios(:,j);
            
        end
    end
end

shortRangePlotY=[];
mediumRangePLotY=[];
longRangePlotY=[];
totalPlotY=[];
legendNames=[];

for i=1:NRobots 
    shortRangePlotY= [shortRangePlotY mean(shortRange(i,:,2:end),3)'];
    mediumRangePLotY= [mediumRangePLotY mean(mediumRange(i,:,2:end),3)'];
    longRangePlotY= [longRangePlotY mean(longRange(i,:,2:end),3)'];
    
    totalPlotY=[totalPlotY mean([shortRangePlotY(:,end) mediumRangePLotY(:,end) longRangePlotY(:,end)],2)];
    
    legendNames=[legendNames; strcat(num2str(i),' robots')];
end

plotRangeX = (0:50:2900)';

shortFigure= figure();
plot(plotRangeX,shortRangePlotY)
title('shortRange RandomField comparison')
ylabel('RMSE')
xlabel('Distance (m)')
grid on
legend(legendNames)
saveas(shortFigure,'shortRange','png')

mediumFigure=figure();
plot(plotRangeX,mediumRangePLotY)
title('mediumRange RandomField comparison')
ylabel('RMSE')
xlabel('Distance (m)')
grid on
legend(legendNames)
saveas(mediumFigure,'mediumRange','png')


longFigure=figure();
plot(plotRangeX,longRangePlotY)
title('longRange RandomField comparison')
ylabel('RMSE')
xlabel('Distance (m)')
grid on
legend(legendNames)
saveas(longFigure,'longRange','png')

totalFigure=figure();
plot(plotRangeX,totalPlotY)
title('total comparison')
ylabel('RMSE')
xlabel('Distance (m)')
grid on
legend(legendNames)
saveas(totalFigure,'total','png')


end