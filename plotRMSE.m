
function plotRMSE(nRobots)
close all
fileList=dir('./results/simulationResultJob_*.mat');
shortRange=cell(nRobots);
mediumRange=cell(nRobots);
longRange=cell(nRobots);
%assume L.currentMSE is a column vector

for i=1:length(fileList)
    L=load(strcat('./results/', fileList(i).name));
    
    if mod(L.jobID, 3)== 1
        shortRange{L.nRobots}(:, size(shortRange{L.nRobots}, 2)+1)= L.currentRMSE;
    elseif mod(L.jobID, 3)== 2
        mediumRange{L.nRobots}(:, size(mediumRange{L.nRobots}, 2)+1)= L.currentRMSE;
    else
        longRange{L.nRobots}(:, size(longRange{L.nRobots}, 2)+1)= L.currentRMSE;
    end
end

shortRangePlotY=[];
mediumRangePLotY=[];
longRangePlotY=[];
totalPlotY=[];
legendNames=[];

for i=1:nRobots
    shortRangePlotY= [shortRangePlotY mean(shortRange{i},2)];
    mediumRangePLotY= [mediumRangePLotY mean(mediumRange{i},2)];
    longRangePlotY= [longRangePlotY mean(longRange{i},2)];
    
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