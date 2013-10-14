function test_simulation(jobID, nMaxRobots, nSimulations)
% main function
% compute a path for the mobile robot in function of the kriging error
close all;

% a path of nWPpath waypoints is recomputed each time the robot reaches
% nWayPoints waypoints
nWayPoints=5;% number of waypoints to reach before recomputing the path
nWPpath=10;% number of waypoints in the path

% choose randomly a random field
fieldNum= randi([1 100]);

if mod(jobID, 3)== 1
    field=load(['./RandomFields/RandField_LR_No' num2str(200+fieldNum) '.csv']);
elseif mod(jobID,3)== 2
    field=load(['./RandomFields/RandField_IR_No' num2str(100+fieldNum) '.csv']);
else
    field=load(['./RandomFields/RandField_SR_No' num2str(fieldNum) '.csv']);
end

%we want a field of 200x200 m. (not 199x199)
field(:,end+1)=field(:,end);
field(end+1,:)=field(end,:);
[Ly,Lx]=size(field);

%position of the stations (static sensors)
stations=[];
stations(:,1)=[24 34 14 94 134 74 94 166 186 174];
stations(:,2)=[166 94 22 14 86 66 174 174 106 34];

%grid used to save the measures
h = 1;
x = 0:h:Lx-1;
y = 0:h:Ly-1;
lx=length(x);
ly=length(y);

%interpolation points (points where the kriging error is computed)
delta=5;
x_=0:delta:Lx-1;
y_=0:delta:Ly-1;

%allowable waypoints position
ph = 20;
px = 0:ph:Lx-1;
py = 0:ph:Ly-1;
lpx=length(px);
lpy=length(py);

RMSEValuesAllScenarios= [];

for nRobots= 1: nMaxRobots
    allRMSESameNumberOfRobots= [];
    for currentSimulation= 1: nSimulations
        %vector containing the sampling points (nan at position non-sampled)
        sVec=ones(lx*ly,1)*nan;
        kVec=1:lx*ly;%indices
        
        sVec=addSamplingPoints(sVec,stations,field,x,y,lx);
        
        %parameter of each robot
        parameters = repmat(struct('speedHeli',3.7,'measPeriod',3,'trendOrder',0, 'pos', 0, 'pathFollowed', [], 'sampledPoints', [] ), nRobots, 1 );
        %make sure to generate different initial positions
        possibilities=lpx*lpy;
        r=randperm(possibilities);
        r=r(1:nRobots);
        
        for i=1:nRobots
            parameters(i).pos= r(i);
        end
        
        iter=1;
        RMSE=[];%roots mean square error
        distances= [];
        distances(1:nRobots,1)= 0;
        h0(1:nRobots,1)=0;
        
        %% loop while the travelded distance is smaller than 7000 m
        while sum(distances(:,end))<3000
            for currentRobot=1:nRobots
                %----------------------------------------------------------------------
                %get sampling points
                X=[];
                indK=kVec(~isnan(sVec)); % get the indices of the sampling points
                indY = ceil(indK/lx);
                indX = indK -(indY-1)*lx;
                Y=sVec(indK);   % sampled values
                X(:,1)=x(indX); % sampling points position
                X(:,2)=y(indY); % sampling points position
                
                %----------------------------------------------------------------------
                %standardisation of the data
                meanV=mean(Y);
                stdV=std(Y);
                Y_=(Y-meanV)/stdV;
                
                %----------------------------------------------------------------------
                %compute variogram
                % the variogam is computed at discrete positions separated by a
                % distance called lag
                lag=10;
                range=130; % range considered to compute the variogram
                [fittedModel,fittedParam]=variogram(X,Y_,lag,range,1);
                
                %----------------------------------------------------------------------
                %kriging interpolation
                [val,krigE_]=kriging(X,Y_,stdV,meanV,fittedModel,fittedParam,parameters(currentRobot).trendOrder,x_,y_);
                if currentRobot== nRobots
                    RMSE(iter) = sqrt(mean(mean((val-field(1:delta:lx,1:delta:ly)).^2)));
                end
                
                %----------------------------------------------------------------------
                %compute path
                [Dh,Dv,Ddu,Ddd]=distanceMatrix(x_,y_,krigE_,px,py);
                [path,~]=findBestPath(px,py,parameters(currentRobot).pos,Dh,Dv,Ddu,Ddd,nWPpath,0.6,4.4);
                
                nP=min(nWayPoints+1,length(path));
                path=path(1:nP,:);
                
                parameters(currentRobot).pathFollowed= [parameters(currentRobot).pathFollowed; path];
                %parameters(currentRobot).pathFollowed=
                %unique(parameters(currentRobot).pathFollowed,
                %'rows','stable');     commented because on seg is not
                %working
                
                indX = floor(path(end,1)/ph)+1;
                indY = floor(path(end,2)/ph)+1;
                parameters(currentRobot).pos= indX + (indY-1)*lpx;
                
                [Pts2visit,distances(currentRobot,iter+1),h0(currentRobot,1)] = findPtsAlongPath(path, parameters(currentRobot).speedHeli,parameters(currentRobot).measPeriod,distances(currentRobot, iter),h0(currentRobot,1));
                
                parameters(currentRobot).sampledPoints= [parameters(currentRobot).sampledPoints; Pts2visit];
                
                %----------------------------------------------------------------------
                %add new sampling points
                %         if iter==1
                %             sVec=addSamplingPoints(sVec,[0,0],field,x,y,lx);
                %         end
                sVec=addSamplingPoints(sVec,Pts2visit,field,x,y,lx);
                %----------------------------------------------------------------------
            end
            iter=iter+1;
        end
        
        dist=sum(distances(:,1:end-1), 1);
        
        RMSE_=interp1(dist,RMSE,0:50:2900);
        allRMSESameNumberOfRobots=[allRMSESameNumberOfRobots  RMSE_'];
    end
        
    RMSEValuesAllScenarios=[RMSEValuesAllScenarios mean(allRMSESameNumberOfRobots,2)];
end
if ~exist('./results', 'dir')
    mkdir('./results');
end
FileName= strcat('./results/simulationResultJob_', num2str(jobID), '.mat');
save( FileName, 'RMSEValuesAllScenarios', 'jobID');

end


%==========================================================================
function [Dh,Dv,Ddu,Ddd]=distanceMatrix(x_,y_,krigE,x,y)
% create the matrices that contains the disance for the ACO algorithm
% input: x_,y_          grid positions
%        krigE          Value of the krigE at the grid position
%        x,y            allowable waypoints positions

Dh=zeros(length(y),length(x)-1);
Dv=zeros(length(y)-1,length(x));
Ddu=zeros(length(y)-1,length(x)-1);
Ddd=zeros(length(y)-1,length(x)-1);

n=10;
dh=(x(2)-x(1))/n;
%val=bilinInterp(f,x,y,x0,y0,lx_1,ly_1);
lx_1=length(x_)-1;
ly_1=length(y_)-1;
for i=1:length(y)
    for j=1:length(x)-1
        %Dh(i,j)=1.5-(krigE(i,j)+krigE(i,j+1));
        for k=0:n
            Dh(i,j) = Dh(i,j) + bilinInterp(x_,y_,krigE,x(j)+k*dh,y(i),lx_1,ly_1);
        end
        %Dh(i,j) = 1.5 - (Dh(i,j)/5);
        Dh(i,j) = Dh(i,j)/(n+1);
    end
end


for i=1:length(y)-1
    for j=1:length(x)
        %Dv(i,j)=1.5-(krigE(i,j)+krigE(i+1));
        for k=0:n
            Dv(i,j) = Dv(i,j) + bilinInterp(x_,y_,krigE,x(j),y(i)+k*dh,lx_1,ly_1);
        end
        %Dv(i,j) = 1.5 - (Dv(i,j)/5);
        Dv(i,j) = Dv(i,j)/(n+1);
    end
end

for i=1:length(y)-1
    for j=1:length(x)-1
        %Dv(i,j)=1.5-(krigE(i,j)+krigE(i+1));
        for k=0:n
            Ddu(i,j) = Ddu(i,j) + bilinInterp(x_,y_,krigE,x(j)+k*dh,y(i)+k*dh,lx_1,ly_1);
            Ddd(i,j) = Ddd(i,j) + bilinInterp(x_,y_,krigE,x(j)+k*dh,y(i+1)-k*dh,lx_1,ly_1);
        end
        %Dv(i,j) = 1.5 - (Dv(i,j)/5);
        Ddu(i,j) = Ddu(i,j)/(n+1);
        Ddd(i,j) = Ddd(i,j)/(n+1);
    end
end

%norm(p(I,:)-p(J,:))*(0.2+E/30);
end

function plotmap(x,y,map)
% 2d plot
figure
[~, ch]=contourf(x,y,map,30);
set(ch,'edgecolor','none');
set(gca,'FontSize',16)
%xlabel('X [m]','FontSize',14)
%ylabel('Y [m]','FontSize', 14)
%colorbar
axis('equal')
axis([-3 203 -3 203])
end

function sVec=addSamplingPoints(sVec,X,field,x,y,lx)
% addSamplingPoints
% input: sVec           vector the contain the sampling values
%        X              position to  be sampled
%        field          field value
%        x,y            field positions
%        lx             length of the field

for i=1:length(X(:,1))
    indX = round(X(i,1))+1;
    indY = round(X(i,2))+1;
    indK = indX+(indY-1)*lx;
    sVec(indK) = field(y(indY)+1,x(indX)+1);
end
end

function val=bilinInterp(x,y,f,x0,y0,lx_1,ly_1)
% perform bilinear interpolation
% input: x,y            grid positions
%        f              values at grid points
%        x0,y0          position of interest
%        lx_1,ly-1      length of the grid - 1
% ouput: val            interpolation value at the position of interest

h=x(2)-x(1);
Ix1 = round(x0/h);
Iy1 = round(y0/h);
if (Ix1 < 1)
    Ix1 = 1;
elseif (Ix1 > lx_1)
    Ix1 = lx_1;
end

if (Iy1 < 1)
    Iy1 = 1;
elseif (Iy1 > ly_1)
    Iy1 = ly_1;
end

Ix2=Ix1+1;
Iy2=Iy1+1;
DeltaX1 = x0 - x(Ix1);
DeltaX2 = x0 - x(Ix2);
DeltaY1 = y0 - y(Iy1);
DeltaY2 = y0 - y(Iy2);
val=1/(h*h)*(f(Iy1,Ix1)*DeltaX2*DeltaY2 - f(Iy1,Ix2)*DeltaX1*DeltaY2 - f(Iy2,Ix1)*DeltaX2*DeltaY1 + f(Iy2,Ix2)*DeltaX1*DeltaY1);

end