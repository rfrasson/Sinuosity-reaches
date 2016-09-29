%River='Sacramento'; 
River='PoDS'; %downstream section of the Po River
%River='Po'; %upstream section of the Po River
pathtodata=['./RawData/' River '/'];
MinReachLen=5; 
tcritReach=2;
numbregresspts=10;
%sets default configurations for Po and Sacramento rivers
if strfind(River,'Sacramento')
    SWATHboundaries= [0 116927.0128 151838.8427]/1000; %division by 1000 to translate to km
    ReferenceDay='170'; %day used for the definition of reach boundaries and for the estimation of A0
    %Actual day user wants to compute reach-averaged values
    Day='170'; %23 is the lowest flow from the series. 86 and 170 are intermediate flow and 65 is the highest flow, 128 is high flow, 
    %available days are: 2 23 44 53 65 86 107 128 149 170 
    m=50;    
    W=160;  %average river width for high flow
    filenameReference=[pathtodata River 'Day' ReferenceDay '.mat'];
    filenameDataset=[pathtodata River 'Day' Day '.mat'];   
end
if strfind(River,'Po')
    ReferenceDay='220'; %day used to trace the center line and to estimate A0
    RefOverpass='560';
    SWATHboundaries=[];
    m=10;    
    W=580;  %average river width for high flow
    %For overpass 560: available days are: '136' '157' '178' '199' '220' '241' '262' '283' '304' '325' '346' '367' '388' '409'
    %For overpass 211: available days are: '145' '166' '187' '208' '229' '250' '271' '292' '313' '334' '355' '376' '397' '418' '439' '460' '481'
    Day='157';
    Overpass='560';
    filenameReference=[pathtodata River 'Day' ReferenceDay '-' RefOverpass '.mat'];
    filenameDataset=[pathtodata River 'Day' Day '-' Overpass '.mat'];
end
lambda=m*W; %Kiel's thesis says this should be 10 * bankfull width

Makeplots=1;
OutputPath='./output/';
OutFileName=[River 'AutomatedSinuosity' Day '.mat'];
SaveResults=1;
SmoothData=1;
VariableSmoothingWindow=1; 
load(filenameReference);
RefRiverObs=RiverObs;
RefTrue=True;

[ReachBoundaries,ReachLength]=FindSinuosityReaches(RefRiverObs,SWATHboundaries,lambda,MinReachLen,tcritReach,numbregresspts,Makeplots);

load(filenameDataset,'RiverObs','True');
[Reach,RiverData,Metadata,ReachTrue,Nodes,NodesTrue]=ReachAveraging(ReachBoundaries, RiverObs,True,RefRiverObs,RefTrue,ReachLength,Day,SaveResults,SmoothData,VariableSmoothingWindow,OutputPath, OutFileName,Makeplots);