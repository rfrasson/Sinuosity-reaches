function [ReachBoundaries,ReachLength]=FindSinuosityReaches(RiverObs,SWATHboundaries,lambda,MinReachLen,tcritReach,numbregresspts,Makeplots)
%This function defines reach boundaries based on river sinuosity and SWOT
%boundaries.
%List of inputs
%RiverObs nodes    : Structure containing RiverObs nodes (see explanation
%                    below)
%SWATHboundaries   : Location of the intersections between SWATH boundaries
%                    and the river centerline as flow distances in km.
%lambda            : Distance between endpoints for the calculation of
%                    sinuosity in meters. (see explanation below) 
%MinReachLen       : Minimum reach length in km
%tcritReach        : Number os Standard Deviations used to create
%                    dS/dx confidense intervals. Tested value = 2
%numbregresspts    : Number of regression points used to estimate dS/dx
%                    tested value = 10
%Makeplots         : Generate reach boundary plots if Makeplots =1

%list of outputs:
%ReachBoundaries   : Indices of the nodes that correspond to reach
%                    boundaries
%ReachLength       : Length of each reach


%Clarification on the format of RiverObs nodes.
%The variable passed as RiverObs must be a structure containing the
%following elements (as vectors with dimenstion eaqual to Number of nodes,1:
%x          : flow distance measured in meters
%H          : Node elevation in meters
%W          : River width in meters
%Easting    : Projected coordinate of the node in m
%Northing   : Projected coordinate of hte node in m
%Lat        : Latitude of the node in decimal degrees
%Lon        : Longitude of the node in decimal degrees
%xtrack     : cross-track position of the node

%Calculation of Lambda : lambda=m*W where M should be 10 () and W is the
%bankfull width in meters. This number is too small if the river is narrow
%(~150 m as the Sacramento) for the node spacing that we tried (250m).
%~30 nodes is the minimum we recomend (for the Sacramento, m = 50 and W =
%160)

    X=RiverObs.Easting;  
    Y=RiverObs.Northing; 
    %% distance calculation
    n=length(X);
    D=nan(size(X));
    D(1)=0;

    for i=2:n,
        D(i)=D(i-1)+ sqrt( (X(i)-X(i-1))^2 + (Y(i)-Y(i-1))^2 );
    end


    %% sinuosity calculation
    S=zeros(size(D));
    for i=1:n,

        jfwd=D>=D(i) & D<D(i)+lambda/2;
        jbak=D<D(i) & D>D(i)-lambda/2;
        j=jfwd | jbak ;

        x=X(j); y=Y(j); d=D(j);

        De=sqrt( (x(1)-x(end))^2 + (y(1)-y(end))^2 );
        Dr=d(end)-d(1);

        S(i)=Dr/De;
    end
    %seach for inflection points in sinuosity vs flow distance
    [ReachBoundaries,~,~,ReachLength]=DetectCurvatureLinear(D/1000,S,0,1,numbregresspts,tcritReach);
    
    FlowDist=D/1000; %converts flow distance to km for further use

    %find nodes that are closest to the swath boundaries
    SwathBoundIDs=zeros(size(SWATHboundaries));
    NodeID=1;
    for countSwath=1:length(SWATHboundaries);
        countnode=NodeID;
        while countnode<length(FlowDist)&& SWATHboundaries(countSwath)>FlowDist(countnode)
            countnode=countnode+1;
        end
        SwathBoundIDs(countSwath)=countnode;
        NodeID=countnode;
    end

    %now impose minimum reach length
    %sort reaches by length
    [MinLength,ReachID ]=min(ReachLength);
    while MinLength(1)<MinReachLen
        %merge reaches until all are larger than the minimum length
        %find the smaller reach. merge it with the surroundings until it is
        %larger than the minimum length. Once that reach is large enough, than
        %sort again and proceed with merging    
        if ReachID>1
        %the reach is not the first in the series, so check if its
        %sinuisity is closer to the up or downstream's
            if ReachID+2<=length(ReachBoundaries)
                %then there is an upstream reach, otherwhise, ReachID is
                %the last reach
                AveSinCurr=mean(S(ReachBoundaries(ReachID):ReachBoundaries(ReachID+1)));
                AveSdownstr=mean(S(ReachBoundaries(ReachID+1):ReachBoundaries(ReachID+2)));
                AveSupstr=mean(S(ReachBoundaries(ReachID-1):ReachBoundaries(ReachID)));
                if abs(AveSdownstr-AveSinCurr)<abs(AveSupstr-AveSinCurr)
                    %current is more similar to the downstream
                    ReachBoundaries=[ReachBoundaries(1:ReachID);ReachBoundaries(ReachID+2:length(ReachBoundaries))];
                    ReachLength(ReachID)=FlowDist(ReachBoundaries(ReachID+1))-FlowDist(ReachBoundaries(ReachID));
                    ReachLength=[ReachLength(1:ReachID);ReachLength(ReachID+2:length(ReachLength))];
                else
                    %current is more similar to the upstream
                    ReachBoundaries=[ReachBoundaries(1:ReachID-1);ReachBoundaries(ReachID+1:length(ReachBoundaries))];
                    ReachLength=[ReachLength(1:ReachID-1);ReachLength(ReachID+1:length(ReachLength))];
                    ReachID=ReachID-1;
                    ReachLength(ReachID)=FlowDist(ReachBoundaries(ReachID+1))-FlowDist(ReachBoundaries(ReachID));                                                            
                end
            else
                %ReachID is the last reach, so merge with upstream
                ReachBoundaries=[ReachBoundaries(1:ReachID-1);ReachBoundaries(ReachID+1)];
                ReachID=ReachID-1;
                ReachLength=ReachLength(1:ReachID);
                ReachLength(ReachID)=FlowDist(ReachBoundaries(ReachID+1))-FlowDist(ReachBoundaries(ReachID));
            end
        else
            %Since it is the first reach, merge first and second reaches
            ReachBoundaries=[ReachBoundaries(1);ReachBoundaries(3:length(ReachBoundaries))];
            ReachLength=ReachLength(2:length(ReachLength));
            ReachLength(1)=FlowDist(ReachBoundaries(2))-FlowDist(ReachBoundaries(1));
        end
        [MinLength,ReachID ]=min(ReachLength);
    end
    %add boundaries at the Swath edge
    NodeID=1;
    for countSwath=1:length(SwathBoundIDs)
        count=NodeID;
        while count<length(ReachBoundaries)&&ReachBoundaries(count)<SwathBoundIDs(countSwath)
            count=count+1;
        end
        if ReachBoundaries(count)==SwathBoundIDs(countSwath)
            %nothing to do. the SwathBoundary is already a reach boundary
        else
            %add SwathBoundary to the reach boundary
            %ReachBoundaries=[ReachBoundaries(1:count-1);SwathBoundIDs(countSwath);ReachBoundaries(count:length(ReachBoundaries))];
            ReachBoundaries(count-1)=SwathBoundIDs(countSwath);
            %update reach lenght
            ReachLength(count-1)=FlowDist(ReachBoundaries(count))-FlowDist(ReachBoundaries(count-1));
            if count<length(ReachLength)
                %update the downstream reach length
                ReachLength(count)=FlowDist(ReachBoundaries(count+1))-FlowDist(ReachBoundaries(count));
                if ReachLength(count)<MinLength
                    %merge with downstream reach
                    ReachBoundaries=[ReachBoundaries(1:count);ReachBoundaries(count+2:length(ReachBoundaries))];
                    ReachLength(count)=FlowDist(ReachBoundaries(count+1))-FlowDist(ReachBoundaries(count));
                    ReachLength=[ReachLength(1:count);ReachLength(count+2:length(ReachLength))];
                end
            end
        end
        NodeID=count;
    end
    if Makeplots ==1
        Dbound=D(ReachBoundaries)/1000;
        %%1st derivative of sinuosity
        Dsdr=zeros(size(S));
        Dsdr(1)=(S(2)-S(1))/(D(2)-D(1));
        for i=2:n-1
            Dsdr(i)=(S(i+1)-S(i-1))/(D(i+1)-D(i-1));
        end
        Dsdr(n)=(S(n)-S(n-1))/(D(n)-D(n-1));

        figure(1)
        plot(X./1000,Y./1000)
        hold on
        plot(X(ReachBoundaries)./1000,Y(ReachBoundaries)./1000,'o',...
                           'MarkerEdgeColor','b',...
                           'MarkerFaceColor','b',...
                           'MarkerSize',6) 
        axis equal
        legend('River centerline','Reach boundaries');
        xlabel('Easting (km)')
        ylabel('Northing (km)')

        figure(2)
        plot(D./1000,S)
        hold on
        plot(Dbound,S(ReachBoundaries),'o',...
                           'MarkerEdgeColor','b',...
                           'MarkerFaceColor','b',...
                           'MarkerSize',6)

        ylabel('Sinuosity') % left y-axis
        xlabel('Flow distance in km')
        legend('Sinuosity','Reach boundaries', 'Location', 'northeast')
    end
end