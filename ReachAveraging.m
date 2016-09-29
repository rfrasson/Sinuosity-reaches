function [Reach,RiverData,Metadata,ReachTrue,Nodes,NodesTrue]=ReachAveraging(ReachBoundaries, RiverObs,True,RefRiverObs,RefTrue,ReachLength,Day,SaveResults,SmoothData,VariableSmoothingWindow,OutputPath, OutFileName,MakePlots)
%this function calculates reach-averaged quantities. 
%List of inputs
%ReachBoundaries     : Node ids of reach boundaries
%RiverObs nodes      : Structure containing RiverObs nodes of the day the
%                      user intends to calculate reach-averaged quantities
%                      (see explanation of the structure below)

%True                : True values of water surface elevation and width for
%                      the evaluation of the reach-averaged quantities
%                      Same format as RiverObs nodes

%RefRiverObs nodes   : Structure containing nodes from the day used to 
%                      define the centerline. Cross-sectional area changes 
%                      are computed with reference to this day(same format as RiverObs)

%RefTrue             : true values during the reference overpass. Same
%                      format as RiverObs
%

%ReachLength         : length of each reach in km

%SmoothData          : Smooths data using a moving average Gaussian window.
%                      Further configurations follows

%VariableSmoothingWindow  : %0 will not use gaussian window
%                           1 will use variable size window 
%                           2 will use the standard window size of 10km

%OutputPath          :Path to the output directory (including final /)

%OutFileName         :output file name

%MakePlots           : 1 makes plots

%list of outputs:
%Reach               : Structure holding the reach-averaged products
%                      (simpler format)
%RiverData           : Structure holding the Reach-averaged products as
%                      explained in River Data Products drocument
%Metadata            : description of the fields inside RiverData
%ReachTrue           : True values of reach averaged products
%Nodes               : Structure holding node data products
%NodesTrue           : True values at the nodes



%Clarification on the format of RiverObs nodes.
%The variable passed as RiverObs must be a structure containing the
%following elements (as vectors with dimenstion eaqual to Number of nodes,1:
%x          : flow distance measured in meters
%H          : Node elevation in meters
%Easting    : Projected coordinate of the node in m
%Northing   : Projected coordinate of hte node in m
%Lat        : Latitude of the node in decimal degrees
%Lon        : Longitude of the node in decimal degrees
%xtrack     : cross-track position of the node





    %% Reach averaging and discharge calculations
    WindowSize=10;
    sigma=2;
    calculateQ=1; %1 estimates Q for each reach


    x=RiverObs.x/1000;
    y=RiverObs.H;
    w=RiverObs.W;
    xtrue=True.x/1000;
    ytrue=True.H;
    wtrue=True.W;
    trueA=True.Area;
    refw=RefRiverObs.W;
    refy=RefRiverObs.H;        
    refTrueA=RefTrue.Area;
    refTrueW=RefTrue.W;
    refTrueH=RefTrue.H;
    
    %checks if there are more than 1 crossection at the same x location
    counter=1;
    for count=1:length(xtrue)-1
        if xtrue(count+1)==xtrue(count)
            xtrue(count+1)=xtrue(count+1)+0.0001;
            counter=counter+1;
        end
    end
    ytrueInt=interp1(xtrue,ytrue,x);
    wtrueInt=interp1(xtrue,wtrue,x);
    trueAInt=interp1(xtrue,trueA,x);
    TrueQInt=interp1(xtrue,True.Q,x);
    refTrueAInt=interp1(xtrue,refTrueA,x);
    refwtrueInt=interp1(xtrue,refTrueW,x);
    refytrueInt=interp1(xtrue,refTrueH,x); 

    if x(1)<xtrue(1)
        ytrueInt(1)=(ytrue(1));
    end
    ytrueInt=ytrueInt(1:length(x));
    len=length(x);
    if isnan(ytrueInt(len))
        coef=polyfit(x(len-2:len-1),ytrueInt(len-2:len-1),1);
        ytrueInt(len)=coef(1)*x(len)+coef(2);
    end
    
    if x(1)<xtrue(1)
        wtrueInt(1)=(wtrue(1));
    end
    
    wtrueInt=wtrueInt(1:length(x));

    if isnan(wtrueInt(len))
        coef=polyfit(x(len-2:len-1),wtrueInt(len-2:len-1),1);
        wtrueInt(len)=coef(1)*x(len)+coef(2);
    end

    trueAInt(1)=(trueA(1));
    trueAInt=trueAInt(1:length(x));

     if isnan(trueAInt(len))
        coef=polyfit(x(len-2:len-1),trueAInt(len-2:len-1),1);
        trueAInt(len)=coef(1)*x(len)+coef(2);
    end

    TrueQInt(1)=(True.Q(1));
    TrueQInt=TrueQInt(1:length(x));

    if isnan(TrueQInt(len))
        coef=polyfit(x(len-2:len-1),TrueQInt(len-2:len-1),1);
        TrueQInt(len)=coef(1)*x(len)+coef(2);
    end    

    refTrueAInt(1)=(refTrueA(1));
    refTrueAInt=refTrueAInt(1:length(x));

    if isnan(refTrueAInt(len))
        coef=polyfit(x(len-2:len-1),refTrueAInt(len-2:len-1),1);
        refTrueAInt(len)=coef(1)*x(len)+coef(2);
    end  

    refwtrueInt(1)=(refTrueW(1));
    refwtrueInt=refwtrueInt(1:length(x));

    if isnan(refwtrueInt(len))
        coef=polyfit(x(len-2:len-1),refwtrueInt(len-2:len-1),1);
        refwtrueInt(len)=coef(1)*x(len)+coef(2);
    end 

    refytrueInt(1)=(refTrueH(1));
    refytrueInt=refytrueInt(1:length(x));

    if isnan(refytrueInt(len))
        coef=polyfit(x(len-2:len-1),refytrueInt(len-2:len-1),1);
        refytrueInt(len)=coef(1)*x(len)+coef(2);
    end

    Nodes.x=x;
    Nodes.y=y;
    Nodes.w=w;
    NodesTrue.x=x;
    NodesTrue.y=ytrueInt;
    NodesTrue.w=wtrueInt;
    
    
    yave=zeros(size(y));
    wave=zeros(size(y));
    yaveRef=zeros(size(y));
    waveRef=zeros(size(y));
    slope=zeros(size(y));
    trueslope=zeros(size(ytrue));
    counter=1;
    countery=0;
    counterSlope=1;
    ReachNum=length(ReachLength);
    if SmoothData==1
        originaly= y;
        if VariableSmoothingWindow==1
            for countReaches=1:ReachNum
                %Smooth y and w
                %for each reach, I'm using ReachLength as the size of the
                %moving average window. That does NOT limit the averaging
                %window to reach boundaries, only how far the window extends
                r1=ReachBoundaries(countReaches); %index of the first node belonging to the reach that is being smoothed
                r2=ReachBoundaries(countReaches+1);%index of the last node belonging to the reach that is being smoothed
                if countReaches>1
                    re1=ReachBoundaries(countReaches-1); %begining of the upstream reach. those are used in the averaging process, but not saved.
                else
                    re1=r1;
                end
                if countReaches<ReachNum
                    re2=ReachBoundaries(countReaches+2); %end of the downstream reach. those are used in the averaging process, but not saved.
                else
                    re2=r2;
                end
                WS = ReachLength(countReaches);%window size equal to the current reach's length
%                 if WS>10
%                     WS=10;
%                 end
                NewSigma=WS/5; %for a 10km window, we want a 2km standard deviation for the Gaussian averaging window. Keep mininum at 1km, which is about 3 nodes
                if NewSigma < 1
                    NewSigma =1;
                end
                [ynew, wnew] = GaussianAveraging(x(re1:re2),y(re1:re2),w(re1:re2),WS,NewSigma);            
                yave(r1:r2)=ynew(r1-re1+1:r2-re1+1);%saves the averaged ys for the current reach only
                wave(r1:r2)=wnew(r1-re1+1:r2-re1+1);
                [ynew, wnew] = GaussianAveraging(x(re1:re2),refy(re1:re2),refw(re1:re2),WS,NewSigma);
                yaveRef(r1:r2)=ynew(r1-re1+1:r2-re1+1);
                waveRef(r1:r2)=wnew(r1-re1+1:r2-re1+1);
            end
        else
            if VariableSmoothingWindow==2
                %10km window            
                [yave, wave] = GaussianAveraging(x, y,w,WindowSize,sigma);
                [yaveRef, waveRef] = GaussianAveraging(x, refy,refw,WindowSize,sigma);
            else
                %no averaging!!
                yave=y;
                wave=w;
            end
        end
    else
        yave=y;
        wave=w;
        yaveRef=refy;
        waveRef=refw;
    end
    
   

    %RMSE for height
    sqdiff=(y-ytrueInt).^2;
    rmseReachHeight=sqrt(mean(sqdiff));
    Output=sprintf('RMSE for water surface elevation = %0.5g m',rmseReachHeight);
    display(Output)

    %calculate slope for each node
    slope(1)=-(yave(2)-yave(1))/(x(2)-x(1));
    trueslope(1)=-(ytrueInt(2)-ytrueInt(1))/(x(2)-x(1));
    for count=2:length(x)
        slope(count)=-(yave(count)-yave(count-1))/(x(count)-x(count-1));
        trueslope(count)=-(ytrueInt(count)-ytrueInt(count-1))/(x(count)-x(count-1));
    end


    %xreach=FlowDist(ReachBoundaries);
    xreach=x(ReachBoundaries);
    sizexreach=size(xreach);
    sizexreach(1)=sizexreach(1)-1;
    xmid=zeros(sizexreach);
    trueReachslope=zeros(sizexreach);
    Reachslope=zeros(sizexreach);
    Reachy=zeros(sizexreach);
    Reachw=zeros(sizexreach);
    ReachA=zeros(sizexreach);
    ReachL=zeros(sizexreach);
    ReachyRef=zeros(sizexreach);
    ReachwRef=zeros(sizexreach);
    TrueReachw=zeros(sizexreach);
    TrueReachy=zeros(sizexreach);
    TrueReachA=zeros(sizexreach);
    TrueReachQ=zeros(sizexreach);
    TrueReachN=zeros(sizexreach);
    RefTrueReachy=zeros(sizexreach);
    RefTrueReachw=zeros(sizexreach);
    TrueReachA0=zeros(sizexreach);

    first=1;
    for count=1:length(xreach)-1
        xmid(count)=(xreach(count)+xreach(count+1))/2;
        last = first;
        while x(last)< xreach(count+1) && last <length(x)
            last = last +1;
        end
        last=last-1;

        trueReachslope(count)=-(ytrueInt(last)-ytrueInt(first))/(x(last)-x(first));
        TrueReachy(count)=mean(ytrueInt(first:last));
        TrueReachw(count)=mean(wtrueInt(first:last));
        TrueReachA(count)=mean(trueAInt(first:last));
        TrueReachA0(count)=mean(refTrueAInt(first:last));
        RefTrueReachy(count)=mean(refytrueInt(first:last)); %reference day
        RefTrueReachw(count)=mean(refwtrueInt(first:last));
        if calculateQ==1
            TrueReachQ(count)=mean(TrueQInt(first:last));
            TrueReachN(count)=TrueReachA(count)^(5/3)*TrueReachw(count)^(-2/3)*(trueReachslope(count)/1000)^(1/2)/TrueReachQ(count);
        end

        Reachslope(count)=-(yave(last)-yave(first))/(x(last)-x(first));
        Reachy(count)=mean(yave(first:last));
        Reachw(count)=mean(wave(first:last));
        ReachyRef(count)=mean(yaveRef(first:last));
        ReachwRef(count)=mean(waveRef(first:last));
        ReachL(count)=x(last)-x(first);
        first=last+1;
    end

    diffslope=(100*Reachslope-100*trueReachslope);
    rmseReachSlope=sqrt(mean(diffslope.^2));
    Output=sprintf('RMSE for reach slope = %0.5g cm/km',rmseReachSlope);
    display(Output)
    if MakePlots==1
        figure(3)    
        plot(x,y,'--b')
        hold on
        plot(x,yave,'-r')
        plot(xtrue,ytrue,'--k')
        title('Water surface profile')
        xlabel('Flow distance (km)')
        ylabel('Elevation (m)');
        legend('Simulated nodes', 'Smoothed surface', 'True WSP')
        
        figure
        plot(x,slope,'-b')
        hold on
        plot(x,trueslope,'--r');
        title('Water surface slope')
        xlabel('Flow distance (km)')
        ylabel('Slope (m/km)');
        legend('Node level slope', 'True slope')
        
        figure
        plot(xmid,Reachslope,'-b')
        hold on
        plot(xmid,trueReachslope,'--r');
        ylabel('Reach averaged slope m/km')
        xlabel('Flow distance')
    end
    


    if calculateQ==1
        ReachA=TrueReachA0+(Reachy-ReachyRef).*(ReachwRef+Reachw)/2;
        %ReachA=TrueReachA0+(Reachy-RefTrueReachy).*(RefTrueReachw+Reachw)/2;
        ReachQ=TrueReachN.^(-1).*ReachA.^(5/3).*Reachw.^(-2/3).*(Reachslope/1000).^0.5;
        diffA=(ReachA-TrueReachA);
        diffQ=(ReachQ-TrueReachQ);
        diffRQ=(100*(ReachQ-TrueReachQ)./TrueReachQ);
        diffW=(Reachw-TrueReachw);
        diffY=(Reachy-TrueReachy);
        rmseReachA=sqrt(mean(diffA.^2));
        rmseReachQ=sqrt(mean(diffQ.^2));
        rmseReachW=sqrt(mean(diffW.^2));
        rmseReachY=sqrt(mean(diffY.^2));
        rmseReachRQ=sqrt(mean(diffRQ.^2));
        Output=sprintf('RMSE for reach elevation = %0.5g m',rmseReachY);
        display(Output)
        Output=sprintf('RMSE for reach width = %0.5g m',rmseReachW);
        display(Output)
        Output=sprintf('RMSE for reach Area = %0.5g m^2',rmseReachA);
        display(Output)
        Output=sprintf('RMSE for reach Discharge = %0.5g m^2/s',rmseReachQ);
        display(Output)
        Output=sprintf('RMSE for reach relative Discharge = %0.5g %',rmseReachRQ);
        display(Output)
        Reach.A=ReachA';
        Reach.Q=ReachQ';
        Reach.W=Reachw';
        Reach.H=Reachy';
        Reach.N=TrueReachN';
        Reach.S=Reachslope';
        Reach.L=ReachL';
        Reach.day=Day;
        Reach.ReachBoundaries=ReachBoundaries;
        Reach.CoordBoundLat=RiverObs.Lat(ReachBoundaries);
        Reach.CoordBoundLon=RiverObs.Lon(ReachBoundaries);
        %RiverData has all the information according to the River Data products
        %document version 2
        RiverData.ReachID=1:length(Reach.A); %in the future, ReachID will hold all the  
        RiverData.Length=Reach.L; %Pass-based reach length
        RiverData.LandsatLength=Reach.L; %need to implement this. Probably calculate from the GRWL centerline
        RiverData.Height=Reach.H;
        StandardErrorH=zeros(size(Reach.H));
        StandardErrorSlope=zeros(size(Reach.H));
        ReachInundArea=zeros(size(Reach.H));
        count=1;
        %compute reach-averaged height and slope uncertainties, inundated area,
        %inundated area uncertainty (waiting on this one)
        for countreach=1:length(Reach.A)
           xcoord=x(ReachBoundaries(countreach):ReachBoundaries(countreach+1))/1000;
           yun=RiverObs.H(ReachBoundaries(countreach):ReachBoundaries(countreach+1));
           y=yave(ReachBoundaries(countreach):ReachBoundaries(countreach+1));
           w=wave(ReachBoundaries(countreach):ReachBoundaries(countreach+1));
           PolCoeff=polyfit(xcoord,yun,1);
           PredY = PolCoeff(1)*xcoord+PolCoeff(2); 
           %standard error for the reach averaged height
           StandardErrorH(countreach)=std(PredY-yun);
           %standard error for the slope
           %find the predicted points
           
           %Standard error for the slope using the regression error
           %estimates (if that is the intended mode, uncomment the
           %following two lines)
           
           %s=sqrt(1/(length(x)-2)*sum((yun-PredY).^2));
           %StandardErrorSlope(countreach)=s/sqrt(sum((x-mean(x)).^2))*100;
           
           %Standard error computed using the idealized formula. 
           Sigmah=4.4;% cm from SWOT error budget
           L=Reach.L(countreach);
           W=Reach.W(countreach)/1000;
           StandardErrorSlope(countreach)=2*pi()*sqrt(1/12*Sigmah^2./(W*L.^3));
  
           %calculating the inundated area by making trapezoids with
           %successive width measurements
           area=0;
           for counter=2:length(w);
               area=area+(x(counter)-x(counter-1))*(w(counter)+w(counter-1))/2;
           end
           ReachInundArea(count)=area/1000; %dividing by 1000 as width is in meters and length in km and I want result in km2
           count=count+1;
        end
        RiverData.HeightUncertainty=StandardErrorH;
        RiverData.Width=Reach.W;
        RiverData.Slope=Reach.S*100; %multiplying by 100 to transform to cm/km
        RiverData.SlopeUncertainty=StandardErrorSlope;
        RiverData.ReachInundArea=ReachInundArea';
        RiverData.ReachInundAreaUncertainty=zeros(size(RiverData.ReachInundArea)); %place holder
        RiverData.ReachDischarge=Reach.Q;
        RiverData.ReachDischargeUncertainty=zeros(size(RiverData.ReachDischarge)); %place holder
        RiverData.ReachAveragedLowFlowXsection=TrueReachA0';
        RiverData.ManningN=Reach.N;
        RiverData.IslandFlag=zeros(size(Reach.N)); %place holder
        RiverData.PlanformClass=ones(size(Reach.N)); %place holder for Planform classification
        RiverData.Flags=ones(size(Reach.N)); %place holder
        RiverData.Corrections.WetTroposphere=ones(size(Reach.N)); %place holder
        RiverData.Corrections.Ionosphere=ones(size(Reach.N)); %place holder
        RiverData.Corrections.Roll=ones(size(Reach.N)); %place holder

        %at this point my river is a single channel river, so connectivity is
        %straightforward
        RiverData.UpstreamReachID=zeros(size(RiverData.ReachID));
        RiverData.DownstreamReachID=zeros(size(RiverData.ReachID));
        RiverData.UpstreamReachID(1)=0; %nothing upstream
        RiverData.DownstreamReachID(length(RiverData.ReachID))=0; %nothing downstream
        for count=2:length(RiverData.ReachID)-1
            RiverData.UpstreamReachID(count)=RiverData.ReachID(count-1);
            RiverData.DownstreamReachID(count)=RiverData.ReachID(count+1);
        end
        RiverData.DistancefromOutlet=zeros(size(Reach.N)); %place holder
        for countreach=1:length(Reach.A)
            RiverData.CenterCoordinates(countreach).Lat=mean(RefRiverObs.Lat(ReachBoundaries(countreach):ReachBoundaries(countreach+1)));
            RiverData.CenterCoordinates(countreach).Lon=mean(RefRiverObs.Lon(ReachBoundaries(countreach:countreach+1)));
            RiverData.CenterCoordinates(countreach).Northing=mean(RefRiverObs.Northing(ReachBoundaries(countreach):ReachBoundaries(countreach+1)));
            RiverData.CenterCoordinates(countreach).Easting=mean(RefRiverObs.Easting(ReachBoundaries(countreach):ReachBoundaries(countreach+1)));
            RiverData.CenterCoordinates(countreach).CrosstrackDistance=mean(RefRiverObs.xtrack(ReachBoundaries(countreach):ReachBoundaries(countreach+1)));
        end
        RiverData.LayoverFlag=zeros(size(Reach.N));
        Metadata={'ReacID holds an unique ID for the reach', 'Length is the reach length in km', 'LandsatLength reach length calculated from landsat','Height is the elevation in meters - figure out the coordinate system',...
            'HeightUncertainty Reach averaged height uncertainty in m', 'Slope is the reach averaged slope in cm/km', 'Slope Uncertainty is the reach averaged slope in cm/km', ...
            'ReachInundArea is the reach inundated area in km2', 'ReachDischarge is the reach averaged discharge in m3/s', 'ReachDischargeUncertainty is the reach averaged discharge in m3/s',...
            'ReachAveragedLowFlowXsection is the reach averaged cross-sectional area during the lowest observed flow', 'ManningN is the reach averaged roughness coefficient', 'Islandflag 0 - islands are present, 1 islands might be present, 2 no islands, 3 unkown',...
            'River Planform classification single chanel, simple multi-channel, braided or anastomosing, river with floodplain interaction', 'Flags- snow, ice, high precipitation impacting signal', 'Corrections three fields, one for wet troposphere, one for ionosphere, one for roll correction',...
            'UpstreamReachID - ID of the upstream reach or reaches (in case of multi-channel) that connect to the current reach', 'DownstreamReachID - ID of the downstream reach or reaches (in case of multi-channel) that connect to the current reach',...
            'Centercoordinates - Latitude, longitude, easting, northing, and crosstrack distance for the center of the reach'};

        ReachTrue.A=TrueReachA';
        ReachTrue.Q=TrueReachQ';
        ReachTrue.W=TrueReachw';
        ReachTrue.H=TrueReachy';
        ReachTrue.N=TrueReachN';
        ReachTrue.S=trueReachslope';
        
        if MakePlots==1
            figure
            plot(xmid,ReachTrue.Q,'-r')
            hold on
            plot(xmid,Reach.Q,'-b')
            ylabel('Reach discharge m3/s')
            xlabel('Flow distance')
        end

        outfilename=[OutputPath OutFileName];
        if SaveResults==1
            save(outfilename,'Reach','RiverData','Metadata', 'ReachTrue','Nodes','NodesTrue')
        end
    end
end