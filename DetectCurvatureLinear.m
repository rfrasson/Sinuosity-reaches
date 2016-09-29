function [ReachBoundaries,Curvature,ReachNum,ReachLength] = DetectCurvatureLinear(x,y,ImposeMinLength,MinReachLen,numbregresspts,tcrit)
%This function searches for inflection points on the curve described by the
%x and y input variables.
%
%List of Inputs:
%x                : flow distance in km.
%y                : Water surface elevation (m) or sinuosity depending on the method
%ImposeMinLength  : 1 = imposes a minimum reach length, removing up and
%                   downstream boundaries if reach length < minimum. Not
%                   suitable for sinuosity calculations
%MinReachLen      : Minimum reach length in km
%numbregresspts   : number of regression points used in the calculation of
%                   dy/dx, which are derived from the slope of regression
%                   lines fitted to the up- and downstream 
%tcrit            : number of standard deviations used in the computation
%                   confidence interfals for dy/dx

%List of outputs:
%ReachBoundaries  : Location of the reach boundaries expressed in terms of
%                   the indices of the vector x
%Curvature        : Estimated value of d2y/dx2 at each node
%ReachNum         : Number of reaches
%ReachLength      : Length of each reach.

    first=1;
    current=numbregresspts;
    last=current+numbregresspts-1;
    SlopeUp=zeros(2,1); %confidence interval for the slope
    SlopeDown=zeros(2,1);
    Curvature=zeros(length(x),1);
    incrementstep=ceil(numbregresspts/2);
    MaxNumbAttempts=5;
    while last<=length(x)&&current<=length(x)
        overlap=1;%1 means the down- and upstream slope CI overlap 0 means it doesn't. Initially set to 1 to enter the loop at least once
        increment=0;%
        Xup=x(first:current);
        Yup=y(first:current); 
        Xdown=x(current:last);
        Ydown=y(current:last);
        NumbAttemps=0;
        while overlap==1&&NumbAttemps<=MaxNumbAttempts     
            PolCoeff=polyfit(Xup,Yup,1);
            %find the predicted points
            PredY = PolCoeff(1)*Xup+PolCoeff(2);
            s=sqrt(1/(current-first-1)*sum((Yup-PredY).^2)); %regression standard error. Page 566 from basic practice of statistics
            meanXup=mean(Xup);%average X for the interval. Saved in a variable as it is used in the calculation of the curvature
            SEb=s/sqrt(sum((Xup-meanXup).^2));
            SlopeUp(1)=PolCoeff(1)-tcrit*SEb;
            SlopeUp(2)=PolCoeff(1)+tcrit*SEb;
            PolCoeff=polyfit(Xdown,Ydown,1);
            %find the predicted points
            PredY = PolCoeff(1)*Xdown+PolCoeff(2);
            s=sqrt(1/(current-first-1)*sum((Ydown-PredY).^2)); %regression standard error. Page 566 from basic practice of statistics
            meanXdown=mean(Xdown);%average X for the interval. Saved in a variable as it is used in the calculation of the curvature
            SEb=s/sqrt(sum((Xdown-meanXdown).^2));
            SlopeDown(1)=PolCoeff(1)-tcrit*SEb;
            SlopeDown(2)=PolCoeff(1)+tcrit*SEb;
            if (SlopeDown(1)>=SlopeUp(1)&&SlopeDown(1)<=SlopeUp(2)||SlopeDown(2)>=SlopeUp(1)&&SlopeDown(2)<=SlopeUp(2))||...
                    (SlopeUp(1)>=SlopeDown(1)&&SlopeUp(1)<=SlopeDown(2)||SlopeUp(2)>=SlopeDown(1)&&SlopeUp(2)<=SlopeDown(2))
                %two intervals overlap, increase the size of the interval
                increment=increment+incrementstep;
                if first-increment>0
                    %add a point downstream
                    Xup=x(first-increment:current);
                    Yup=y(first-increment:current);
                    if last+increment<=length(x);
                        Xdown=x(current:last+increment);
                        Ydown=y(current:last+increment);
                    end
                else
                    if last+increment<=length(x);
                        Xdown=x(current:last+increment);
                        Ydown=y(current:last+increment);
                    else
                        %cannot increase X and Y sizes anymore best I can
                        %do is use the slopes to compute the curvature.
                        overlap=0;                        
                    end
                end
            else
                %the two slopes are different, compute curvature
                overlap=0;                
            end
            NumbAttemps=NumbAttemps+1; %tracks  how many times we attempted to find the curvature. If the maximum number is reached, we just use what we have.
        end
        Curvature(current)=(mean(SlopeDown)-mean(SlopeUp))/(meanXdown-meanXup);
        first=first+1;
        current=current+1;
        last=last+1;
    end
    %since I cannot compute the curvature at the edges of the river, assume
    %they are the same as the first computed curvature
    Curvature(1:numbregresspts-1)=Curvature(numbregresspts);    
    Curvature(length(Curvature)-numbregresspts+2:length(Curvature))=Curvature(length(Curvature)-numbregresspts+1);
    
    ReachBoundaries=zeros(length(Curvature),1);
    initial=1;
    ReachBoundaries(1)=initial;
    ReachNum=1;

    %Break up reaches (at this point no restriction on length)
    ReachBoundaries(1)=1;
    for count=2:length(Curvature)-1
        if Curvature(count)*Curvature(count-1)<0
            %sign reversal
            ReachNum=ReachNum+1;
            ReachBoundaries(ReachNum)=count;%begining of reach           
        end
    end
    ReachNum=ReachNum+1;
    ReachBoundaries(ReachNum)=length(x); %last reach ends on the last node.
    ReachBoundaries=ReachBoundaries(1:ReachNum);
    ReachNum=ReachNum-1; %now this is the number of reaches.
    %impose minimum reach leangh
    if ImposeMinLength==1
        count=2;
       while count<=length(ReachBoundaries)
           if x(ReachBoundaries(count)-1)-x(ReachBoundaries(count-1))<MinReachLen
               %merge reaches
               if count==2
                   %first reach, all that needs to be done is merge 1st and
                   %second
                   Curvature(1:ReachBoundaries(count)-1)=Curvature(ReachBoundaries(count));
                   ReachBoundaries=[ReachBoundaries(1:count-1);ReachBoundaries(count+1:length(ReachBoundaries))];
                   ReachNum=ReachNum-1;
                                      
               else
                   if count<length(ReachBoundaries)
                       Curvature(ReachBoundaries(count-1):ReachBoundaries(count)-1)=Curvature(ReachBoundaries(count-1)-1);
                       ReachBoundaries=[ReachBoundaries(1:count-2);ReachBoundaries(count+1:length(ReachBoundaries))];
                       ReachNum=ReachNum-2;                       
                   else
                       %merging last with the one before last
                       Curvature(ReachBoundaries(count-1):ReachBoundaries(count)-1)=Curvature(ReachBoundaries(count-1)-1);
                       ReachBoundaries=ReachBoundaries(1:count-1);
                       ReachNum=ReachNum-1;
                   end
               end 
           else
               count=count+1;
           end
       end
    end
    %calculate lengths
    ReachLength=zeros(ReachNum,1);
    for count=2:ReachNum+1
        ReachLength(count-1)=x(ReachBoundaries(count))-x(ReachBoundaries(count-1));
    end
end

