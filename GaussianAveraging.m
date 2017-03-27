function [ yave, wave ] = GaussianAveraging( x,y,w,WindowSize,sigma)
%This function performs smoothing using Gaussian weigths
%List of inputs:
%x          : Flow distances associated with the nodes
%y          : Water surface elevations
%w          : River widths at the nodes
%WindowSize : Size of the Window in the same units as x
%sigma      : Standard deviation of the normal distribution used to assign
%             the averaging weights

%List of outputs:
%yave       : smoothed elevations
%wave       : smoothed widths
    yave=zeros(size(y));
    wave=zeros(size(y));
    for count=1:length(x)
        first=count;
        while x(first)> x(count)-WindowSize/2 && first >1
            first = first -1;
        end
        last = count;
        while x(last)< x(count)+WindowSize/2 && last <length(x)
            last = last +1;
        end
        sumweight=0;
        accy=0;
        accw=0;
        for count2=first:last
            %average
            weight=normpdf(x(count2),x(count),sigma);
            if ~isnan(y(count2)) %if y at that point is missing, don't use it in the averaging
                sumweight=sumweight+weight;
                accy=accy+weight*y(count2);
                accw=accw+weight*w(count2);
            end
        end
        yave(count)=accy/sumweight; %smoothed y, not the reach average
        wave(count)=accw/sumweight;
    end
end

