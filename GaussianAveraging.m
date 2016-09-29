function [ y, w ] = GaussianAveraging( x,y,w,WindowSize,sigma)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
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
            sumweight=sumweight+weight;
            accy=accy+weight*y(count2);
            accw=accw+weight*w(count2);
        end
        yave(count)=accy/sumweight; %smoothed y, not the reach average
        wave(count)=accw/sumweight;
    end
    y=yave; %smoothed y replaces the riverobs height
    w=wave;
end

