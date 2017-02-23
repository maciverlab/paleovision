function [interpVisualRange,C0RangeSymExtendedNew] =interpolateContrastRange(C0Range,visualRangeSolns)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interpolate contrast to include the 0 point at which the visual range will be 0
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [C0RangeSym,indSym]=sort([C0Range,-C0Range(C0Range<0 & C0Range>-0.3),-C0Range(C0Range>0& C0Range<0.3)]);
    indp=find(C0Range<0 &C0Range>-0.3); indm=find(C0Range>0 &C0Range<0.3);
    [C0RangeSymExt,indExt]=sort([C0RangeSym,0]);
   
    C0RangeSymExtendedNew=linspace(min(C0Range),max(C0Range),101);    
    interpVisualRange=zeros(length(C0RangeSymExtendedNew),...
        size(visualRangeSolns,2),...
        size(visualRangeSolns,3));
    
    for k=1:size(visualRangeSolns,3)
        for j=1:size(visualRangeSolns,2)
            temp=visualRangeSolns(:,j,k);
            temp=[temp;
                visualRangeSolns(indp(1):indp(end),j,k);
                visualRangeSolns(indm(1):indm(end),j,k)]; temp=temp(indSym); 
            temp=[temp;0]; temp=temp(indExt);
            
            interpRange=interp1(C0RangeSymExt,temp,C0RangeSymExtendedNew,'pchip');
            interpVisualRange(:,j,k)=interpRange;
        end
    end
end
