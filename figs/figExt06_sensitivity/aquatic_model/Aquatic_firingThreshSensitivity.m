function Aquatic_firingThreshSensitivity
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of visual range based on optical properties with different water types
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global BIGEYEROOT
%% Intialize variables   
    run ParametersSensitivity.m
    load('ParametersSensitivity.mat');
    
    conditions={'HighTurbidity','Clear','AbsDom','ScatDom'};
    viewing={'Upward','Horizontal'};
    waterDepth=8; 
    epsilon=5e-4;

    %% For each water condition and viewing direction calculate visual range based on optical properties
    DrangeAquatic=linspace(minpupil,maxpupil,30);   
    visualRangeAquaticSensitivity=zeros(length(DrangeAquatic),length(conditions),length(viewing));
    rAll=[3,10,5,2.3;1.5,10,3,1.5];
    for c=1:length(conditions)
        aValue=a.(conditions{c}); 
        bValue=b.(conditions{c});
        pAbsorbValue=pAbsorb.(conditions{c});
        for v=1:length(viewing)
            fprintf('condition: %s\nviewing: %s\n',conditions{c}, viewing{v});
            KValue=K.(conditions{c}).(viewing{v});
            LValue=L.(conditions{c}).(viewing{v});           
            
            r=rAll(v,c);            
            for j=1:length(DrangeAquatic)
                    D=DrangeAquatic(j);
                    delta=10^(floor(log10(r))-5);
                    
                    P=10;
                    while abs(P-1)>epsilon
                        P=firingThresh(waterDepth,lambda,...
                             pAbsorbValue,aValue,bValue,KValue,LValue,...
                             r,D,X,Dt,q,d,k,len,T,M,R);
                         if P>1
                             r=r-delta;
                         else
                             r=r+delta;
                         end
                         %clc;
                         %fprintf('condition: %s\nviewing: %s\npupil diameter: %f\nrange: %f\nerror: %f\n',...
                         %    conditions{c}, viewing{v}, D, r,abs(P-1) );
                    end
                    visualRangeAquaticSensitivity(j,c,v)=r;
            end
        end
    end

save([BIGEYEROOT 'figExt06_sensitivity/aquatic_model/meteoSensitivityAquatic.mat'],'visualRangeAquaticSensitivity','DrangeAquatic');
     
function  solution=firingThresh(depth,lambda,photoreceptorAbsorption,aAll,bAll,KAll,LAll,r,A,X,Dt,q,d,k,len,T,M,R) 
    L=LAll(:,depth);
    [~,ind]=max(L);
    a=aAll(:); a=a(ind);
    K=KAll(:,depth); K=K(ind);
    b=bAll(:); b=b(ind);
    
    Nfalse=((T*M*A)/(2*r*d))^2*X*Dt;
    %Supplementary Appendix, pg 3 for explanation.
    %Supplementary Appendix, Eq 7,8 pg. 5
    
    Rh=L.*lambda.*(1-exp(-k*photoreceptorAbsorption*len));
    Ro=Rh.*(1-exp((K-a-b)*r));
    
    %Integral
    Bh=trapz(lambda,Rh);
    Bo=trapz(lambda,Ro);
    if isinf(Bo)
        Bo=0;
    end
    
    Nh=((pi/4)^2)*(A^2)*((T/r)^2)*q*Dt*Bh;
    No=((pi/4)^2)*(A^2)*((T/r)^2)*q*Dt*Bo;
    %Supplementary Appendix, Eq. 9 pg. 5
    solution=(R*sqrt(No+Nh+2*Nfalse))/(abs(No-Nh));
