function Aquatic_firingThresh
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of visual range based on optical properties
%%
%% Title                : A massive increase in visual range preceded the origin of terrestrial vertebrates
%% Authors              : Ugurcan Mugan, Malcolm A. MacIver
%% Authors' Affiliation : Northwestern University
%% This work is licensed under the Creative Commons Attribution 4.0 International License. 
%% To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
%% January 2017
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global BIGEYEROOT
%% Variable Initialization
    run Parameters.m
    load('Parameters.mat')
    conditions={'Daylight','Moonlight','Starlight'};
    viewing={'Upward','Horizontal'};
    DrangeAquatic=linspace(minpupil,maxpupil,30); 
    depth=8; epsilon=5e-4;

%% Given condition, viewing solve for range for all pupil diameter 
    visualRangeAquatic=zeros(length(Drange.Aquatic),length(conditions),length(viewing));
    rInit.Upward=[5,2,1];
    rInit.Horizontal=[3,1,.3];
    for c=1:length(conditions);
        % get values that are only condition dependent and viewing
        % independent
        a=aAquatic.(conditions{c}); 
        b=bAquatic.(conditions{c});
        C0=C0Aquatic.(conditions{c});
        pAbsorb=pAbsorbAquatic.(conditions{c});
        for j=1:length(viewing);
            % get values that are both condition and viewing dependent
            K=KAquatic.(conditions{c}).(viewing{j});
            L=LAquatic.(conditions{c}).(viewing{j});
            
            r=rInit.(viewing{j})(c);
            %for all pupil values find corresponding range that solves the
            %firing threshold
            for i=1:length(Drange.Aquatic)
                D=DrangeAquatic(i);
                delta=10^(floor(log10(r))-4);
                
                P=10;
                while abs(P-1)>epsilon
                    P=firingThresh(depth,lambda,...
                        pAbsorb,a,b,K,L,...
                        r,D,XAquatic,DtAquatic,qAquatic,d,k,len,T,M,R,C0);
                    if P>1
                        r=r-delta;
                    else
                        r=r+delta;
                    end
                    clc;
                    fprintf('pupil diameter: %f\nrange: %f\nerror: %f\n',...
                        D,r,abs(P-1));
                end
                visualRangeAquatic(i,c,j)=r;
            end
        end
    end
    
%% save results
save([BIGEYEROOT 'fig04_visualrange/aquatic_model/meteoAquatic.mat'],'visualRangeAquatic','DrangeAquatic');

function  solution=firingThresh(depth,lambda,photoreceptorAbsorption,a,b,KAll,LAll,r,D,X,Dt,q,d,k,len,T,M,R,C0)
    lambda1=lambda(1); lambda2=lambda(end);
    K=KAll(:,depth);
    L=LAll(:,depth);
   
    alphaInterp=@(l) interp1(lambda,photoreceptorAbsorption,l,'pchip');
    aInterp=@(l) interp1(lambda,a,l,'pchip');
    bInterp=@(l) interp1(lambda,b,l,'pchip');
    LInterp=@(l) interp1(lambda,L,l,'pchip');
    KInterp=@(l) interp1(lambda,K,l,'pchip');
    
    Nfalse=((T*M*D)/(2*r*d))^2*X*Dt; %Supplementary Appendix, pg 3 for explanation.
    %Supplementary Appendix, Eq 7,8 pg. 5
    %Integrand
    RhFunc=@(l) LInterp(l).*l.*(1-exp(-k*alphaInterp(l)*len));
    RoFunc=@(l) LInterp(l).*l.*(1-exp(-k*alphaInterp(l)*len)).*(1+(C0*exp((KInterp(l)-aInterp(l)-bInterp(l))*r)));
    %Integral
    Rh=integral(RhFunc,lambda1,lambda2);
    Ro=integral(RoFunc,lambda1,lambda2);
    Nh=((pi/4)^2)*(D^2)*((T/r)^2)*q*Dt*Rh;
    No=((pi/4)^2)*(D^2)*((T/r)^2)*q*Dt*Ro;
    %Supplementary Appendix, Eq. 9 pg. 5
    solution=(R*sqrt(No+Nh+2*Nfalse))/(abs(No-Nh));    