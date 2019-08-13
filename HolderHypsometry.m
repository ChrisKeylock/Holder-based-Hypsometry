function output=HolderHypsometry(Z,alpha,numbins,plotflag)

%This function is an implementation of the methodology described in

% C. J. Keylock, A. Singh, P. Passalacqua, and E. Foufoula-Georgiou. 2019.
% Hölder-conditioned Hypsometry: A Refinement To A Classical Approach For 
% The Characterization of Topography. Water Resource Research

%It takes as input two rectangular arrays Z (elevations) and alpha (Holder exponents) of the same size. Returns the
%hypsometry and the Holder-conditioned hypsometric values. (If the fraclab
%toolbox has been used to generate alpha then these will be square arrays
%of dimension 2^J where J is an integer).

%The number of bins to discretise may be specified in numbins. If this is
%not provided, an estimate of the correct number is used.

%Output is a structure containing:
%output.Ihyp - classic hypsometric integral 
%output.Ialpha - The equivalent to Ihyp for the alpha values rather than Z
%  values (the integral of the HECAS defined below)
%output.Ihypcondalpha - Holder conditioned hypsometry
%output.SigmaIhypcondalpha - the standard deviation of Ihypcondalpha
%output.ThetaIhypcondalpha - the slope to the Ihypcondalpha values as a
%  function of alpha
%output.Ialphacondz - Holder exponent catchment area scaling (HECAS)
%output.alphavals - For each bin, the minimum (colums 1 and 4), 
%  median (columns (2 and 5) and maximum (columns 3 and 6) alpha values 
%  based on the rescaled (columns 1 -3) and original (columns 4-6) values.
%output.Zvals - For each bin, the minimum (colums 1 and 4), 
%  median (columns (2 and 5) and maximum (columns 3 and 6) elevation values 
%  based on the rescaled (columns 1 -3) and original (columns 4-6) values.

%Depending on plotting options chosen additional fields appear in output
%corresponding to the graphical elements generated.

%plotflag is a vector containing from 1 to 3 values ranging from 1 to 3.
%These specify the figures you wish to generate. If this is not provided,
%no figures are generated.

%The types of figures are as follows:
%1 - Two panels: The upper gives the overall hypsometric curves for Z and
%    alpha. It resembles Fig. 3c in the manuscript.
%    The lower gives the Holder conditioned hypsometry and the HECAS
%    The horizontal lines are the hypsometric integral and the Holder
%    equivalent. It resembles Fig. 6 in the manuscript.
%2 - Two panels: The upper panel shows the hypsometric functions 
%    conditioned on alpha for those closest to the 5%, 25% 50%, 75%, 95% 
%    values. The lower panel shows the Holder exponent catchment area 
%    scaling functions conditioned on z for those closest to 5%, 25% 50%, 
%    75%, 95% values. It resembles Fig. 5 in the manuscript.
%3 - The joint probability distribution functions for elevation (y-axis) 
%    and alpha (x-axis). The left hand panel contours the probabilities. 
%    The right-hand panel contours the log_10 probabilities  

if nargin<4
    plotflag=0;
else
    if length(plotflag)>3
        plotflag=plotflag(1:3);
    end
    flag=0;
    for loop1=1:3
        if plotflag(loop1)<1 || plotflag(loop1)>4
            flag=1;
        end
    end
    if flag==1
        disp('Plotting options inappropriately formatted')
        return
    end
end

minalpha=min(alpha(:));
maxalpha=max(alpha(:));
minZ=min(Z(:));
maxZ=max(Z(:));

sizer=size(Z);
sizeralpha=size(alpha);
if sum(sizer-sizeralpha)~=0
    disp('The elevation and Holder exponent arrays are not the same size')
    return
end

if nargin<3 || numbins < 25 || numbins>(sizer(1)*sizer(2)/10);
    numbins=-1;
end

%Make unique values
Z=Z+randn(sizer(1),sizer(2)).*eps;
alpha=alpha+randn(sizer(1),sizer(2)).*eps;
origZ=Z;origa=alpha;
Z=(Z-minZ)./(maxZ-minZ);
alpha=(alpha-minalpha)./(maxalpha-minalpha);

%Integrals calculated on individual points (no binning)
n{1}=sort(Z(:));
n{2}=sort(alpha(:));
n{3}=[1:length(n{1})]./length(n{1});

Ihyp=trapz(n{1},n{3});
Ialpha=trapz(n{2},n{3});

%Conditional integrals - determine binning
if numbins==-1
    numbins=max([25 round((log2(sizer(1)*sizer(2))/2)^2)]);
end

quantile_listZ=round(linspace(min(Z(:)),max(Z(:)),numbins+1).*(sizer(1)*sizer(2)));
quantile_listalpha=round(linspace(min(alpha(:)),max(alpha(:)),numbins+1).*(sizer(1)*sizer(2)));

%convert Z and alpha into quantiles
qZ=reshape(tiedrank(Z(:)),sizer);
qa=reshape(tiedrank(alpha(:)),sizer);

%zfuncalpha and alphafuncz are stored in case one wishes to inspect
for loop1=1:numbins
    temp=find(qa>=quantile_listalpha(loop1) & qa<=quantile_listalpha(loop1+1));
    alphavals(loop1,1)=min(alpha(temp));
    alphavals(loop1,2)=median(alpha(temp));
    alphavals(loop1,3)=max(alpha(temp));
    alphavals(loop1,4)=min(origa(temp));
    alphavals(loop1,5)=median(origa(temp));
    alphavals(loop1,6)=max(origa(temp));
    zfuncalpha{loop1}=sort(Z(temp));
    zfuncalpha{loop1}(length(temp)+1)=1;
    nc=[1:length(temp)+1]./(length(temp)+1);
    Ihypcondalpha(loop1)=trapz(zfuncalpha{loop1},nc);

    temp=find(qZ>=quantile_listZ(loop1) & qZ<=quantile_listZ(loop1+1));
    Zvals(loop1,1)=min(Z(temp));
    Zvals(loop1,2)=median(Z(temp));
    Zvals(loop1,3)=max(Z(temp));
    Zvals(loop1,4)=min(origZ(temp));
    Zvals(loop1,5)=median(origZ(temp));
    Zvals(loop1,6)=max(origZ(temp));
    alphafuncz{loop1}=sort(alpha(temp));
    alphafuncz{loop1}(length(temp)+1)=1;
    nc=[1:length(temp)+1]./(length(temp)+1);
    Ialphacondz(loop1)=trapz(alphafuncz{loop1},nc);
end

temp=polyfit(alphavals(:,2)',Ihypcondalpha-Ihyp,1);

output.Ihyp=Ihyp; 
output.Ialpha=Ialpha;
output.Ihypcondalpha=Ihypcondalpha;
output.SigmaIhypcondalpha=std(Ihypcondalpha);
output.ThetaIhypcondalpha=temp(1);
output.Ialphacondz=Ialphacondz;
output.alphavals=alphavals;
output.Zvals=Zvals;

%Discretized functions as outputs
for loop1=1:numbins
    temp=find(qZ>=quantile_listZ(loop1) & qZ<=quantile_listZ(loop1+1));
    graphIhyp(loop1,1)=median(Z(temp));
    graphA(loop1,1)=length(temp)./(sizer(1)*sizer(2));

    temp=find(qa>=quantile_listalpha(loop1) & qa<=quantile_listalpha(loop1+1));
    graphIalpha(loop1,1)=median(alpha(temp));
end

%set A to the midpoint
graphA(1)=graphA(1)/2;
graphA=cumsum(graphA);
%Take to the end of the axes
graphA(length(graphA)+1)=1;

graphIalpha(length(graphIalpha)+1)=1;
graphIhyp(length(graphIhyp)+1)=1;

output.graphA=graphA;
output.graphIalpha=graphIalpha;
output.graphIhyp=graphIhyp;

if plotflag(1)~=0
for plot_loop=1:length(plotflag)

   if plotflag(plot_loop)==1
   %The first panel
   %Plots conventional hypsometry (black) and 
   %Holder catchment scaling relation grey) in the top panel
   
   %The second panel
   %Plots the conditioned integrals:
   %Ihyp conditioned on alpha in black
   %Ialpha conditioned on z in grey
   
   figure
   subplot(2,1,1)
   plot(graphIhyp,graphA,'-k','linewidth',1)
   hold on
   plot(graphIalpha,graphA,'-','color',[0.55 0.55 0.55],'linewidth',1)
   xlim([0 1]); ylim([0 1.01])
   legend('Hypsometry','Holder scaling')
   xlabel('\alpha*,  z*')
   ylabel('P[A*(z*)],  P[A*(\alpha*)]')
   subplot(2,1,2)
   p1=plot(graphIalpha(1:length(Ihypcondalpha)),Ihypcondalpha,'-xk','linewidth',1)
   hold on
   p2=plot(graphIhyp(1:length(Ialphacondz)),Ialphacondz,'-x','color',[0.55 0.55 0.55],'linewidth',1)
   xlim([0 1]); ylim([0 1])
   plot([0 1],[Ihyp Ihyp],'--k','linewidth',1)
   plot([0 1],[Ialpha Ialpha],'--','color',[0.55 0.55 0.55],'linewidth',1)
   legend([p1 p2],'I_{hyp|\alpha*}','I_{\alpha*|z*}')
   xlabel('\alpha*,  z*')
   ylabel('Conditional measures')
   elseif plotflag(plot_loop)==2
    %The upper panel
    %The hypsometric functions conditioned on alpha 
    %for those closest to 5%, 25% 50%, 75%, 95% values
    
    %The lower panel
    %The Holder exponent catchment area scaling functions conditioned on z
    %for those closest to 5%, 25% 50%, 75%, 95% values
    target=[0.025 0.25 0.5 0.75 0.975];
    tmin=target-1/(2*numbins);
    tmax=target+1/(2*numbins);
    figure
    for loop1=1:length(target)
        temp=find(qa>=tmin(loop1)*sizer(1)*sizer(2) & qa<=tmax(loop1)*sizer(1)*sizer(2));
        graphzfuncalpha=sort(Z(temp));
        nc=[1:length(temp)]./length(temp);
        nc(length(nc)+1)=1;
        graphzfuncalpha(length(graphzfuncalpha)+1)=1;
        output.Ihypcondalpha_graph(loop1)=trapz(graphzfuncalpha,nc);
        output.graphzfuncalpha{loop1}=graphzfuncalpha;
        output.graphzfuncalpha_medianalpha(loop1,1)=median(alpha(temp));
        output.graphzfuncalpha_medianalpha(loop1,2)=median(origa(temp));
        
        subplot(2,1,1)
        hold on
        if loop1==1
            plot(graphzfuncalpha,nc,'-.','color',[0.55 0.55 0.55],'linewidth',1)
        elseif loop1==2
            plot(graphzfuncalpha,nc,'-','color',[0.55 0.55 0.55],'linewidth',1)   
        elseif loop1==3
            plot(graphzfuncalpha,nc,'--k','linewidth',1)   
        elseif loop1==4
            plot(graphzfuncalpha,nc,'-.k','linewidth',1)   
        else
            plot(graphzfuncalpha,nc,'-k','linewidth',1)   
        end
       
        temp=find(qZ>=tmin(loop1)*sizer(1)*sizer(2) & qZ<=tmax(loop1)*sizer(1)*sizer(2));
        graphalphafuncz=sort(alpha(temp));
        nc=[1:length(temp)]./length(temp);
        nc(length(nc)+1)=1;
        graphalphafuncz(length(graphalphafuncz)+1)=1;
        xlim([0 1]); ylim([0 1.01])
        output.Ialphacondz_graph(loop1)=trapz(graphalphafuncz,nc);
        output.graphalphafuncz{loop1}(:,1)=graphalphafuncz;
        output.graphalphafuncz_medianZ(loop1,1)=median(Z(temp));
        output.graphalphafuncz_medianZ(loop1,2)=median(origZ(temp));
        if loop1==5
            legend('\alpha*_{0.025}','\alpha*_{0.25}','\alpha*_{0.5}','\alpha*_{0.75}','\alpha*_{0.975}')
            xlabel('z*')
            ylabel('P[A*(z*|\alpha*)]')
        end
        subplot(2,1,2)
        hold on
        if loop1==1
            plot(graphalphafuncz,nc,'-.','color',[0.55 0.55 0.55],'linewidth',1)
        elseif loop1==2
            plot(graphalphafuncz,nc,'-','color',[0.55 0.55 0.55],'linewidth',1)   
        elseif loop1==3
            plot(graphalphafuncz,nc,'--k','linewidth',1)   
        elseif loop1==4
            plot(graphalphafuncz,nc,'-.k','linewidth',1)   
        else
            plot(graphalphafuncz,nc,'-k','linewidth',1)   
        end
        xlim([0 1]); ylim([0 1.01])
        clear nc graphalphafuncz graphzfuncalpha
        if loop1==5
            legend('z*_{0.025}','z*_{0.25}','z*_{0.5}','z*_{0.75}','z*_{0.975}')
            xlabel('\alpha*')
            ylabel('P[A*(\alpha*|z*)]')
        end
    end
   elseif plotflag(plot_loop)==3
       %2D histogram
       f1=figure;
       H=histogram2(origZ,origa,numbins,'Normalization','pdf');
       valsH=H.Values;
       binsy=H.XBinEdges;
       binsx=H.YBinEdges;
       close(f1)
       filter=ones(5)/25;
       smoothedhist=conv2(valsH,filter,'same');
       smoothedhist=smoothedhist./sum(smoothedhist(:));
       figure
       subplot(1,2,1)
       contourf(smoothedhist,[0.001:0.0005:0.005])
       axis image
       xlabel('\alpha')
       ylabel('z(m)')
       colorbar
       set(gca,'Position',[0.08 0.33 0.34 0.36]);
       xt=xticks;
       yt=yticks;
       xticklabels(round(binsx(xt),2));
       yticklabels(round(binsy(yt),2,'significant'));
       subplot(1,2,2)
       contourf(log10(smoothedhist),[-5:0.5:-2.5])
       axis image
       xlabel('\alpha')
       ylabel('z(m)')
       colorbar
       set(gca,'Position',[0.56 0.33 0.34 0.36]);
       xt=xticks;
       yt=yticks;
       xticklabels(round(binsx(xt),2));
       yticklabels(round(binsy(yt),2,'significant'));
   end
end
end
    