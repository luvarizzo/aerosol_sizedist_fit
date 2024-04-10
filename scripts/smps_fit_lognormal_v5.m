% smps_fit_lognormal_v8

% Fit one to three lognormal curves to submicrometer aerosol size distribution samples.
% Assumes the following order of modes: nucleation, Aitken, accumulation.
% If only two modes are fitted, they will be considered as Aitken and accumulation.
% Can extrapolate to greater diameters, if needed.
% Can interpolate to a common diameter sequence, if needed.
% It can be run in the fully automatic mode, so that the algorithm decides
% between one, two or three modes.

% Version:
% Identical to v5, but with ffitgaussian3_notoolbox.m, which does not require the matlab Curve fitting toolbox.

% Written by Luciana Rizzo (lrizzo@usp.br)

% Funcionalities: 
%    * option of fully automatic mode, without the need of user intervention
%    * option of manual mode, in which the user marks the start guess of mode diameters in the size distribution plot
%    * option to interpolate (shape-preserving interpolation)
%    * option to extrapolate (based on the lognormal fits)
%    * automatically correct fit order if Dp1>Dp2 or Dp2>Dp3
%    * ffitgaussian3 uses the previous fitting parameters as a first guess
%    * automatically classifies fitted modes into: nucleation(Dpg<30nm); Aitken(30<=Dpg<=90nm); Accum1(90<Dpg<160nm); Accum2(Dpg>=160nm)
%    * saves the sample indexes that did not obtain a satisfactory fit in the automatic mode

% Main changes from v4:
%   * Basic quality standard based on Nfit/Ndata instead of RMSE
%   * Using ffitgaussian3, which restricts Dpg ultrafine<25nm, Dpg Aitken>30nm, Dpg Accum>90nm
%   * Option to fit 3 modes without the nucleation mode, i.e., it is possible to have 1 Aitken (Dpg>30nm) and 2 Accum modes (Dpg>90nm)

% Required functions:
% ffitgaussian3_notoolbox.m
% rsquared.m
% rmse (built-in matlab function since version R2022b) 

% INPUT: data (1st column=index or timeline; 1st row=diam; rows=dN/dlogDp)
%        extradiam (extra diameter array to extrapolate smps data - optional)   
% OUTPUT: extrasmps (extrapolated dN/dlogDp values for the defined extradiam array)
%         smpsnew (if interpolation was made, this output holds the new interpolated dN/dlogDp values)
%         fits col 1-7: matlab time and datevec
%              col 8-16: fit parameters (N1,N2,N3,Dpg1,Dpg2,Dpg3,sg1,sg2,sg3)
%              col 17-21: goodness of fitting (sse,rsquare,dfe,adjrsquare,rmse)           
%              col 22-23: Ndata | Nfit (total concentration comparison)
%         fits_into_5modes col 1-7: matlab time and datevec
%              col 8-22: fit parameters (Nnucl,Naitken1,Naitken2,Naccum1,Naccum2,Dpg1,Dpg2,Dpg3,Dpg4,Dpg5,sg1,sg2,sg3,sg4,sg5)
%              col 23-27: goodness of fitting (sse,rsquare,dfe,adjrsquare,rmse)           
%              col 28-29: Ndata | Nfit (total concentration comparison)
%         rever: sample indexes that did not obtain a satisfactory fit  

% -----------------------------------------------------------------------
% @ User options
% -----------------------------------------------------------------------
onebyone=false; % show plot and prompt for every size distribution fit
showall=true; % show plot for every size distribution fit, without pauses
nonstop=false; % if false, the program prompts for user decision whenever it cannot find a good fit; if true, the program gives up the bad fits and moves to the next size distribution, without prompts

% Do you want to extrapolate the size distribution towards larger diameters?
extra=1; % (yes=1/no=0)
if (extra==1)
    % list the diamters you want to extrapolate
    extradiam=log10([450.0 480.0 500.0 550.0 600.0 ]);
else
    extradiam=0;
end

% Do you want to interpolate the size distribution to a particular diameter sequence?
inter=0; % (yes=1/no=0)
if(inter==1)
    % list the diameter sequence you want to interpolate
    diamnew=[10.7 12.4 14.5 16.9 19.7 23.0 26.9 31.5 36.9 43.3 50.9 60.0 70.8 83.8 99.5 118.8 142.4 171.8 208.7 255.7];
    smpsnew=zeros(size(data,1),length(diamnew));
    smpsnew(1,:)=diamnew;
end

% -----------------------------------------------------------------------
% @ Parameters
%   (may be tuned by the user to improve the algorithm's decisions)
% -----------------------------------------------------------------------
    Nultraf=3; % mininum percentage N10-30/Ntotal to fit three modes
    Naccum=1; % mininum percentage N>100/Ntotal to fit three modes
    rtailmax=15; %10 % maximum percentage dN(last_diam)/dNmax to consider the fit ok
    NfitNdata=0.05; % to consider the fit ok, Nfit/Ndata must be >1-NfitNdata and <1+NfitNdata
    largedisp=2.6;
    superp=0.6; % ratio between DpNucl and DpAitken above which Nucleation and Aitken are considered superposed %0.5
%    maxrmse=400; % maximum rmse to consider the fit ok  % disabled
      
% -----------------------------------------------------------------------
% @ Reading data files 
% -----------------------------------------------------------------------
disp('Reading data files')

% file path of input data
caminho1 = 'C:\...\';
% file path of output data
caminho2 = 'C:\...\';

% file name
nome = 'smps_example.csv'; 
file = readtable([caminho1,nome], "ReadVariableNames",false);
timeline = table2array(file(2:height(file),2)); % the timeline is not used in the calculations, it is just a reference of time
diam = table2array(file(1,3:width(file))); % diameter sequence in nm
sdist = table2array(file(2:height(file),3:width(file))); % normalized size distribution, dN/dlogDp
data=cat(1,log10(diam),sdist);

if(inter==1)
    smpsnew=zeros(size(data,1),length(diamnew));
    smpsnew(1,:)=diamnew;
end

% -----------------------------------------------------------------------
% @ Program starts here
% -----------------------------------------------------------------------           
    fits=zeros(size(data,1),16); 
    extrasmps=zeros(size(data,1),size(extradiam,2));
    extrasmps(1,:)=10.^extradiam;
    rever=0;
    for i=2:size(data,1)
        twomodes=false; force=false; bad=false;
        if (sum(data(i,2:size(data,2)))==0 || any(isnan(data(i,:))))  % desconsider size dist with NaNs
            disp('Sample contains NaNs. Skipped.')
            close; continue
        end   
        disp(i)
        if(showall)
            figure(1)
        end
        plot(data(1,:),data(i,:),'-b.'),title(num2str(i))
        set(gca,'XTick',[1 1.5 1.8 2.0 2.15 2.3 2.5 2.7])
        set(gca,'XTickLabel',floor([10^1 10^1.5 10^1.8 10^2.0 10^2.15 10^2.3 10^2.5 10^2.7]))
        xlabel('Dp (nm)'),ylabel('dN/dlogDp (cm-3)')
        set(1,'WindowStyle','docked')
        hold on        
        dlogdp(2:size(data,2))=diff(data(1,:));
        dlogdp(1)=dlogdp(2);
        aux=find(data(1,:)<log10(30));
        concuf=100*(data(i,aux)*dlogdp(aux)')/(data(i,:)*dlogdp'); % percent concentration of ultrafine (10-40mn)
        aux=find(data(1,:)>log10(100));
        concacc=100*(data(i,aux)*dlogdp(aux)')/(data(i,:)*dlogdp'); % percent concentration of accumulation (>100mn)        
        Ndata=data(i,:)*dlogdp'; % total concentration of measured size distribution
        disp(['%N10-30  ',num2str(concuf),'    %N>100  ',num2str(concacc) ])        
        if(concuf>Nultraf && concacc>Naccum)
            choosefit=3;
        else
            choosefit=2;
            %choosefit=input('Fit one gaussian (1) or fit two gaussian (2) or fit three gaussian (3) or skip fitting (0)?   ');
        end
        clear aux dlogdp
        if(choosefit==0) % Disabled. 
            close; continue
        else
            manual=false;
            diam=data(1,1:size(data,2));
            smps=data(i,1:size(data,2));               
            [aux1,aux2,aux3]=ffitgaussian3_notoolbox(choosefit,manual,diam,smps,extradiam,fits(i-1,1:9)); % Call fit gaussian function
            fits(i,1:14)=aux1;  
            dlogdpfit(2:size(diam,2))=diff(diam(1,:));
            dlogdpfit(1)=dlogdpfit(2);
            Nfit=aux2*dlogdpfit'; % total concentration of fitted size distribution
            fits(i,15)=Ndata; fits(i,16)=Nfit;
            plot(diam,aux2,'-r')
            if(extra==1)
                extrasmps(i,:)=aux3;
                plot(extradiam,extrasmps(i,:),'ro')
            end
            legend('data',num2str(choosefit))
            hold off 
            disp([num2str(choosefit),' Gaussian fit'])
            disp(floor([fits(i,1),fits(i,2),fits(i,3)]))
            disp(floor([10^fits(i,4),10^fits(i,5),10^fits(i,6)]))
            disp(floor([10^fits(i,7)*100,10^fits(i,8)*100,10^fits(i,9)*100]))
            disp(['AdjR2 = ',num2str(fits(i,13)),' RMSE = ',num2str(fits(i,14))])
            disp(['Nfit/Ndata = ',num2str(fits(i,16)/fits(i,15))])
            if(extra==1)
                rtail=100*extrasmps(i,size(extrasmps,2))/max(smps);
            else
                rtail=100*smps(1,size(smps,2))/max(smps);
            end
            disp(['Right tail = ',num2str(rtail),' %'])
            clear aux1 aux2 aux3
            % Conditions to decide whether the fitting is ok
            ok=1;
            direto=false;        
            % Basic quality standard: rtail e Nfit/Ndata 
            if(rtail<rtailmax && Nfit/Ndata>1-NfitNdata && Nfit/Ndata<1+NfitNdata) 
                if(choosefit==3 && 10^(fits(i,4)-fits(i,5))>superp) % 1st criteria: Nucleation and Aitken are superposed, so there is not an nucleation mode  
                        disp('*** Warning: Nucleation and Aitken superposed --> trying 2 modes fitting')
                        hold on
                        manual=false;
                        Nfitantes=Nfit;
                        [aux1,aux2,aux3]=ffitgaussian3_notoolbox(2,manual,diam,smps,extradiam,fits(i-1,1:9)); % Call fit gaussian function
                        fits(i,1:14)=aux1; 
                        dlogdpfit(2:size(diam,2))=diff(diam(1,:));
                        dlogdpfit(1)=dlogdpfit(2);
                        Nfit=aux2*dlogdpfit'; % total concentration of fitted size distribution
                        fits(i,15)=Ndata; fits(i,16)=Nfit;
                        plot(diam,aux2,'-g')
                        if(extra==1)
                            extrasmps(i,:)=aux3;
                            plot(extradiam,extrasmps(i,:),'ro')
                        end
                        legend('data','3','2')
                        hold off
                        if(extra==1)
                            rtail=100*extrasmps(i,size(extrasmps,2))/max(smps);
                        else
                            rtail=100*smps(1,size(smps,2))/max(smps);
                        end
                        disp('2 Gaussian fit')
                        disp(floor([fits(i,1),fits(i,2),fits(i,3)]))
                        disp(floor([10^fits(i,4),10^fits(i,5),10^fits(i,6)]))
                        disp(floor([10^fits(i,7)*100,10^fits(i,8)*100,10^fits(i,9)*100]))
                        disp(['AdjR2 = ',num2str(fits(i,13)),' RMSE = ',num2str(fits(i,14))])
                        disp(['Nfit/Ndata = ',num2str(fits(i,16)/fits(i,15))])
                        disp(['% right tail =',num2str(rtail)])
                        if(rtail<rtailmax && Nfit/Ndata>1-NfitNdata && Nfit/Ndata<1+NfitNdata) %fits(i,14)<maxrmse)
                            if(abs(1-Nfit/Ndata)<abs(1-Nfitantes/Ndata)+0.02)
                                ok=1; % Test if the 2 mode fitting is better compared to the 3 mode fitting
                            else % If not, returns to the 3 mode fitting and may prompt for a user decision
                                disp('*** Tryed 2 modes fitting (Nucleation and Aitken superposed), but 3 modes looks better')
                                disp(['Nfit/Ndata(3modes)=',num2str(Nfitantes/Ndata),'; Nfit/Ndata(2modes)=',num2str(Nfit/Ndata)])
                                manual=false;
                                Nfitantes=Nfit;
                                [aux1,aux2,aux3]=ffitgaussian3_notoolbox(3,manual,diam,smps,extradiam,fits(i-1,1:9)); % Call fit gaussian function
                                fits(i,1:14)=aux1;
                                dlogdpfit(2:size(diam,2))=diff(diam(1,:));
                                dlogdpfit(1)=dlogdpfit(2);
                                Nfit=aux2*dlogdpfit'; % total concentration of fitted size distribution
                                fits(i,15)=Ndata; fits(i,16)=Nfit;
                                if(extra==1)
                                    rtail=100*extrasmps(i,size(extrasmps,2))/max(smps);
                                else
                                    rtail=100*smps(1,size(smps,2))/max(smps);
                                end
                                if(nonstop)
                                    ok=0;
                                else
                                    ok=2;
                                    choosefit=3;
                                end
                            end
                        else
                            disp('*** 2 modes did not work, return to 3 modes and move on')
                            disp(['Nfit/Ndata(3modes)=',num2str(Nfitantes/Ndata),'; Nfit/Ndata(2mode)=',num2str(Nfit/Ndata)])
                            choosefit=3;
                            manual=false;
                            Nfitantes=Nfit;
                            [aux1,aux2,aux3]=ffitgaussian3_notoolbox(2,manual,diam,smps,extradiam,fits(i-1,1:9)); % Call fit gaussian function
                            fits(i,1:14)=aux1;
                            dlogdpfit(2:size(diam,2))=diff(diam(1,:));
                            dlogdpfit(1)=dlogdpfit(2);
                            Nfit=aux2*dlogdpfit'; % total concentration of fitted size distribution
                            fits(i,15)=Ndata; fits(i,16)=Nfit;
                            if(extra==1)
                                extrasmps(i,:)=aux3;
                            end
                            if(extra==1)
                                rtail=100*extrasmps(i,size(extrasmps,2))/max(smps);
                            else
                                rtail=100*smps(1,size(smps,2))/max(smps);
                            end
                            ok=1;                            
                        end                        
                end
                if(choosefit==2) 
                    if (fits(i,3)/fits(i,2)>10 || fits(i,3)/fits(i,2)<0.1)   % 2nd criteria: fit only one mode?    
                        Nfitantes=Nfit;
                        disp('*** Warning: One mode might be more appropriate --> trying 1 mode fitting')
                        hold on
                        manual=false;
                        [aux1,aux2,aux3]=ffitgaussian3_notoolbox(1,manual,diam,smps,extradiam,fits(i-1,1:9)); % Call fit gaussian function
                        fits(i,1:14)=aux1; 
                        dlogdpfit(2:size(diam,2))=diff(diam(1,:));
                        dlogdpfit(1)=dlogdpfit(2);
                        Nfit=aux2*dlogdpfit'; % total concentration of fitted size distribution
                        fits(i,15)=Ndata; fits(i,16)=Nfit;
                        plot(diam,aux2,'-g')
                        if(extra==1)
                            extrasmps(i,:)=aux3;
                            plot(extradiam,extrasmps(i,:),'go')
                        end
                        legend('data','2','1')
                        hold off
                        if(extra==1)
                            rtail=100*extrasmps(i,size(extrasmps,2))/max(smps);
                        else
                            rtail=100*smps(1,size(smps,2))/max(smps);
                        end
                        disp('1 Gaussian fit')
                        disp(floor([fits(i,1),fits(i,2),fits(i,3)]))
                        disp(floor([10^fits(i,4),10^fits(i,5),10^fits(i,6)]))
                        disp(floor([10^fits(i,7)*100,10^fits(i,8)*100,10^fits(i,9)*100]))
                        disp(['AdjR2 = ',num2str(fits(i,13)),' RMSE = ',num2str(fits(i,14))])
                        disp(['Nfit/Ndata = ',num2str(fits(i,16)/fits(i,15))])
                        disp(['% right tail =',num2str(rtail)])
                        choosefit=1;
                        if(rtail<rtailmax && Nfit/Ndata>1-NfitNdata && Nfit/Ndata<1+NfitNdata) %fits(i,14)<maxrmse)
                            if(abs(1-Nfit/Ndata)<abs(1-Nfitantes/Ndata)+0.01)
                                ok=1; % Test if the 1 mode fitting is better compared to the 2 mode fitting
                            else % If not, returns to the 2 mode fitting
                                disp('*** Tryed 1 mode fitting, but 2 modes is better')
                                disp(['Nfit/Ndata(2modes)=',num2str(Nfitantes/Ndata),'; Nfit/Ndata(1mode)=',num2str(Nfit/Ndata)])
                                choosefit=2;
                                manual=false;
                                Nfitantes=Nfit;
                                [aux1,aux2,aux3]=ffitgaussian3_notoolbox(2,manual,diam,smps,extradiam,fits(i-1,1:9)); % Call fit gaussian function
                                fits(i,1:14)=aux1;
                                dlogdpfit(2:size(diam,2))=diff(diam(1,:));
                                dlogdpfit(1)=dlogdpfit(2);
                                Nfit=aux2*dlogdpfit'; % total concentration of fitted size distribution
                                fits(i,15)=Ndata; fits(i,16)=Nfit;

                                if(extra==1)
                                    extrasmps(i,:)=aux3;
                                     plot(extradiam,extrasmps(i,:),'ro')
                                end
                                 legend('data',num2str(2))
                                 hold off
                                if(extra==1)
                                    rtail=100*extrasmps(i,size(extrasmps,2))/max(smps);
                                else
                                    rtail=100*smps(1,size(smps,2))/max(smps);
                                end                              
                                ok=1;
                            end
                        else
                            disp('*** 1 mode does not work, return to 2 modes and move on')
                            disp(['Nfit/Ndata(2modes)=',num2str(Nfitantes/Ndata),'; Nfit/Ndata(1mode)=',num2str(Nfit/Ndata)])
                            choosefit=2;
                            manual=false;
                            Nfitantes=Nfit;
                            [aux1,aux2,aux3]=ffitgaussian3_notoolbox(2,manual,diam,smps,extradiam,fits(i-1,1:9)); % Call fit gaussian function
                            fits(i,1:14)=aux1;
                            dlogdpfit(2:size(diam,2))=diff(diam(1,:));
                            dlogdpfit(1)=dlogdpfit(2);
                            Nfit=aux2*dlogdpfit'; % total concentration of fitted size distribution
                            fits(i,15)=Ndata; fits(i,16)=Nfit;
                            if(extra==1)
                                extrasmps(i,:)=aux3;
                            end
                            if(extra==1)
                                rtail=100*extrasmps(i,size(extrasmps,2))/max(smps);
                            else
                                rtail=100*smps(1,size(smps,2))/max(smps);
                            end
                            ok=1;
                        end                        
                    end
                end                              
                if(ok==1)
                    disp('ok!')
                    if(fits(i,4)>fits(i,5) || fits(i,5)>fits(i,6))                        
                        %disp('*** Warning: Dp1>Dp2 or Dp2>Dp3')
                        aux=fits(i,:);
                        if(fits(i,4)>fits(i,5))
                            fits(i,1)=aux(1,2);fits(i,2)=aux(1,1);
                            fits(i,4)=aux(1,5);fits(i,5)=aux(1,4);
                            fits(i,7)=aux(1,8);fits(i,8)=aux(1,7);
                        end
                        if(fits(i,5)>fits(i,6))
                            fits(i,2)=aux(1,3);fits(i,3)=aux(1,2);
                            fits(i,5)=aux(1,6);fits(i,6)=aux(1,5);
                            fits(i,8)=aux(1,9);fits(i,9)=aux(1,8);
                        end
                        clear aux
                    end
                    if(~direto && onebyone)
                        disp('Press any key to move to the next size distribution fit')
                        pause
                    end          
                end
            else
                disp('*** Fit do not reach quality standards')
                if(rtail>rtailmax)
                    disp(['% right tail > ',num2str(rtailmax)])
                    disp('Forcing right tail to zero')
                    force=true;
                    plot(data(1,:),data(i,:),'-b.'),title(num2str(i))
                    set(gca,'XTick',[1 1.5 1.8 2.0 2.15 2.3 2.5 2.7])
                    set(gca,'XTickLabel',floor([10^1 10^1.5 10^1.8 10^2.0 10^2.15 10^2.3 10^2.5 10^2.7]))
                    xlabel('Dp (nm)'),ylabel('dN/dlogDp (cm-3)')
                    set(1,'WindowStyle','docked')
                    hold on
                    manual=false;
                    diam_=diam; diam_(length(diam)+1)=extradiam(length(extradiam));
                    smps_=smps; smps_(length(smps)+1)=0;
                    [aux1,aux2,aux3]=ffitgaussian3_notoolbox(choosefit,manual,diam_,smps_,extradiam,fits(i-1,1:9)); % Call fit gaussian function
                    clear diam_ smps_
                    fits(i,1:14)=aux1; 
                    dlogdpfit(2:size(diam,2))=diff(diam(1,:));
                    dlogdpfit(1)=dlogdpfit(2);
                    Nfit=aux2(1:length(aux2)-1)*dlogdpfit'; % total concentration of fitted size distribution
                    fits(i,15)=Ndata; fits(i,16)=Nfit;
                    plot(diam,aux2(1:length(aux2)-1),'-r')
                    if(extra==1)
                        extrasmps(i,:)=aux3;
                        plot(extradiam,extrasmps(i,:),'ro')
                    end
                    legend('data',num2str(choosefit))
                    hold off
                    disp([num2str(choosefit),' Gaussian fit'])
                    if(extra==1)
                        rtail=100*extrasmps(i,size(extrasmps,2))/max(smps);
                    else
                        rtail=100*smps(1,size(smps,2))/max(smps);
                    end
                    disp(floor([fits(i,1),fits(i,2),fits(i,3)]))
                    disp(floor([10^fits(i,4),10^fits(i,5),10^fits(i,6)]))
                    disp(floor([10^fits(i,7)*100,10^fits(i,8)*100,10^fits(i,9)*100]))
                    disp(['AdjR2 = ',num2str(fits(i,13)),' RMSE = ',num2str(fits(i,14))])
                    disp(['Nfit/Ndata = ',num2str(fits(i,16)/fits(i,15))])
                    disp(['% right tail =',num2str(rtail)])                    
                    if(fits(i,4)>fits(i,5) || fits(i,5)>fits(i,6))
                        aux=fits(i,:);
                        if(fits(i,4)>fits(i,5))
                            fits(i,1)=aux(1,2);fits(i,2)=aux(1,1);
                            fits(i,4)=aux(1,5);fits(i,5)=aux(1,4);
                            fits(i,7)=aux(1,8);fits(i,8)=aux(1,7);
                        end
                        if(fits(i,5)>fits(i,6))
                            fits(i,2)=aux(1,3);fits(i,3)=aux(1,2);
                            fits(i,5)=aux(1,6);fits(i,6)=aux(1,5);
                            fits(i,8)=aux(1,9);fits(i,9)=aux(1,8);
                        end
                        clear aux
                    end                           
                end
                if(rtail>rtailmax || Nfit/Ndata<1-NfitNdata || Nfit/Ndata>1+NfitNdata) %fits(i,14)<maxrmse)
                    bad=true;
                    disp(['Nfit/Ndata = ',num2str(Nfit/Ndata)])
                end
                if(choosefit==2 && bad) % If the fit is bad with 2 modes, try 3 modes and ask user if the result is fine
                    disp('*** Fitting is bad with 2 modes, trying 3 modes instead')
                    if(concuf<2)
                        choosefit=4; % 3 modes without nucleation mode
                    else
                        choosefit=3; % 3 modes with nucleation mode
                    end
                    plot(data(1,:),data(i,:),'-b.'),title(num2str(i))
                    set(gca,'XTick',[1 1.5 1.8 2.0 2.15 2.3 2.5 2.7])
                    set(gca,'XTickLabel',floor([10^1 10^1.5 10^1.8 10^2.0 10^2.15 10^2.3 10^2.5 10^2.7]))
                    xlabel('Dp (nm)'),ylabel('dN/dlogDp (cm-3)')
                    set(1,'WindowStyle','docked')
                    hold on
                    manual=false;
                    if(force)
                        disp('Forcing right tail to zero')
                        diam_=diam; diam_(length(diam)+1)=extradiam(length(extradiam));
                        smps_=smps; smps_(length(smps)+1)=0;
                        [aux1,aux2,aux3]=ffitgaussian3_notoolbox(choosefit,manual,diam_,smps_,extradiam,fits(i-1,1:9)); % Call fit gaussian function
                        clear diam_ smps_
                        aux2(length(aux2))=[];
                    else    
                        [aux1,aux2,aux3]=ffitgaussian3_notoolbox(choosefit,manual,diam,smps,extradiam,fits(i-1,1:9)); % Call fit gaussian function
                    end    
                    fits(i,1:14)=aux1;
                    dlogdpfit(2:size(diam,2))=diff(diam(1,:));
                    dlogdpfit(1)=dlogdpfit(2);
                    Nfit=aux2*dlogdpfit'; % total concentration of fitted size distribution
                    fits(i,15)=Ndata; fits(i,16)=Nfit;
                    plot(diam,aux2,'-r')
                    if(extra==1)
                        extrasmps(i,:)=aux3;
                        plot(extradiam,extrasmps(i,:),'ro')
                    end
                    legend('data',num2str(3))
                    hold off
                    disp('3 Gaussian fit without nucleation mode')
                    if(extra==1)
                        rtail=100*extrasmps(i,size(extrasmps,2))/max(smps);
                    else
                        rtail=100*smps(1,size(smps,2))/max(smps);
                    end
                    disp(floor([fits(i,1),fits(i,2),fits(i,3)]))
                    disp(floor([10^fits(i,4),10^fits(i,5),10^fits(i,6)]))
                    disp(floor([10^fits(i,7)*100,10^fits(i,8)*100,10^fits(i,9)*100]))
                    disp(['AdjR2 = ',num2str(fits(i,13)),' RMSE = ',num2str(fits(i,14))])
                    disp(['Nfit/Ndata = ',num2str(fits(i,16)/fits(i,15))])
                    disp(['% right tail =',num2str(rtail)])                    
                    if(fits(i,4)>fits(i,5) || fits(i,5)>fits(i,6))
                        aux=fits(i,:);
                        if(fits(i,4)>fits(i,5))
                            fits(i,1)=aux(1,2);fits(i,2)=aux(1,1);
                            fits(i,4)=aux(1,5);fits(i,5)=aux(1,4);
                            fits(i,7)=aux(1,8);fits(i,8)=aux(1,7);
                        end
                        if(fits(i,5)>fits(i,6))
                            fits(i,2)=aux(1,3);fits(i,3)=aux(1,2);
                            fits(i,5)=aux(1,6);fits(i,6)=aux(1,5);
                            fits(i,8)=aux(1,9);fits(i,9)=aux(1,8);
                        end
                        clear aux
                    end
                end
                disp('Increased number of modes from 2 to 3')
                if(rtail<rtailmax && Nfit/Ndata>1-NfitNdata && Nfit/Ndata<1+NfitNdata) %fits(i,14)<maxrmse) 
                    ok=1;
                else
                    ok=2;
                    if(nonstop)
                        ok=0;
                    else
                        ok=input('Is the fitting ok?(1) Try again?(2) Give up this fitting?(0)');
                    end
                end                
            end
            %
            if(ok==0)
                if(extra==1);extrasmps(i,:)=0; end
                trying=false;
                fits(i,:)=0;
                if(nonstop)
                    rever=cat(1,rever,i); % keep track of the samples that need review
                end
                continue;
            end
            if (ok==1 || ok==11) 
                trying=false; 
%                  disp('Press any key to continue')
%                  pause
            end 
            if (ok==2 || ok==22)
                trying=true; manual=true;                
                while(trying)        
%                    plot(data(1,:),data(i,:),'-b.')
%                    title(num2str(i))
                    hold on
 %                   set(gca,'XTick',[1 1.5 1.8 2.0 2.15 2.3 2.5 2.7])
 %                   set(gca,'XTickLabel',floor([10^1 10^1.5 10^1.8 10^2.0 10^2.15 10^2.3 10^2.5 10^2.7]))
 %                   xlabel('Dp (nm)'),ylabel('dN/dlogDp (cm-3)')
                    if(twomodes)
                        nmodes=2;
                    else
                        disp('Fit one gaussian (1) or fit two gaussian Dpg2>30nm (2)')
                        disp(' or fit three gaussian with nucleation mode Dpg1<25nm(3)')
                        nmodes=input('or fit three gaussian without nucleation mode Dpg1>30nm (4)?');
                    end
                    twomodes=false;
                    [aux1,aux2,aux3]=ffitgaussian3_notoolbox(nmodes,manual,diam,smps,extradiam,zeros(1,9)); % Call fit gaussian function
                    fits(i,1:14)=aux1;  
                    dlogdpfit(2:size(diam,2))=diff(diam(1,:));
                    dlogdpfit(1)=dlogdpfit(2);
                    Nfit=aux2*dlogdpfit'; % total concentration of fitted size distribution
                    fits(i,15)=Ndata; fits(i,16)=Nfit;
                    plot(diam,aux2,'-c')
                    if(extra==1)
                        extrasmps(i,:)=aux3;
                        plot(extradiam,extrasmps(i,:),'ro')
                    end
                    legend('data',num2str(nmodes))
                    hold off
                    if(extra==1)
                        rtail=100*extrasmps(i,size(extrasmps,2))/max(smps);
                    else
                        rtail=100*smps(1,size(smps,2))/max(smps);
                    end                    
                    disp([num2str(nmodes),' Gaussian fit'])
                    disp(floor([fits(i,1),fits(i,2),fits(i,3)]))
                    disp(floor([10^fits(i,4),10^fits(i,5),10^fits(i,6)]))
                    disp(floor([10^fits(i,7)*100,10^fits(i,8)*100,10^fits(i,9)*100]))
                    disp(['AdjR2 = ',num2str(fits(i,13)),' RMSE = ',num2str(fits(i,14))])
                    disp(['Nfit/Ndata = ',num2str(fits(i,16)/fits(i,15))])
                    disp(['% right tail =',num2str(rtail)])
                    clear aux1 aux2 aux3
                    % Conditions to decide wether the fitting is ok
                    if(rtail<rtailmax && Nfit/Ndata>1-NfitNdata && Nfit/Ndata<1+NfitNdata) %fits(i,14)<maxrmse) 
                        if(fits(i,6)<log10(90) && nmodes==3)
                            ok=2;
                            disp('*** Accum Diam < 90 nm --> choose 2 modes fitting')
                        elseif(fits(i,4)>fits(i,5) || fits(i,5)>fits(i,6))
                            %ok=2;
                            %disp('*** Dp1>Dp2 or Dp2>Dp3')
                            aux=fits(i,:);
                            if(fits(i,4)>fits(i,5))
                                fits(i,1)=aux(1,2);fits(i,2)=aux(1,1);
                                fits(i,4)=aux(1,5);fits(i,5)=aux(1,4);
                                fits(i,7)=aux(1,8);fits(i,8)=aux(1,7);
                            end
                            if(fits(i,5)>fits(i,6))
                                fits(i,2)=aux(1,3);fits(i,3)=aux(1,2);
                                fits(i,5)=aux(1,6);fits(i,6)=aux(1,5);
                                fits(i,8)=aux(1,9);fits(i,9)=aux(1,8);
                            end
                            clear aux
                        elseif(fits(i,1)/fits(i,2)>0.75 &&  10^(fits(i,4)-fits(i,5))>superp)
                            ok=2;
                            disp('*** Warning: Nucleation and Aitken superposed')
                        elseif(fits(i,7)>log10(largedisp) || fits(i,8)>log10(largedisp) || fits(i,9)>log10(largedisp))
                            ok=2;
                            disp(['*** Warning: Dispersion of one of the modes is > ',num2str(largedisp)])
                        else
                            disp('ok!')
                            disp(['Right tail = ',num2str(rtail),'% '])
                        end
                    else
                        disp('*** Fit do not reach quality standards')
                        if(rtail>rtailmax)
                            disp(['% right tail > ',num2str(rtailmax)])
                        end
                        if(Nfit/Ndata>1-NfitNdata && Nfit/Ndata<1+NfitNdata) %fits(i,14)<maxrmse) 
                            disp(['Nfit/Ndata out of range',num2str(Nfit/Ndata)])
                        end
                    end
                    ok=input('Is the fitting ok?(1) Try again?(2) Give up this fitting?(0)');
                    if(ok==0)
                        if(extra==1);extrasmps(i,:)=0; end                        
                        trying=false;
                        fits(i,:)=0;
                    end
                    if (ok==1 || ok==11); trying=false; end                    
                end                
            end
        end    
        if(inter==1 && ok>0)  % interpolate            
            smpsnew(i,1:length(diamnew))=interp1(data(1,:),data(i,:),log10(diamnew),'pchip');
            a=find(smpsnew(i,:)<0 | isnan(smpsnew(i,:))); % Replace eventual negative interpolation values per 0
            smpsnew(i,a)=0;                          
        end  
% Saves progress every 3000 samples
%         if(mod(i,3000)==0)
%             save(['fits',num2str(i),'.txt'],'fits','-ASCII')
%             if(inter==1);save(['smpsnew',num2str(i),'.txt'],'smpsnew','-ASCII');end
%             save(['extrasmps',num2str(i),'.txt'],'extrasmps','-ASCII')
%         end
    end  
    fits(1,:)=[];%
    fits=cat(2,timeline,datevec(timeline),fits);    
%
% convert a1,a2,... fitting parameters to N1,N2,N3,Dpg1,Dpg2,Dpg3,s1,s2,s3
    fits(:,8)=fits(:,8).*fits(:,14)*sqrt(pi); %N1
    fits(:,9)=fits(:,9).*fits(:,15)*sqrt(pi); %N2
    fits(:,10)=fits(:,10).*fits(:,16)*sqrt(pi); %N3
    fits(:,11)=10.^fits(:,11); %Dpg1
    fits(:,12)=10.^fits(:,12); %Dpg2
	fits(:,13)=10.^fits(:,13); %Dpg3
    fits(:,14)=10.^(fits(:,14)/sqrt(2)); %sg1
	fits(:,15)=10.^(fits(:,15)/sqrt(2)); %sg2
	fits(:,16)=10.^(fits(:,16)/sqrt(2)); %sg3
%	
%   Classifies fitted modes into: nucleation(Dpg<30nm);
%   Aitken1(30<=Dpg<60nm); Aitken2(60<=Dpg<90nm)
%   Accum1(90<=Dpg<160nm); Accum2(Dpg>=160nm) 
    fitsnew=cat(2,fits(:,1:7),nan(size(fits,1),15),fits(:,17:23));
    for i=1:size(fits,1)
        disp(i)
        if(fits(i,11)<30 && fits(i,8)>0) % Nucleation
            fitsnew(i,8)=fits(i,8);
            fitsnew(i,13)=fits(i,11);
            fitsnew(i,18)=fits(i,15);
        elseif(fits(i,11)>=30&& fits(i,11)<60) % Aitken1            
            fitsnew(i,9)=fits(i,8);
            fitsnew(i,14)=fits(i,11);
            fitsnew(i,19)=fits(i,15);
        elseif(fits(i,11)>=60&& fits(i,11)<90) % Aitken2
            fitsnew(i,10)=fits(i,8);
            fitsnew(i,15)=fits(i,11);
            fitsnew(i,20)=fits(i,15);
        elseif(fits(i,11)>=90 && fits(i,11)<160) % Accum1
            fitsnew(i,11)=fits(i,8);
            fitsnew(i,16)=fits(i,11);
            fitsnew(i,21)=fits(i,15);
        elseif(fits(i,11)>=160) % Accum2
            fitsnew(i,12)=fits(i,8);
            fitsnew(i,17)=fits(i,11);
            fitsnew(i,22)=fits(i,15);
        end
        if(fits(i,12)<30 && fits(i,9)>0) % Nucleation
            fitsnew(i,8)=fits(i,9);
            fitsnew(i,13)=fits(i,12);
            fitsnew(i,18)=fits(i,15);
        elseif(fits(i,12)>=30 && fits(i,12)<60) % Aitken1
            fitsnew(i,9)=fits(i,9);
            fitsnew(i,14)=fits(i,12);
            fitsnew(i,19)=fits(i,15);
        elseif(fits(i,12)>=60 && fits(i,12)<90) % Aitken2
            fitsnew(i,10)=fits(i,9);
            fitsnew(i,15)=fits(i,12);
            fitsnew(i,20)=fits(i,15);
        elseif(fits(i,12)>=90 && fits(i,12)<160) % Accum1
            fitsnew(i,11)=fits(i,9);
            fitsnew(i,16)=fits(i,12);
            fitsnew(i,21)=fits(i,15);
        elseif(fits(i,12)>=160) % Accum2
            fitsnew(i,12)=fits(i,9);
            fitsnew(i,17)=fits(i,12);
            fitsnew(i,22)=fits(i,15);
        end
        if(fits(i,13)<30 && fits(i,10)>0) % Nucleation
            fitsnew(i,8)=fits(i,10);
            fitsnew(i,13)=fits(i,13);
            fitsnew(i,18)=fits(i,16);
        elseif(fits(i,13)>=30 && fits(i,13)<60) % Aitken1
            fitsnew(i,9)=fits(i,10);
            fitsnew(i,14)=fits(i,13);
            fitsnew(i,19)=fits(i,16);
        elseif(fits(i,13)>=60 && fits(i,13)<90) % Aitken2
            fitsnew(i,10)=fits(i,10);
            fitsnew(i,15)=fits(i,13);
            fitsnew(i,20)=fits(i,16);        
        elseif(fits(i,13)>=90 && fits(i,13)<160) % Accum1
            fitsnew(i,11)=fits(i,10);
            fitsnew(i,16)=fits(i,13);
            fitsnew(i,21)=fits(i,16);
        elseif(fits(i,13)>=160) % Accum2
            fitsnew(i,12)=fits(i,10);
            fitsnew(i,17)=fits(i,13);
            fitsnew(i,22)=fits(i,16);
        end
    end  

% -----------------------------------------------------------------------
% @ Save output files
% -----------------------------------------------------------------------
    
save([caminho2,'fits_',nome(1:length(nome)-4),'.txt'],'fits','-ASCII')
save([caminho2,'fits_',nome(1:length(nome)-4),'.mat'],'fits','-mat')
save([caminho2,'fits_into_5modes_',nome(1:length(nome)-4),'.txt'],'fitsnew','-ASCII')
save([caminho2,'fits_into_5modes_',nome(1:length(nome)-4),'.mat'],'fitsnew','-mat')
if (inter==1 && extra==1)
    sdistnew=cat(2,cat(1,0,timeline),smpsnew,extrasmps);
    save([caminho2,'newsizedist.txt'],'sdistnew','-ASCII')
elseif (inter==1 && extra==0)
    save([caminho2,'interpsmps.txt'],'smpsnew','-ASCII')
elseif (inter==0 && extra==1)
    sdistnew=cat(2,cat(1,0,timeline),data,extrasmps);
    sdistnew(1,2:length(diam)) = 10.^sdistnew(1,2:length(diam));
    save([caminho2,'extrapsizedist.txt'],'sdistnew','-ASCII')
    save([caminho2,'extrapsizedist.mat'],'sdistnew','-mat')
end
save([caminho2,'sample_indexes_bad_fit.txt'],'rever','-ASCII')

% Loop to verify the interpolation and extrapolation results
%     figure
%     for i=2:size(data,1)
%         semilogx(10.^diam,data(i,:),'b*'),grid on, hold on
%         semilogx(sdistnew(1,2:size(sdistnew,2)),sdistnew(i,2:size(sdistnew,2)),'-ro'),hold off
%         title(num2str(i))
%         pause
%     end
