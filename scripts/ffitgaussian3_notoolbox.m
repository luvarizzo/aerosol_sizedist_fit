% Function FitGaussian v3 - no Curve Fitting Toolbox version
% Fit 1, 2 or 3 Gaussian - fittype('a1*exp(-((x-b1)/c1)^2)+a2*exp(-((x-b2)/c2)^2)')
% Update from FitGaussian: uses the last fitting parameters as a first guess
% Written by Luciana Rizzo (lrizzo@usp.br)

% input: ngauss (integer) (number of gaussians you want to fit)
%        manual (logical) (manual choice of fitting parameter starting points)
%        fdiam
%        fsmps
%        extdiam (vector) (extrapolation diameters - use 0 if you don't want extrapolation)
%        prev (previous fitted parameters) (a1 | a2 | a3 | b1 | b2 | b3 | c1 | c2 | c3 )
%
% output: ffits (a1 | a2 | a3 | b1 | b2 | b3 | c1 | c2 | c3 | sse | r^2 | dfe | adjr^2 | rmse)         
%         curva (fitted size dist)
%         curvaext (fitted extrapolated size dist)
%
% fittype('a1*exp(-((x-b1)/c1)^2)+a2*exp(-((x-b2)/c2)^2)+a3*exp(-((x-b3)/c3)^2)')
%
% restrictions: Dpg ultrafine <30nm
%               Dpg Aitken >30nm
%               Dpg Accum >90nm
%               sg between 1.2 and 2.1 for all modes (Hussein et al., 2005)
%               (se o ajuste for manual, não há restrições)  
% changes from v2: 
%               Another option for 3 lognormal fitting, without Nucleation mode (Dpg1>30nm)       
%
%
    function [ffits,curva,curvaext] = ffitgaussian3_notoolbox(ngauss,manual,fdiam,fsmps,extdiam,prev)
    maxDpg1=30; % max Dpg for ultrafine mode, in nm
    minDpgAitken=30; % min Dpg for Aitken mode, in nm
    minDpgAccum=90; % min Dpg for Accum mode, in nm
    maxsg=log10(2.1); % max sg for all modes
    minsg=log10(1.2); % min sg for all modes
%   @Fit 1 Gaussian ----------------------------
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% without curve fitting toolbox
    if (ngauss==1 || ngauss==11)
        y = @(b,fdiam) (b(1).*exp(-((fdiam-b(2))./b(3)).^2));
        OLS = @(b) sum((y(b,fdiam) - fsmps).^2); % Ordinary Least Squares cost function
        opts = optimset('MaxFunEvals',50000, 'MaxIter',10000);
        if(~manual)
            if(prev(6)>1)
                start = [prev(3); prev(6); prev(9);];
            else
                faux=find(fsmps==max(fsmps));
                start(2)=fdiam(faux(1)); %Dpg start point at the sizedist peak
                start(1)=max(fsmps)/4;
                start(3)=log10(2);
                clear faux
            end
            low3=start(2)-0.3;
            up3=start(2)+0.3;
            low = [10; low3; minsg;];
            up = [max(fsmps); up3; maxsg;]; 

        else  % if manual
            disp('Enter StartPoint for Dpg1')
            aux=ginput(1);
            start3=aux(1)
            low3=start3-0.2;   
            up3=start3+0.3;
            start = [max(fsmps)/4; start3; log10(2);];
            low = [10; low3; minsg;];
            up = [max(fsmps); up3; maxsg;];
        end      
        B = fminsearchbnd(OLS, start, low, up, opts);
        ffits=[NaN; NaN; B(1); NaN; NaN; B(2); NaN; NaN; B(3);];
        curva = y(B,fdiam);
        ffits(10) = NaN;  % sum of squares due to error
        [ffits(11),ffits(13)] = rsquared(fsmps,curva,length(B));  % coefficient of determination % degree of freedom adjusted rsquare
        ffits(12) = length(fsmps)-length(B); % degrees of freedom
        ffits(14) = rmse(fsmps,curva);
        if(extdiam(1)>0)
            curvaext = y(B,extdiam);
        else
            curvaext=0;
        end
    end        
 %   
 %  @Fit 2 Gaussian ----------------------------
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% without curve fitting toolbox
    if (ngauss==2 || ngauss==22)
        y = @(b,fdiam) (b(1).*exp(-((fdiam-b(2))./b(3)).^2)) + (b(4).*exp(-((fdiam-b(5))./b(6)).^2));        
        OLS = @(b) sum((y(b,fdiam) - fsmps).^2); % Ordinary Least Squares cost function
        opts = optimset('MaxFunEvals',50000, 'MaxIter',10000);
        if(~manual)
            if(prev(5)>1 && prev(6)>1) % if previous fit parameters are available
                start = [prev(2); prev(5); prev(8); prev(3); prev(6); prev(9);];                    
                if(start(2)<log10(minDpgAitken));start(2)=log10(minDpgAitken);end
                if(start(5)<log10(minDpgAccum));start(5)=log10(minDpgAccum);end
                low2=log10(minDpgAitken); up2=start(2)+0.2;
                low3=log10(minDpgAccum); up3=start(5)+0.2;                
            else % previous fit parameters are not available
                faux=find(fsmps==max(fsmps));
                start = [max(fsmps)/4; fdiam(faux(1)); log10(2); ...   % first shot on Dpg1 and Dpg2
                    max(fsmps)/4; 1.8; log10(2);]; 
                clear faux                
                low2=log10(minDpgAitken);
                low3=log10(minDpgAccum);
                if(start(2)<2) % if Dp(max smps) < 100 nm, conc(Aitken) > conc(Accum)  (also have used Dp<150 = 2.18)
                    up2=start(2)+0.2;
                    start(5)=log10(minDpgAccum)+0.2; up3=max(fdiam); 
                else
                    start(5)=start(2); up3=max(fdiam);  % calibrate better the initial shots
                    start(2)=start(5)-0.2; up2=start(2)+0.2;  
                end               
                if(start(2)<log10(minDpgAitken));start(2)=log10(minDpgAitken);end
                if(start(5)<log10(minDpgAccum));start(5)=log10(minDpgAccum);end
            end
            low = [10; low2; minsg; 10; low3; minsg;];
            up = [max(fsmps); up2; maxsg; max(fsmps); up3; maxsg;];                             
        else  % if manual
            disp('Enter StartPoint for Dpg1')
            aux=ginput(1); 
            start2=aux(1)
            low2=start2-0.2;   
            up2=start2+0.2;
            disp('Enter StartPoint for Dpg2')
            aux=ginput(1);   
            start3=aux(1)
            low3=start3-0.2;   
            up3=start3+0.2;   %up3=max(fdiam);
            start = [max(fsmps)/4; start2; log10(2); ...
                max(fsmps)/4; start3; log10(2);];
            low = [10; low2; minsg; 10; low3; minsg;];
            up = [max(fsmps); up2; maxsg; max(fsmps); up3; maxsg;];
        end  
        B = fminsearchbnd(OLS, start, low, up, opts);
        ffits=[NaN; B(1); B(4); NaN; B(2); B(5); NaN; B(3); B(6);];
        curva = y(B,fdiam);
        ffits(10) = NaN;  % sum of squares due to error
        [ffits(11),ffits(13)] = rsquared(fsmps,curva,length(B));  % coefficient of determination % degree of freedom adjusted rsquare
        ffits(12) = length(fsmps)-length(B); % degrees of freedom
        ffits(14) = rmse(fsmps,curva);
        if(extdiam(1)>0)
            curvaext = y(B,extdiam);
        else
            curvaext=0;
        end
    end    
%   @Fit 3 Gaussian with ultrafine mode (Dpg<25nm) ----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% without curve fitting toolbox
    if (ngauss==3 || ngauss==33)
        % y = @(b,fdiam) b(1).*exp(-((fdiam-b(2))./b(3)).^2);
        % y = @(b,fdiam) (b(1).*exp(-((fdiam-b(2))./b(3)).^2)) + (b(4).*exp(-((fdiam-b(5))./b(6)).^2));
        y = @(b,fdiam) b(1).*exp(-((fdiam-b(2))./b(3)).^2) + b(4).*exp(-((fdiam-b(5))./b(6)).^2) + b(7).*exp(-((fdiam-b(8))./b(9)).^2);             % Objective function
        OLS = @(b) sum((y(b,fdiam) - fsmps).^2); % Ordinary Least Squares cost function
        opts = optimset('MaxFunEvals',50000, 'MaxIter',10000);
        %
        if(~manual)
            if(prev(4)>0 && prev(5)>1 && prev(6)>1) % if previous fit parameters are available
                start = [prev(1); prev(4); prev(7); prev(2); prev(5); prev(8); prev(3); prev(6); prev(9);];                         
            else % previous fit parameters are not available
                start = [max(fsmps)/4; 1.3; log10(2); ...
                    max(fsmps)/4; 1.8; log10(2); ...
                    max(fsmps)/4; 2.3; log10(2);];
            end     
            if(start(2)>log10(maxDpg1));start(2)=log10(maxDpg1);end
            if(start(5)<log10(minDpgAitken));start(5)=log10(minDpgAitken);end
            if(start(8)<log10(minDpgAccum));start(8)=log10(minDpgAccum);end            
            low1=start(2)-0.4; up1=log10(maxDpg1);  % !!! restriction: Dpg ultrafine >8nm & <25nm
            low2=log10(minDpgAitken); up2=start(5)+0.2; 
            low3=log10(minDpgAccum); up3=max(fdiam);     
            low = [10; low1; minsg; 10; low2; minsg; 10; low3; minsg;];
            up = [max(fsmps); up1; maxsg; max(fsmps); up2; maxsg; max(fsmps); up3; maxsg;];                                                
            %
        else  % if manual
            disp('Enter StartPoint for Dpg1')
            aux=ginput(1);
            start1=aux(1)
            low1=start1-0.4;   
            up1=log10(maxDpg1); % restriction for max Dpg ultrafine
            if(start1>log10(maxDpg1))
                start1=log10(maxDpg1);
            end
            disp('Enter StartPoint for Dpg2')
            aux=ginput(1);
            start2=aux(1)
            low2=start2-0.2;   
            up2=start2+0.2;
            disp('Enter StartPoint for Dpg3')
            aux=ginput(1);
            start3=aux(1)
            low3=start3-0.2;   
            up3=max(fdiam);   %up3=max(fdiam); up3=start3+0.2;
            start = [max(fsmps)/4; start1; log10(2); ...
                max(fsmps)/4; start2; log10(2); ...
                max(fsmps)/4; start3; log10(2);];
            low = [10; low1; minsg; 10; low2; minsg; 10; low3; minsg;];
            up = [max(fsmps); up1; maxsg; max(fsmps); up2; maxsg; max(fsmps); up3; maxsg;];        
        end
        B = fminsearchbnd(OLS, start, low, up, opts);
        ffits=[B(1); B(4); B(7); B(2); B(5); B(8); B(3); B(6); B(9)];
        curva = y(B,fdiam);
        ffits(10) = NaN;  % sum of squares due to error
        [ffits(11),ffits(13)] = rsquared(fsmps,curva,length(B));  % coefficient of determination % degree of freedom adjusted rsquare
        ffits(12) = length(fsmps)-length(B); % degrees of freedom
        ffits(14) = rmse(fsmps,curva);
        if(extdiam(1)>0)
            curvaext = y(B,extdiam);
        else
            curvaext=0;
        end
    end   
    %
%   @Fit 3 Gaussian without ultrafine mode (Dpg1>30nm) ----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% without curve fitting toolbox
    if (ngauss==4 || ngauss==44)
        y = @(b,fdiam) b(1).*exp(-((fdiam-b(2))./b(3)).^2) + b(4).*exp(-((fdiam-b(5))./b(6)).^2) + b(7).*exp(-((fdiam-b(8))./b(9)).^2);             % Objective function
        OLS = @(b) sum((y(b,fdiam) - fsmps).^2); % Ordinary Least Squares cost function
        opts = optimset('MaxFunEvals',50000, 'MaxIter',10000);        
        if(~manual)
            if(prev(4)>0 && prev(5)>1 && prev(6)>1) % if previous fit parameters are available
                start = [prev(1); prev(4); prev(7); prev(2); prev(5); prev(8); prev(3); prev(6); prev(9);];                          
            else % previous fit parameters are not available
                start = [max(fsmps)/4; 1.7; log10(2); ... %50nm
                    max(fsmps)/4; 2.0; log10(2); ...      %100nm
                    max(fsmps)/4; 2.3; log10(2);];        %200nm 
            end 
            if(start(2)>log10(minDpgAitken));start(2)=log10(minDpgAitken);end
            if(start(5)<log10(minDpgAccum));start(5)=log10(minDpgAccum);end
            if(start(8)<log10(minDpgAccum));start(8)=log10(minDpgAccum);end            
            low1=log10(minDpgAitken); up1=start(2)+0.2;   % !!! restriction: Dpg Aitken >30nm
            low2=log10(minDpgAccum); up2=start(5)+0.2; 
            low3=log10(minDpgAccum); up3=max(fdiam);     
            low = [10; low1; minsg; 10; low2; minsg; 10; low3; minsg;];
            up = [max(fsmps); up1; maxsg; max(fsmps); up2; maxsg; max(fsmps); up3; maxsg;];         
        else  % if manual
            disp('Enter StartPoint for Dpg1')
            aux=ginput(1);
            start1=aux(1)
            low1=log10(minDpgAitken); % !!! restriction: Dpg Aitken >30nm
            up1=start1+0.2;
            if(start1<log10(minDpgAitken))
                start1=log10(minDpgAitken);
            end
            disp('Enter StartPoint for Dpg2')
            aux=ginput(1);
            start2=aux(1)
            low2=start2-0.2;   
            up2=start2+0.2;
            disp('Enter StartPoint for Dpg3')
            aux=ginput(1);
            start3=aux(1)
            low3=start3-0.2;   
            up3=max(fdiam);   %up3=max(fdiam); up3=startb3+0.2;
            start = [max(fsmps)/4; start1; log10(2); ...
                max(fsmps)/4; start2; log10(2); ...
                max(fsmps)/4; start3; log10(2);];
            low = [10; low1; minsg; 10; low2; minsg; 10; low3; minsg;];
            up = [max(fsmps); up1; maxsg; max(fsmps); up2; maxsg; max(fsmps); up3; maxsg;];                
        end
        B = fminsearchbnd(OLS, start, low, up, opts);
        ffits=[B(1); B(4); B(7); B(2); B(5); B(8); B(3); B(6); B(9)];
        curva = y(B,fdiam);
        ffits(10) = NaN;  % sum of squares due to error
        [ffits(11),ffits(13)] = rsquared(fsmps,curva,length(B));  % coefficient of determination % degree of freedom adjusted rsquare
        ffits(12) = length(fsmps)-length(B); % degrees of freedom
        ffits(14) = rmse(fsmps,curva);
        if(extdiam(1)>0)
            curvaext = y(B,extdiam);
        else
            curvaext=0;
        end
    end
