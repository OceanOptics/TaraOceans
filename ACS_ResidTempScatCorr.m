
function [ap_TSalScatCorr] = ACS_ResidTempScatCorr(ap_uncorr, cp, wavel)
    
    wavel = wavel(:)';
    
    [nspect, nwavel] = size(ap_uncorr);

    if size(ap_uncorr)~=size(cp)
        error('ap,cp size mismatch')
    end
    if length(wavel)~=nwavel
        error('Invalid lambda')
    end

    tmp = xlsread('Sullivan_etal_2006_instrumentspecific.xls');
    psiT = interp1(tmp(:,1),tmp(:,2), wavel);
    
    opts = optimset('fminsearch');      
    %opts = optimset(opts,'NonlEqnAlgorithm', 'gn'); 
    opts = optimset(opts,'MaxIter',20000000); 
    opts = optimset(opts,'MaxFunEvals',20000); 
    opts = optimset(opts,'TolX',1e-8);
    opts = optimset(opts,'TolFun',1e-8);
    
    NIR = find(wavel>=710 & wavel<=750);  % spectral srange for optimization (710 to 750nm)
    ref = find(wavel>=730, 1,'first');

    ap_TSalScatCorr = nan(size(ap_uncorr));
    deltaT = nan(size(ap_uncorr,1),1);
    fiterr = nan(size(ap_uncorr,1),1);
   
    for k = 1:nspect

        if all(isfinite(ap_uncorr(k,:)))
            % guess for paramters (beamc at lambda0, beamc slope)
            delT = 0;

            % minimization routine
            [deltaT(k), fiterr(k)] = fminsearch(@f_TS, 0, opts, ap_uncorr(k,:), cp(k,:), psiT, NIR, ref);

            bp = cp(k,:) - ap_uncorr(k,:);
            ap_TSalScatCorr(k,:) = ap_uncorr(k,:) - psiT.*deltaT(k) - ...
                    ((ap_uncorr(k,ref)-psiT(ref).*deltaT(k))./bp(ref)).*bp;
            
            disp(['  Residual T,S,Scat correction: k=' ...
                num2str(k) ' deltaT=' num2str(deltaT(k)) ' err=' num2str(fiterr(k))]);
			
        end

    end




return






function costf = f_TS(delT, ap, cp, psiT, NIR, ref);

    bp = cp - ap;

    costf = sum(abs(   ap(NIR) - psiT(NIR).*delT - ((ap(ref)-psiT(ref).*delT)./bp(ref)).*bp(NIR)    ));

return
    