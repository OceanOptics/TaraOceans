function [msec,c,a,cwl,awl] = Read_ACSBin043_0909(thisfile)

    disp(['Reading ' thisfile ' using Read_ACSBin...']);
    
    err = unix(['wine PREPACS.exe acs043_0909.dev ' thisfile ' prepacs.tmp +nul']);	
	
    
	[dat] = textread('prepacs.tmp','','delimiter','\t', 'headerlines',1);

	if ~isempty(dat)
		msec = dat(:,1);
		c = dat(:,2:89);
		a = dat(:,90:177);
	else
		[msec,a,c] = deal([])
    end
		
    % s/n 043_0909
    cwl = [400.7	403.9	407.6	410.9	414.4	418.2	421.9	426.0	429.9	433.6	437.5	441.2	445.3	449.3	453.6	457.8	461.7	465.9	470.1	474.3	478.5	482.9	487.1	491.3	495.3	499.5	503.6	507.8	512.1	516.3	520.4	524.6	528.7	532.5	536.5	540.3	544.4	548.1	552.1	555.8	559.6	563.2	566.8	570.5	574.2	577.8	581.5	585.4	589.3	593.2	597.2	601.1	605.3	609.3	613.4	617.5	621.5	625.5	629.5	633.5	637.7	641.8	645.8	649.9	654.0	658.0	662.3	666.4	670.3	674.2	678.1	681.9	685.6	689.5	693.1	696.8	700.5	703.8	707.3	710.7	714.0	717.3	720.8	724.0	727.2	730.4	733.2	736.4];
    awl = [401.4	404.6	408.7	411.7	415.3	419.0	423.0	426.9	430.6	434.5	438.4	442.3	446.3	450.4	454.5	458.5	462.6	466.8	471.0	475.2	479.4	483.6	487.8	491.8	496.0	500.2	504.3	508.7	512.8	517.2	521.3	525.3	529.2	533.2	537.2	541.1	545.1	548.8	552.8	556.5	560.1	562.8	566.4	570.2	573.9	577.6	581.2	585.1	589.0	592.7	596.6	600.8	605.0	609.0	612.9	617.0	621.2	625.2	629.2	633.2	637.3	641.3	645.3	649.4	653.4	657.7	661.6	665.5	669.6	673.6	677.5	681.4	685.1	689.0	692.6	696.2	699.9	703.4	706.8	710.2	713.7	716.9	720.3	723.3	726.8	729.7	732.9	735.8];
    
    
	%[msec,c,a,cwl,awl] = Read_ACSDat043('prepacs.tmp');
		
return