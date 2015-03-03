function [msec,c,a,cwl,awl] = Read_ACSDat007_0112(fname)

	[dat] = textread(fname,'','delimiter','\t', 'headerlines',21);

	if ~isempty(dat)
		msec = dat(:,1);
		c = dat(:,2:85);
		a = dat(:,86:169);
	else
		[msec,a,c] = deal([])
    end
		
    % s/n 007_0112
    cwl = [399.1 403.7 408.2 412.5 416.9 421.7 426.9 431.5 435.7 440.5 444.9 450.0 454.8 459.4 464.0 468.7 473.4 478.3 483.8 488.1 492.7 497.2 501.9 506.8 511.6 516.6 521.6 526.3 531.0 535.3 539.9 544.4 549.2 553.8 558.4 562.8 567.3 571.2 575.4 579.5 583.9 588.3 592.7 597.2 601.8 606.5 611.2 616.0 620.7 625.2 629.7 634.2 638.7 643.3 647.8 652.7 657.3 661.9 666.5 670.9 675.7 679.9 684.2 688.4 692.3 696.5 700.2 704.1 707.9 711.5 715.2 718.9 722.2 725.9 729.2 732.6 735.8 738.8 741.9 744.9 747.6 750.6 753.1 755.5];	
    awl = [398.0 402.8 407.3 411.9 416.2 420.8 426.0 430.8 434.8 439.4 444.4 448.8 453.9 458.5 462.9 467.6 472.6 477.5 482.5 487.4 491.8 496.3 500.7 505.5 510.4 515.6 520.6 525.4 530.1 534.7 539.2 543.7 548.5 553.1 557.9 562.1 566.6 572.0 576.3 580.5 584.6 588.8 593.5 598.1 602.8 607.3 612.3 617.0 621.7 626.2 630.7 635.3 639.8 644.3 648.9 653.7 658.5 663.1 667.7 672.1 676.5 680.9 685.0 689.0 693.1 697.1 701.0 704.7 708.4 712.3 715.8 719.3 722.9 726.2 729.7 733.0 735.9 739.4 742.3 745.3 747.9 750.5 753.3 756.5];
return

