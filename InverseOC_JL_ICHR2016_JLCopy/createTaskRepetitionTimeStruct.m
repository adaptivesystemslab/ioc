function segmentsOfInterest = createTaskRepetitionTimeStruct() 
%Added by PC. Defines struct with time information of segments of interest %in single arm sequences     
    segmentsOfInterest = struct;     

% Values specified using original sampling rate of 200Hz     
    segmentsOfInterest(1).name = "Subject01_SingleArmBaseLine";     
    segmentsOfInterest(1).start = 19140; % 1 sec before picking takes place
    segmentsOfInterest(1).pickTarget = 19340;     
    segmentsOfInterest(1).startDisplacement = 19380;
    segmentsOfInterest(1).endDisplacement = 19550;
    segmentsOfInterest(1).placementTarget = 19650;    
    segmentsOfInterest(1).end = 19880; % 1 sec after placement ended
    
    segmentsOfInterest(2).name = "Subject01_SingleArm80Time";     
    segmentsOfInterest(2).start = 9230; % 1 sec before picking takes place
    segmentsOfInterest(2).pickTarget = 9430;     
    segmentsOfInterest(2).startDisplacement = 9480;
    segmentsOfInterest(2).endDisplacement = 9640;
    segmentsOfInterest(2).placementTarget = 9700;    
    segmentsOfInterest(2).end = 10070; % 1 sec after placement ended
    
    segmentsOfInterest(3).name = "Subject01_SingleArm60Time";     
    segmentsOfInterest(3).start = 8440; % 1 sec before picking takes place
    segmentsOfInterest(3).pickTarget = 8460;     
    segmentsOfInterest(3).startDisplacement = 8690;
    segmentsOfInterest(3).endDisplacement = 8820;
    segmentsOfInterest(3).placementTarget = 8890;    
    segmentsOfInterest(3).end = 9280; % 1 sec after placement ended
    
    segmentsOfInterest(4).name = "Subject01_SingleArmMinimumTime";     
    segmentsOfInterest(4).start = 8940; % 1 sec before picking takes place
    segmentsOfInterest(4).pickTarget = 9140;     
    segmentsOfInterest(4).startDisplacement = 9180;
    segmentsOfInterest(4).endDisplacement = 9300;
    segmentsOfInterest(4).placementTarget = 9320;    
    segmentsOfInterest(4).end = 9600; % 1 sec after placement ended
    
    segmentsOfInterest(5).name = "Subject02_SingleArmBaseLine";     
    segmentsOfInterest(5).start = 20200; % 1 sec before picking takes place
    segmentsOfInterest(5).pickTarget = 20400;     
    segmentsOfInterest(5).startDisplacement = 20500;
    segmentsOfInterest(5).endDisplacement = 20670;
    segmentsOfInterest(5).placementTarget = 20700;    
    segmentsOfInterest(5).end = 21250; % 1 sec after placement ended

    segmentsOfInterest(6).name = "Subject02_SingleArm80Time";     
    segmentsOfInterest(6).start = 9750; % 1 sec before picking takes place
    segmentsOfInterest(6).pickTarget = 9950;     
    segmentsOfInterest(6).startDisplacement = 10600;
    segmentsOfInterest(6).endDisplacement = 10320;
    segmentsOfInterest(6).placementTarget = 10400;    
    segmentsOfInterest(6).end = 10750; % 1 sec after placement ended
    
    segmentsOfInterest(7).name = "Subject02_SingleArm60Time";     
    segmentsOfInterest(7).start = 6680; % 1 sec before picking takes place
    segmentsOfInterest(7).pickTarget = 6880;     
    segmentsOfInterest(7).startDisplacement = 6940;
    segmentsOfInterest(7).endDisplacement = 7090;
    segmentsOfInterest(7).placementTarget = 7100;    
    segmentsOfInterest(7).end = 7420; % 1 sec after placement ended
    
    segmentsOfInterest(8).name = "Subject02_SingleArmMinimumTime";     
    segmentsOfInterest(8).start = 6030; % 1 sec before picking takes place
    segmentsOfInterest(8).pickTarget = 6230;     
    segmentsOfInterest(8).startDisplacement = 6270;
    segmentsOfInterest(8).endDisplacement = 6350;
    segmentsOfInterest(8).placementTarget = 6350;    
    segmentsOfInterest(8).end = 6660; % 1 sec after placement ended
    
    segmentsOfInterest(9).name = "Subject03_SingleArmBaseLine";     
    segmentsOfInterest(9).start = 10500; % 1 sec before picking takes place
    segmentsOfInterest(9).pickTarget = 10700;     
    segmentsOfInterest(9).startDisplacement = 10790;
    segmentsOfInterest(9).endDisplacement = 10930;
    segmentsOfInterest(9).placementTarget = 11000;    
    segmentsOfInterest(9).end = 11300; % 1 sec after placement ended

    segmentsOfInterest(10).name = "Subject03_SingleArm80Time";     
    segmentsOfInterest(10).start = 8200; % 1 sec before picking takes place
    segmentsOfInterest(10).pickTarget = 8400;     
    segmentsOfInterest(10).startDisplacement = 8470;
    segmentsOfInterest(10).endDisplacement = 8600;
    segmentsOfInterest(10).placementTarget = 8640;    
    segmentsOfInterest(10).end = 9070; % 1 sec after placement ended
    
    segmentsOfInterest(11).name = "Subject03_SingleArm60Time";     
    segmentsOfInterest(11).start = 5460; % 1 sec before picking takes place
    segmentsOfInterest(11).pickTarget = 5660;     
    segmentsOfInterest(11).startDisplacement =5720;
    segmentsOfInterest(11).endDisplacement = 5820;
    segmentsOfInterest(11).placementTarget = 5840;    
    segmentsOfInterest(11).end = 6120; % 1 sec after placement ended
    
    segmentsOfInterest(12).name = "Subject03_SingleArmMinimumTime";     
    segmentsOfInterest(12).start = 4720; % 1 sec before picking takes place
    segmentsOfInterest(12).pickTarget = 4920;     
    segmentsOfInterest(12).startDisplacement = 4990;
    segmentsOfInterest(12).endDisplacement = 5150;
    segmentsOfInterest(12).placementTarget = 5190;    
    segmentsOfInterest(12).end = 5480; % 1 sec after placement ended
    
    segmentsOfInterest(13).name = "Subject04_SingleArmBaseLine";     
    segmentsOfInterest(13).start = 22000; % 1 sec before picking takes place
    segmentsOfInterest(13).pickTarget = 22200;     
    segmentsOfInterest(13).startDisplacement = 22270;
    segmentsOfInterest(13).endDisplacement = 22420;
    segmentsOfInterest(13).placementTarget = 22500;    
    segmentsOfInterest(13).end = 22770; % 1 sec after placement ended
    
    segmentsOfInterest(14).name = "Subject04_SingleArm80Time";     
    segmentsOfInterest(14).start = 16720; % 1 sec before picking takes place
    segmentsOfInterest(14).pickTarget = 16920;     
    segmentsOfInterest(14).startDisplacement = 17020;
    segmentsOfInterest(14).endDisplacement = 17130;
    segmentsOfInterest(14).placementTarget = 17160;    
    segmentsOfInterest(14).end = 17740; % 1 sec after placement ended
    
    segmentsOfInterest(15).name = "Subject04_SingleArm60Time";     
    segmentsOfInterest(15).start = 8900; % 1 sec before picking takes place
    segmentsOfInterest(15).pickTarget = 9100;     
    segmentsOfInterest(15).startDisplacement = 9210;
    segmentsOfInterest(15).endDisplacement = 9340;
    segmentsOfInterest(15).placementTarget = 9380;    
    segmentsOfInterest(15).end = 9640; % 1 sec after placement ended
    
    segmentsOfInterest(16).name = "Subject04_SingleArmMinimumTime";     
    segmentsOfInterest(16).start = 8640; % 1 sec before picking takes place
    segmentsOfInterest(16).pickTarget = 8840;     
    segmentsOfInterest(16).startDisplacement = 8970;
    segmentsOfInterest(16).endDisplacement = 9060;
    segmentsOfInterest(16).placementTarget = 9100;    
    segmentsOfInterest(16).end = 9410; % 1 sec after placement ended
    
% This is the original segment that was analyzed. However, since results
% are confusing, new repetitions (see below) are considered instead
%     segmentsOfInterest(5).name = "Subject05_SingleArmBaseLine";     
%     segmentsOfInterest(5).start = 23630; % 1 sec before picking takes place
%     segmentsOfInterest(5).pickTarget = 23800;     
%     segmentsOfInterest(5).startDisplacement = 23930;
%     segmentsOfInterest(5).endDisplacement = 24035;
%     segmentsOfInterest(5).placementTarget = 24130;    
%     segmentsOfInterest(5).end = 25000; % 1 sec after placement ended

% Option 1. It is included in last 30 sec of full sequence
    segmentsOfInterest(17).name = "Subject05_SingleArmBaseLine";     
    segmentsOfInterest(17).start = 22140; % 1 sec before picking takes place
    segmentsOfInterest(17).pickTarget = 22340;     
    segmentsOfInterest(17).startDisplacement = 22420;
    segmentsOfInterest(17).endDisplacement = 22600;
    segmentsOfInterest(17).placementTarget = 22660;    
    segmentsOfInterest(17).end = 23180; % 1 sec after placement ended
    
    segmentsOfInterest(18).name = "Subject05_SingleArm80Time";     
    segmentsOfInterest(18).start = 15270; % 1 sec before picking takes place
    segmentsOfInterest(18).pickTarget = 15470;     
    segmentsOfInterest(18).startDisplacement = 15520;
    segmentsOfInterest(18).endDisplacement = 15680;
    segmentsOfInterest(18).placementTarget = 15750;    
    segmentsOfInterest(18).end = 16200; % 1 sec after placement ended
    
    segmentsOfInterest(19).name = "Subject05_SingleArm60Time";     
    segmentsOfInterest(19).start = 11200; % 1 sec before picking takes place
    segmentsOfInterest(19).pickTarget = 11400;     
    segmentsOfInterest(19).startDisplacement = 11520;
    segmentsOfInterest(19).endDisplacement = 11600;
    segmentsOfInterest(19).placementTarget = 11700;    
    segmentsOfInterest(19).end = 12080; % 1 sec after placement ended
    
    segmentsOfInterest(20).name = "Subject05_SingleArmMinimumTime";     
    segmentsOfInterest(20).start = 11760; % 1 sec before picking takes place
    segmentsOfInterest(20).pickTarget = 11960;     
    segmentsOfInterest(20).startDisplacement = 12050;
    segmentsOfInterest(20).endDisplacement = 12210;
    segmentsOfInterest(20).placementTarget = 12300;    
    segmentsOfInterest(20).end = 12610; % 1 sec after placement ended

% Option 2. Not included in last 30 sec of full sequence    
%     segmentsOfInterest(5).name = "Subject05_SingleArmBaseLine";     
%     segmentsOfInterest(5).start = 15230; % 1 sec before picking takes place
%     segmentsOfInterest(5).pickTarget = 15520;     
%     segmentsOfInterest(5).startDisplacement = 15590;
%     segmentsOfInterest(5).endDisplacement = 15730;
%     segmentsOfInterest(5).placementTarget = 15750;    
%     segmentsOfInterest(5).end = 16010; % 1 sec after placement ended
    
end