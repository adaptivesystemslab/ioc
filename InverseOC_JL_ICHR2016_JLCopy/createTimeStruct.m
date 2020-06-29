function segmentsOfInterest = createTimeStruct() 
%Added by PC. Defines struct with time information of segments of interest %in single arm sequences     
segmentsOfInterest = struct;     
segmentsOfInterest(1).name = "Subject01_SingleArm60Time";     
segmentsOfInterest(1).middle.start =2600;     
segmentsOfInterest(1).middle.end =2935;     
segmentsOfInterest(1).final.start =4320;    
segmentsOfInterest(1).final.middle = 4503; 
segmentsOfInterest(1).final.end =4665;     

segmentsOfInterest(2).name = "Subject01_SingleArm80Time";     
segmentsOfInterest(2).middle.start =3360;     
segmentsOfInterest(2).middle.end =3650;     
segmentsOfInterest(2).final.start =4655;     
segmentsOfInterest(2).final.middle = 4725; 
segmentsOfInterest(2).final.end = 4850;    

segmentsOfInterest(3).name = "Subject01_SingleArmBaseLine";     
segmentsOfInterest(3).middle.start = 5045;     
segmentsOfInterest(3).middle.end = 5510;     
segmentsOfInterest(3).final.start = 9520;    
segmentsOfInterest(3).final.middle = 10000; 
segmentsOfInterest(3).final.end = 10345;     

segmentsOfInterest(4).name = "Subject01_SingleArmMinimumTime";     
segmentsOfInterest(4).middle.start =1650;    
segmentsOfInterest(4).middle.end =1900;     
segmentsOfInterest(4).final.start =3800;     
segmentsOfInterest(4).final.middle = 3926; 
segmentsOfInterest(4).final.end =4040;    

segmentsOfInterest(5).name = "Subject02_SingleArm60Time";     
segmentsOfInterest(5).middle.start =2450;    
segmentsOfInterest(5).middle.end =2840;     
segmentsOfInterest(5).final.start =4950;     
segmentsOfInterest(5).final.middle = 5075; 
segmentsOfInterest(5).final.end = 5225;     

segmentsOfInterest(6).name = "Subject02_SingleArm80Time";     
segmentsOfInterest(6).middle.start =3100;     
segmentsOfInterest(6).middle.end =3515;     
segmentsOfInterest(6).final.start = 5735;  
segmentsOfInterest(6).final.middle = 6000; 
segmentsOfInterest(6).final.end = 6170;     

segmentsOfInterest(7).name = "Subject02_SingleArmBaseLine";     
segmentsOfInterest(7).middle.start = 5375;     
segmentsOfInterest(7).middle.end = 6175;     
segmentsOfInterest(7).final.start = 10200;   
segmentsOfInterest(7).final.middle = 10430; 
segmentsOfInterest(7).final.end = 10690;     

segmentsOfInterest(8).name = "Subject02_SingleArmMinimumTime";     
segmentsOfInterest(8).middle.start =2275;     
segmentsOfInterest(8).middle.end =2570;    
segmentsOfInterest(8).final.start =4195;   
segmentsOfInterest(8).final.middle = 4405; 
segmentsOfInterest(8).final.end =4550;     

segmentsOfInterest(9).name = "Subject03_SingleArm60Time";     
segmentsOfInterest(9).middle.start =1950;     
segmentsOfInterest(9).middle.end =2235;     
segmentsOfInterest(9).final.start =3750;  
segmentsOfInterest(9).final.middle = 3839; 
segmentsOfInterest(9).final.end =3950;     

segmentsOfInterest(10).name = "Subject03_SingleArm80Time";     
segmentsOfInterest(10).middle.start =1975;     
segmentsOfInterest(10).middle.end =2175;     
segmentsOfInterest(10).final.start =4765;     
segmentsOfInterest(10).final.middle = 4920; 
segmentsOfInterest(10).final.end = 5045;     

segmentsOfInterest(11).name = "Subject03_SingleArmBaseLine";     
segmentsOfInterest(11).middle.start =3035;     
segmentsOfInterest(11).middle.end =3300;     
segmentsOfInterest(11).final.start = 5050;    
segmentsOfInterest(11).final.middle = 5245;  
segmentsOfInterest(11).final.end = 5375;     

segmentsOfInterest(12).name = "Subject03_SingleArmMinimumTime";     
segmentsOfInterest(12).middle.start =2120;     
segmentsOfInterest(12).middle.end =2400;     
segmentsOfInterest(12).final.start =3305;     
segmentsOfInterest(12).final.middle = 3404; 
segmentsOfInterest(12).final.end =3500;     

segmentsOfInterest(13).name = "Subject04_SingleArm60Time";     
segmentsOfInterest(13).middle.start =2510;     
segmentsOfInterest(13).middle.end =2890;     
segmentsOfInterest(13).final.start = 5075;    
segmentsOfInterest(13).final.middle = 5243; 
segmentsOfInterest(13).final.end = 5370;     

segmentsOfInterest(14).name = "Subject04_SingleArm80Time";     
segmentsOfInterest(14).middle.start =3150;     
segmentsOfInterest(14).middle.end =3615;     
segmentsOfInterest(14).final.start = 6975;     
segmentsOfInterest(14).final.middle = 7355; 
segmentsOfInterest(14).final.end = 7890;     

segmentsOfInterest(15).name = "Subject04_SingleArmBaseLine";     
segmentsOfInterest(15).middle.start = 7705;     
segmentsOfInterest(15).middle.end = 8710;     
segmentsOfInterest(15).final.start = 13245;     
segmentsOfInterest(15).final.middle = 13600; 
segmentsOfInterest(15).final.end = 13830;     

segmentsOfInterest(16).name = "Subject04_SingleArmMinimumTime";     
segmentsOfInterest(16).middle.start =2665;     
segmentsOfInterest(16).middle.end =3225;     
segmentsOfInterest(16).final.start =4865; 
segmentsOfInterest(16).final.middle = 4966; 
segmentsOfInterest(16).final.end = 5100;     

segmentsOfInterest(17).name = "Subject05_SingleArm60Time";     
segmentsOfInterest(17).middle.start =3300;     
segmentsOfInterest(17).middle.end =3775;     
segmentsOfInterest(17).final.start = 6000;    
segmentsOfInterest(17).final.middle = 6190; 
segmentsOfInterest(17).final.end = 6400;     

segmentsOfInterest(18).name = "Subject05_SingleArm80Time";     
segmentsOfInterest(18).middle.start =3370;     
segmentsOfInterest(18).middle.end =3800;     
segmentsOfInterest(18).final.start = 6600;     
segmentsOfInterest(18).final.middle = 6941; 
segmentsOfInterest(18).final.end = 7370;     

segmentsOfInterest(19).name = "Subject05_SingleArmBaseLine";     
segmentsOfInterest(19).middle.start = 5150;     
segmentsOfInterest(19).middle.end = 6450;     
segmentsOfInterest(19).final.start = 11925;     
segmentsOfInterest(19).final.middle = 12200; 
segmentsOfInterest(19).final.end = 12490;     

segmentsOfInterest(20).name = "Subject05_SingleArmMinimumTime";     
segmentsOfInterest(20).middle.start =2965;     
segmentsOfInterest(20).middle.end =3465;     
segmentsOfInterest(20).final.start = 6000;     
segmentsOfInterest(20).final.middle = 6156;     
segmentsOfInterest(20).final.end = 6315; 
end