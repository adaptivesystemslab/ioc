function [height, weight, gender] = patientTypeSelData_Healthy1(subjNumber)
    % subject information for the lowerbody exercise set
    
    % height in [m]
    % weight in [kg]
    % gender in 'm' or 'f'
    
    switch subjNumber
        case {1, 21, 22, 23}
            height = 1.7;
            weight = 65.7;
            gender = 'm';
            
        case 2
            height = 1.55;
            weight = 70;
            gender = 'f';
            
        case 3
            height = 1.6;
            weight = 47.72;
            gender = 'f';

        case 4
            height = 1.83;
            weight = 90;
            gender = 'm';
            
        case 5
            height = 1.5;
            weight = 52;
            gender = 'f';
            
        case 6
            height = 1.7;
            weight = 70;
            gender = 'm';
            
        case 7
            height = 1.68;
            weight = 71;
            gender = 'm';
            
        case 8 
            height = 1.75;
            weight = 60;
            gender = 'm';
            
        case 9
            height = 1.575;
            weight = 52;
            gender = 'f';
            
        case 10
            height = 1.77;
            weight = 63.5;
            gender = 'm';
            
        case 11
            height = 1.7;
            weight = 61;
            gender = 'm';
            
        case 12
            height = 1.63;
            weight = 60.7;
            gender = 'f';
            
        case 13
            height = 1.8;
            weight = 70;
            gender = 'm';
            
        case 14
            height = 1.76;
            weight = 68;
            gender = 'm';
            
        case 15
            height = 1.64;
            weight = 63;
            gender = 'f';
            
        case 16
            height = 1.57;
            weight = 55;
            gender = 'f';
            
        case 17
            height = 1.7;
            weight = 82;
            gender = 'm';
            
        case 18
            height = 1.68;
            weight = 56;
            gender = 'f';
            
        case 19
            height = 1.83;
            weight = 71;
            gender = 'm';

        case 20
            height = 1.77;
            weight = 75;
            gender = 'm';            
            
        otherwise
            
    end
end