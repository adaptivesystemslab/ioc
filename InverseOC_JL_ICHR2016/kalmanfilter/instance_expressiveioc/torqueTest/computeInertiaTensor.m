function I = computeInertiaTensor(mass, lenght, width, height, com)
    
    % Within this function we assume all links are represented as simple 3D rectangles
    % The rectangle is aligned with the y-axis, i.e, lenght, width and height are measured
    % along y-axis, x-axis and z-axis respectively.

    % Check whether all parameters have been provided. If not, set default values
    if ~exist('mass', 'var')
        mass = 1;
    end
    
    if ~exist('lenght', 'var')
        lenght = 1e-2;
    end
    
    if ~exist('width', 'var')
        width = 1e-2;
    end
    
    if ~exist('height', 'var')
        height = 1e-2;
    end
        
    Icom = zeros(3,3);
    Icom(1,1) = 1/12*mass*(lenght*lenght + height*height);
    Icom(2,2) = 1/12*mass*(width*width + height*height);
    Icom(3,3) = 1/12*mass*(lenght*lenght + width*width);
    
    I = Icom;
    
    if nargin>4
       S = [0 -com(3) com(2); com(3) 0 -com(1); -com(2) com(1) 0];
       
       Iframe = Icom + mass*(S'*S);       
       I = Iframe;
    end
end
