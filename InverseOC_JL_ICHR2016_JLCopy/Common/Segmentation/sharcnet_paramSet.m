function package = sharcnet_paramSet(package, param)
    
    package.manSeg = param.manSeg;
    switch param.manSeg
        case 'Segmentation_manual'
            package.manSeg = 'Segmentation_manual';
            package.halfSegments = 1;

        case 'Segmentation_manual_annotatedZVC'
            package.manSeg = 'Segmentation_manual_annotatedZVC';
            package.halfSegments = 1;

        case 'Segmentation_manual_annotatedZVCHalfSeg'
            package.manSeg = 'Segmentation_manual_annotatedZVCHalfSeg';
            package.halfSegments = 0;
    end