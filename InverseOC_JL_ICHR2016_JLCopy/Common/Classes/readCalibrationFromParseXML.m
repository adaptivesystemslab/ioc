function calibrationStruct = readCalibrationFromParseXML(calibrationNode)
    alignmentRow = eval(['[' calibrationNode(1).Children.Data '];']);
    sensitivityRow = eval(['[' calibrationNode(2).Children.Data '];']);
    offsetRow = eval(['[' calibrationNode(3).Children.Data '];']);

    calibrationStruct.alignment = reshape(alignmentRow, [3 3]);
    calibrationStruct.sensitivity = reshape(sensitivityRow, [3 3]);
    calibrationStruct.offset = reshape(offsetRow, [3 1]); 
end