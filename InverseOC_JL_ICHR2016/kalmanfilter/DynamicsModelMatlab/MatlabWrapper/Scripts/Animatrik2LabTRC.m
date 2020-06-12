function Animatrik2LabTRC( trcInFile, MarkerRotateList )
    %ANIMATRIK2LABTRC Scales Animatrik TRC file into meters
    %   More importantly, the list passed in as a second parameter will
    %   enforce a -90 degree rotation about x on those markers.
    
    % mm -> m
    TRCPreScale = 0.1 * FbxModel.PreScale;
    
    % Some markers are incorrectly rotated 90degrees in x
    TRCRotFactor = rotx( 90, 'deg' );
    
    % Result FIle
    TrcIn = strcat( trcInFile, '.trc' );
    TrcOut = strcat( trcInFile, '_Corrected.trc' );
    
    % Don't load using read TRC, it will fail because the markers
    % names start with an underscore ('_') !
    text = fileread( TrcIn );
    expr = '(\W+)(_)(\d*)';
    text = regexprep( text, expr, '$1M_$3' );
        fid = fopen( TrcOut, 'w' );
        fwrite( fid, text, 'char*1' );
        fclose( fid );
    clear text;
    
    % Now load the trc and perform corrections
    %  - Assign missing markers
    Take = readTrc( TrcOut );
    names = fieldnames( Take.data );
    for name = { names{:} }
        i = name{1};
        if( size(Take.data.(i),2) > 1 )
            data = TRCPreScale * Take.data.(i);
            
            FixedAlready = strfind( MarkerRotateList, i );
            if( ~isempty( cat( 1, FixedAlready{:} ) ) )
                data = (TRCRotFactor * data')';
            end  
            
            missing = abs( data ) < 1e-8;
            Take.data.(i)(missing) = 1e9;
            Take.data.(i)(~missing) = data(~missing);
        end
    end

    Take.NumFrames = uint32(Take.NumFrames);
    Take.NumMarkers = uint32(Take.NumMarkers);
    Take.OrigDataStartFrame = uint32(Take.OrigDataStartFrame);
    Take.OrigNumFrames = uint32(Take.OrigNumFrames);
    writeTrc(Take, TrcOut);
end