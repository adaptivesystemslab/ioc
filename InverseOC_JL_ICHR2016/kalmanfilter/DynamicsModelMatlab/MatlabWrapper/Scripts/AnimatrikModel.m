classdef AnimatrikModel
    %ANIMATRIKMODEL 
    %   Takes a FbxModel, prunes it to a particular sub-model that you
    %   provide by name and then allows you to export it to XML format.
    %   Note that when you specify a 'ModelName', it appends ':Solving'
    %   to search for 'ModelName:Solving' as a node in the model.
    %
    %   If you don't want to search for 'ModelName:Solving' as the root
    %   node you can instead pass another string to be appended to the
    %   Model Name to search for 'ModelName(param3)'
    %
    %   See also: FBXMODEL, FBXNODE, EXPORTTOXML
    
    properties( Constant )
        RootPostfix = ':Solving';
    end
    properties( Access = private )
        mRootPrefix = 'Armando';
        mRootNode;
        mCModel;
    end
    methods
        function obj = AnimatrikModel( inobj, modelName, postFix )
            if nargin < 3
                postFix = obj.RootPostfix;
            end
            obj.mRootPrefix = modelName;
            obj.mRootNode = inobj.GetModel( strcat(obj.mRootPrefix, postFix) );
        end
        function ExportToXml( obj, outFile )
            fprintf( '\nPreparing to export model to Xml.\n' );
            c = AnimatrikModel.CountNodes( obj.mRootNode );
            fbxToXml( obj.mRootPrefix, obj.mRootNode, outFile, c );
            fprintf( 'Export to Xml completed!\n' );
        end
    end
    methods( Static, Access = private )
        function c = CountNodes( node )
            c = 1;
            for i = 1:length(node.Children)
                c = c + AnimatrikModel.CountNodes( node.Children{i} );
            end
        end
    end
end

function fbxToXml( modelName, node, outFile, count )
    %FBXTOXML Called to export a particular FbxNode rootnode as XML.
    %   Will call FRAMESTOXML to output children node. This function also
    %   writes all the necessary root node XML content, such as Gravity and
    %   World coordinate frame
    
    docNode = com.mathworks.xml.XMLUtils.createDocument('rlmdl');
    docRoot = docNode.getDocumentElement;

    modelElement = docNode.createElement('model');
    docRoot.appendChild( modelElement );

    manufacturer = docNode.createElement('manufacturer');
    manufacturer.setTextContent('Animatrik');
    modelElement.appendChild( manufacturer );

    name = docNode.createElement('name');
    name.setTextContent( modelName );
    modelElement.appendChild( name );

    world = docNode.createElement('world');
    world.setAttribute('id', 'world');

    transformToXml( docNode, world, node.Translation, [0;0;0] );  
    
    gravity = docNode.createElement('g');
    vectorToXml( docNode, gravity, [0;0;-9.80665] );
    
    modelElement.appendChild( world );

    % Recurse through model tree and output XML
    fprintf( 'Joints remaining:\t%8u', count );
    count = framesToXml( docNode, modelElement, world, node.Children{1}, count );

    % Return number of joints remaining
    count = count - 1;
    
    % Print out number of joints remaining (should be 0)
    for i = 1:8
        fprintf( '\b' );
    end
    fprintf( '%8u\n', count );
    
    fprintf( 'Writing out Xml File.\n' );
    xmlwrite( outFile, docNode );
end

function vectorToXml( doc, element, v )
% VECTORTOXML Prints Vector into XML Element of Document
% Inserts the 3d vector v into XML element like so:
%     <element>
%         <x>v(1)</x>
%         <y>v(2)</y>
%         <z>v(3)</z>
%     ...

    xn = doc.createElement('x');
    yn = doc.createElement('y');
    zn = doc.createElement('z');
    
    xn.setTextContent( num2str(v(1),16) );
    yn.setTextContent( num2str(v(2),16) );
    zn.setTextContent( num2str(v(3),16) );

    element.appendChild(xn);
    element.appendChild(yn);
    element.appendChild(zn);
end

function transformToXml( doc, element, transl, rot )
% TRANSFORMTOXML Prints Transform into XML Element of Document
% Calls VECTORTOXML to print both the rotation and translation provided
% as 3D Vectors.

    rotation = doc.createElement('rotation');
    translation = doc.createElement('translation');

    vectorToXml( doc, rotation, rot );
    vectorToXml( doc, translation, transl );

    element.appendChild( rotation );
    element.appendChild( translation );
end

function writeGeneralFrame( doc, parent, frameType, aId )
    preAid = strcat( 'pre:', aId );
    postAid = strcat( 'post:', aId );

    writeFrameDefinition( doc, parent, preAid );
    writeFrameDefinition( doc, parent, postAid );

    frame = doc.createElement(frameType);
    frame.setAttribute( 'id', strcat( preAid, '_to_', postAid ) );

    frameLink = doc.createElement('frame');
    frameLinkA = doc.createElement('a');
    frameLinkA.setAttribute( 'idref', preAid );
    frameLinkB = doc.createElement('b');
    frameLinkB.setAttribute( 'idref', postAid );
    frameLink.appendChild(frameLinkA);
    frameLink.appendChild(frameLinkB);
    frame.appendChild(frameLink);

    parent.appendChild(frame);
end

function frame = writeFixedFrame( doc, parent, aId, bId, transl, rot )
    fixedFrame = doc.createElement('fixed');
    fixedFrame.setAttribute( 'id', strcat( aId, '_to_', bId ) );

    fixedFrameLink = doc.createElement('frame');
    fixedFrameLinkA = doc.createElement('a');
    fixedFrameLinkA.setAttribute( 'idref', aId );
    fixedFrameLinkB = doc.createElement('b');
    fixedFrameLinkB.setAttribute( 'idref', bId );
    fixedFrameLink.appendChild(fixedFrameLinkA);
    fixedFrameLink.appendChild(fixedFrameLinkB);
    fixedFrame.appendChild(fixedFrameLink);
    
    transformToXml( doc, fixedFrame, transl, rot );
    
    frame = fixedFrame;
    parent.appendChild(frame);
end

function frame = writeFrameDefinition( doc, parent, frameName )
    frame = doc.createElement('frame');
    frame.setAttribute( 'id', frameName );
    parent.appendChild(frame);
end

function writePrismatic3DOF( doc, parent, preName, name, postName )
    jXid = strcat( name,':prismX' );
    jYid = strcat( name,':prismY' );
    jZid = strcat( name,':prismZ' );

    writeFixedFrame( doc, parent, preName, strcat('pre:', jXid), ...
                                    zeros([1,3]), [0;90;0] );
    writeGeneralFrame( doc, parent, 'prismatic', jXid );


    writeFixedFrame( doc, parent, strcat('post:', jXid), strcat('pre:', jYid), ...
                                    zeros([1,3]), [-90;0;0] );
    writeGeneralFrame( doc, parent, 'prismatic', jYid );


    writeFixedFrame( doc, parent, strcat('post:', jYid), strcat('pre:', jZid), ...
                                    zeros([1,3]), [90;0;-90] );
    writeGeneralFrame( doc, parent, 'prismatic', jZid );


    writeFixedFrame( doc, parent, strcat('post:',jZid), postName, zeros([1,3]), zeros([1,3]) );
end

function writeRevolute3DOF( doc, parent, preName, name, postName )
    jXid = strcat( name,':rotateX' );
    jYid = strcat( name,':rotateY' );
    jZid = strcat( name,':rotateZ' );

    writeGeneralFrame( doc, parent, 'revolute', jXid );
    writeGeneralFrame( doc, parent, 'revolute', jYid );
    writeGeneralFrame( doc, parent, 'revolute',jZid );
    
%     writeFixedFrame( doc, parent, preName, strcat('pre:', jXid), ...
%                                     zeros([1,3]), [0;90;0] );
% 
% 
%     writeFixedFrame( doc, parent, strcat('post:', jXid), strcat('pre:', jYid), ...
%                                     zeros([1,3]), [-90;0;0] );
% 
% 
%     writeFixedFrame( doc, parent, strcat('post:', jYid), strcat('pre:', jZid), ...
%                                     zeros([1,3]), [90;0;-90] );
%     writeFixedFrame( doc, parent, strcat('post:',jZid), postName, zeros([1,3]), zeros([1,3]) );

    writeFixedFrame( doc, parent, preName, strcat('pre:', jXid), ...
                                    zeros([1,3]), [0;0;0] );


    writeFixedFrame( doc, parent, strcat('post:', jXid), strcat('pre:', jYid), ...
                                    zeros([1,3]), [-90;0;0] );


    writeFixedFrame( doc, parent, strcat('post:', jYid), strcat('pre:', jZid), ...
                                    zeros([1,3]), [0;-90;0] );
    writeFixedFrame( doc, parent, strcat('post:',jZid), postName, zeros([1,3]), [90,90,0]);



end


function n = framesToXml( doc, parent, parentFrame, fbxnode, nleft )
    %FRAMESTOXML Recursively generates XML for a particular frame.
    %   Generates the Pre/Post Frames, Children Frames (recursively) and
    %   the joints of this particular frame.
    
    preMe = strcat( 'pre:', fbxnode.Name );
    postMe = strcat( 'post:', fbxnode.Name );
    writeFrameDefinition( doc, parent, preMe );
    postFrame = writeFrameDefinition( doc, parent, postMe );

    % Generate Children Definitions and Joints first
    for i = 1:length(fbxnode.Children)
        nleft = framesToXml( doc, parent, postFrame, fbxnode.Children{i}, nleft );
    end
    
    pid = char(parentFrame.getAttribute('id'));
    
    % Write fixed frame from parent to me.
    writeFixedFrame( doc, parent, pid, preMe, ...
                     fbxnode.Translation, fbxnode.Rotation );
    disp(['Name:' fbxnode.Name])
    disp('T:')
    disp(fbxnode.Translation)
    disp('R:')
    disp(fbxnode.PreRotation)
    disp(' ');
    % Write joints
    ptag = char(parentFrame.getTagName());
    if( strcmp(ptag, 'world') )
        middleFrameName = strcat('prism2revolute:', fbxnode.Name);
        writeFrameDefinition( doc, parent, middleFrameName );
        writePrismatic3DOF( doc, parent, preMe, fbxnode.Name, middleFrameName );
        writeRevolute3DOF( doc, parent, middleFrameName, fbxnode.Name, postMe );
    elseif( length(fbxnode.Children) >= 1 )
        writeRevolute3DOF( doc, parent, preMe, fbxnode.Name, postMe );
    end
    
    % Return number of joints remaining
    n = nleft - 1;
    
    % Print out number of joints remaining.
    for i = 1:8
        fprintf( '\b' );
    end
    fprintf( '%8u \n', n ); 
end