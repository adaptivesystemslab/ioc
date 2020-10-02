classdef FbxModel < handle
    %FBXMODEL Loads an ASCII FBX File
    %
    %   This class is used to load an FBX ASCII format file in order to
    %   gain a representation (Nodal) required for exporting into other
    %   formats.
    %   
    %   It is normally used in tandem with the AnimatrikModel to export
    %   into the XML format used in the Animation Lab. The general usage is
    %   as follows:
    %
    %   % Load the Model
    %   fbxObj = FbxModel('ModelFile_ASCII.fbx');
    %   % Get some submodel by name
    %   armando = AnimatrikModel( fbxObj, 'Armando' );
    %   % Save as XML
    %   armando.ExportToXml( 'Models\armando.xml' );
    %
    %   See also: ANIMATRIKMODEL, FBXNODE
    
    properties( Constant )
        PreScale = 0.01; % cm -> m
    end
    properties( SetAccess = private, GetAccess = public )
        RootNode;
    end
    methods
        function obj = FbxModel( fbxFilePath )
            %FBXMODEL Loads FBX ASCII model specified by file path
            
            obj.RootNode = FbxModel.ReadFbxFile( fbxFilePath );
            obj.RootNode.BakeScaling( ones([1,3]) );
        end
        
        function node = GetModel( obj, modelName )
            %GETMODEL Searches through FbxNode Hierarchy to find a
            %particular node with the given name and returns that node.
            
            node = FbxModel.FindNodeInTree( modelName, obj.RootNode );
        end
    end
    
    methods( Static, Access = private )
        function nodeResult = FindNodeInTree( nodeName, node )
            nodeResult = FbxNode.Null;
            
            if strcmp( nodeName, node.Name )
                nodeResult = node;
                return
            end
            
            % Current node isn't the right node...recurse!
            for ni = 1:length(node.Children)
                no = FbxModel.FindNodeInTree( nodeName, node.Children{ni} );
                if no ~= FbxNode.Null
                    nodeResult = no;
                    break
                end
            end
        end
        
        function RootNode = ReadFbxFile( file )
            fhandle = fopen( file );
            nodes = { FbxNode.Root };
            
            NodeIDMap = containers.Map( 'KeyType', 'int32', 'ValueType', 'any' );
            NodeIDMap( nodes{1}.ID ) = nodes{1};
            
            % Read the whole file, such that each cell represents a line in
            % the ASCII file.
            C = textscan( fhandle, '%s', 'Delimiter', '');
            fclose( fhandle );
            C = C{1};
            
            percent = 0;
            waitH = waitbar( percent, 'Loading FBX ASCII Model...' );
            lineCount = numel(C);
            for lineNum = 1:lineCount
                tline = C{lineNum};
                
                newPercent = lineNum / lineCount;
                if( abs(newPercent - percent) >= 0.01 )
                    percent = newPercent;
                    waitbar( percent, waitH );
                end
                
                % Empty line, ignore
                if isempty(tline)
                    continue;
                end

                % Attempt to match with some property
                indexOfModel = regexpi( tline, '^model[:](\s+)' );
                indexOfConnection = regexpi( tline, '^C[:](\s+)' );
                [propStart, propEnd] = regexpi( tline, '^P[:](\s+)"' );

                if length(indexOfModel) == 1
                    % Model property, this represents a node in the FBX
                    % hierarchy of nodes
                    
                    [startId, endId] = regexpi( tline, '[0-9]+' );
                    [startName, endName] = regexpi( tline, '"Model::(.+?)"' );

                    nodeID = str2num(strcat('uint64(', tline(startId:endId), ')'));
                    nodeName = tline((startName(1)+8):(endName(1)-1));

                    nodes{ end + 1 } = FbxNode( nodeID, nodeName );
                    NodeIDMap( nodeID ) = nodes{end};

                elseif length(propStart) == 1 && ~isempty(nodes)
                    
                    % This sets a particular property of a model node; we
                    % only care about transformations
                    tline = tline((propEnd+1):end);
                    quotes = strfind( tline, '"' );
                    tname = tline(1:(quotes(1) - 1));

                    isPreRot = strcmp(tname, 'PreRotation');
                    isTrans = strcmp(tname, 'Lcl Translation');
                    isRot = strcmp(tname, 'Lcl Rotation');
                    isScale = strcmp(tname, 'Lcl Scaling');

                    if isTrans || isRot || isPreRot || isScale
                        delim = strfind(tline, ',');
                        tline = tline((delim(length(delim)-2)+1):end);
                        vector = str2num(tline);
                        if isTrans
                            nodes{ end }.Translation = FbxModel.PreScale .* vector;
                        elseif isRot
                            nodes{ end }.Rotation = vector;
                        elseif isPreRot
                            nodes{ end }.PreRotation = vector;
                        elseif isScale
                            nodes{ end }.Scaling = vector;
                        end
                    end

                elseif length(indexOfConnection) == 1

                    % this property (connection) tells us how nodes are
                    % connected, and allows us to construct the tree of FBX
                    % Nodes
                    [startId, endId] = regexpi( tline, '[0-9]+' );
                    id1 = str2num(strcat('uint64(', tline(startId(1):endId(1)), ')'));
                    id2 = str2num(strcat('uint64(', tline(startId(2):endId(2)), ')'));

                    if( NodeIDMap.isKey( id1 ) && NodeIDMap.isKey( id2 ) )
                        child = NodeIDMap( id1 );
                        parent = NodeIDMap( id2 );
                        parent.AddChild( child );
                    end
                end
            end
            close( waitH );
            
            % return the root node, user can traverse the tree structure
            % using this reference
            RootNode = nodes{ 1 };
        end
    end
end