classdef FbxNode < handle
    properties( Constant )
        Null = FbxNode( -1, 'Empty' );
        Root = FbxNode( 0, 'Root' );
    end
    properties
        Name;
        ID;
        Translation = zeros([1,3]);
        Rotation = zeros([1,3]);
        PreRotation = zeros([1,3]);
        Scaling = ones([1,3]);
        Children = {};
        AnimationValues = [];
        AnimationTimes = [];
    end
    methods
        function [] = BakeScaling( obj, scale )
            obj.Translation = obj.Translation .* scale;
            
            scale = scale .* obj.Scaling;
            for c = 1:length(obj.Children)
                child = obj.Children{ c };
                child.BakeScaling( scale );
            end
        end
        function obj = FbxNode( id, name )
            obj.ID = id;
            obj.Name = name;
        end
        function [] = AddChild( obj, c )
            obj.Children{ length(obj.Children) + 1 } = c;
        end
    end
end