function [M, Poses] = testMeasurements(vis, mdl, cam)
    I1 = [];
    C1 = [];
    V1 = [];
    ExpWorld = vertcat( [4.5, -0.13888, 2 - 0.43603, 1], ...
                        [4.5, 0.43085, 2 + 0.08427, 1], ...
                        [4.5, 0.43085, 2 - 0.43603, 1], ...
                        [4.5, 0.43085, 2 + 0.43603, 1], ...
                        [4.5, -0.43085, 2 + 0.43603, 1], ...
                        [4.5, -0.43085, 2 - 0.43603, 1]);
    M = NaN(240, 2 * size(ExpWorld, 1));
    Poses = zeros(240, numel(mdl.joints));
    mdl.position = 0;
    for t = 0:240
        mdl.joints(1).position = sin(2 * pi * t / 240) * (pi/32);
        mdl.joints(2).position = cos(2 * pi * t / 53) * (pi/16);
        mdl.joints(3).position = tanh((t - 120)/120) * pi/16;
        mdl.forwardPosition;
        
        cam.updatePose;

        vis.update;

        I0 = I1;
        C0 = C1;
        V0 = V1;
        I1 = cam.observation;
        [~, V1, ~] = cam.getFeatures;
        
        if t == 0
            continue;
        end

        Poses(t,:) = mdl.position;
        %featureIndices = matchFeatures(C0, C1);
        %matchedV0 = V0(featureIndices(:,1),:);
        %matchedV1 = V1(featureIndices(:,2),:);
        %showMatchedFeatures(I0, I1, matchedV0, matchedV1);

%         imshow(I1);
%         hold on;
%         plot(V1(:,1), V1(:,2), 'r.', 'MarkerSize', 1);
%         hold off;
%         drawnow;

        for iv = 1:size(V1, 1)
            dst = Inf;
            is = NaN;
            for ireal = 1:size(ExpWorld, 1)
                P = cam.projectPoint(ExpWorld(ireal,:)')';
                if norm(P - V1(iv,:)) < dst && norm(P - V1(iv,:)) < 5
                    is  = ireal;
                    dst = norm(P - V1(iv,:));
                end
            end
            if(~isnan(is))
                M(t, (2*is-1):2*is) = V1(iv,:);
            end
        end
    end

end