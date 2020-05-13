function model = updateModelInfo(model)
 length = [0.4140    0.4388    0.0939];
            mass = [ 7.2288   18.5238   10.6926];
            
            com{1} = [ 0.1606;     0.0583;          0];
            com{2} = [ 0.1764;     0.0683;          0];
            com{3} = [ 0.0223;     0.0141;          0];
            
            
            inert{1} = [    1.2387         0         0
                0         0         0
                0         0    1.2387];
            
            inert{2} = [    3.5671         0         0
                0         0         0
                0         0    3.5671];
            
            inert{3} = [  0.0943         0         0
                0         0         0
                0         0    0.0943];
            
            model.transforms(4).t(1:3, 4) = [length(1) 0 0]';
            model.transforms(6).t(1:3, 4) = [length(2) 0 0]';
            model.transforms(8).t(1:3, 4) = [length(3) 0 0]';
            
            for i = 1:3
                model.bodies(i).com = com{i};
                model.bodies(i).I   = inert{i};
                model.bodies(i).m   = mass(i);
            end
end