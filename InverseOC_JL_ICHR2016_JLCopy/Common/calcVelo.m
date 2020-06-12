function velo = calcVelo(pos)
    velo = zeros(size(pos));
    for i = 1:size(pos, 2)
        veloTemp = [0; diff(pos(:, i))];
        velo(:, i) = veloTemp;
    end