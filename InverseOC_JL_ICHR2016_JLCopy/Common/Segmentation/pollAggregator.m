function select = pollAggregator(name)

    select.settingName = name;
    parsing = regexp(name, '_', 'split');

    switch parsing{1}
        case 'None'
            select.name = 'None';

        case 'Boosting'
            select.name = 'Boosting';
            select.iteration = str2num(parsing{2});
    %         Resampling method, ensemble count

        case 'Bagging'
            select.name = 'Bagging';
            select.iteration = str2num(parsing{2});
    %         ensemble count
    end
end