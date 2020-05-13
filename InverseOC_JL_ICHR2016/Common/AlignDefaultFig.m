% In dual-monitor setups, sometimes the figures are displayed
% offscreen. A recommended fix:            
fprintf('SCS: Changing default window position to laptop screen\n');

% check existing position
% get(0,'DefaultFigurePosition')
% get(0,'FactoryFigurePosition')

% set the new position
CurentPosition = get(0,'DefaultFigurePosition')
NewPosition = [30 280 560 420]
% figure('units','pixels','position',NewPosition)
set(0,'defaultfigureposition',NewPosition)