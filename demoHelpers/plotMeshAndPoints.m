function  plotMeshAndPoints( V, F, params )

%--------------------------------------------
% Initialization
%--------------------------------------------
params.null = [];
n = getoptions(params,'n',1);
verInd = getoptions(params,'verInd',1:n);
scattColor = getoptions(params,'scattColor','c');
scattSize = getoptions(params,'scattSize',500);

% check size of verInd compared to n
if n ~= size( verInd, 2 )
    n = size( verInd, 2 );
end
%============================================

%--------------------------------------------
% Take out NaNs
%--------------------------------------------
nanInd = isnan(verInd);
%============================================

%--------------------------------------------
% Extract verteces and faces information
%--------------------------------------------
axis equal; axis off;

patch('vertices',V','faces',F','FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',1);
camlight
hold on
scatter3(V(1,verInd(~nanInd)), V(2,verInd(~nanInd)), V(3,verInd(~nanInd)),scattSize,scattColor,'filled');

axis equal, axis off
%============================================

end




