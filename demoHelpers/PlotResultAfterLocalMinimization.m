function PlotResultAfterLocalMinimization(V1,F1,V2,F2,idx1,idx2,title1,title2)
figure
colorVerts = V1(:,idx1)';
scattColor = bsxfun(@rdivide,bsxfun(@minus,colorVerts,min(colorVerts,[],1)), (max(colorVerts,[],1)-min(colorVerts,[],1)));
figure('Position', [100, 100, 800, 800]);
subplot(1,2,1)
title(title1)
params.scattColor = scattColor;
params.verInd  = idx1;
plotMeshAndPoints( V1, F1, params )
subplot(1,2,2)
title(title2)
params.scattColor = scattColor;
params.verInd  = idx2;
plotMeshAndPoints( V2, F2, params )

end