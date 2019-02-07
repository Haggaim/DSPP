%-------------------------------------------------------------------------
% Demo for "DS++: A flexible, scalable and provably tight relaxation for matching problems"
% using two shapes from: http://www.wisdom.weizmann.ac.il/~ylipman/CPsurfcomp/
%The demo matches n points on one tooth to k<=n points on the second tooth
%---------------------------------------------------------------------------------
% Code written by Haggai Maron and Nadav Dym. Please cite
% "DS++: A flexible, scalable and provably tight relaxation for matching problems"
% questions regarding the code can be sent to nadavdym@gmail.com
%---------------------------------------------------------------------------------
addpath(genpath(pwd));
n=110;
k=100;
useGeodesicDistance=0;
% load shapes
[V1,F1] = read_off('a10_sas.off'); V1=V1';F1=F1';
[V2,F2] = read_off('a13_sas.off'); V2=V2';F2=F2';

if useGeodesicDistance
    adj1 = triangulation2adjacency(F1,V1');
    fullD1 = graphallshortestpaths(adj1,'directed',false);    
    adj2 = triangulation2adjacency(F2,V2');
    fullD2 = graphallshortestpaths(adj2,'directed',false);
else
    % calc euclidean distance matrices
    fullD1 = squareform(pdist(V1));
    fullD2 = squareform(pdist(V2));
end

% sample points
idx1 = chooseFarthestPoints(fullD1,n);
idx2 = chooseFarthestPoints(fullD2,k);


% call DSPP
D1=fullD1(idx1,idx1);
D2=fullD2(idx2,idx2);
%problemType='GW';
problemType='L1';

[X_final,X_opt,orig_obj,opt_obj] = solveDSpp(D1,D2,problemType);


% visualize correspondence
final_idx2 = idx2;
final_idx1 = X_final * idx1';
PlotResultAfterLocalMinimization(V2',F2',V1',F1',final_idx2,final_idx1,'source','target');