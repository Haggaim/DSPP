function idx = chooseFarthestPoints(dist,n)

[~, idx]  = max(dist(1,:));
for ii = 2:n
  [~, idx(ii)]  = max(min(dist(idx,:),[],1));
end

end
