function mergecache(old,new)

% Combine data from two MAT-files
load(old);
load(new);

%% cache file after merging a and b
h=datestr(clock,0);
cacheFile=['cacheMerge_',h(1:11),'_',h(13:14),'_',h(16:17),'_',h(19:20)];
fullcacheFile=[pwd '/' cacheFile '.mat'];
if exist(fullcacheFile,'file')==0
    save(cacheFile);
else
    warning('Overriding cache file');
end
