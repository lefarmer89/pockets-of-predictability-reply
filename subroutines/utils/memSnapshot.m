function rssGB = memSnapshot(label)
% MEMSNAPSHOT  Print and return the calling MATLAB process's resident set
% size in GB. Used to instrument peak-memory checkpoints inside the heavy
% hyperparameter-sweep entry points.
%
%   rssGB = memSnapshot('post-load')
%
% Implementation: shells out to `ps -o rss=` (macOS / Linux). On systems
% where the call fails the function returns NaN and prints a warning-free
% line so the calling code keeps running.

pid = matlabProcessID;
[status, out] = system(sprintf('ps -o rss= -p %d', pid));
if status == 0
    rssKB = str2double(strtrim(out));
    rssGB = rssKB / 1024 / 1024;
else
    rssGB = NaN;
end
fprintf('  [mem] %-36s RSS = %5.1f GB\n', label, rssGB);
end
