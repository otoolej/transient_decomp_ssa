%-------------------------------------------------------------------------------
% add_path_here: add current path will all subdirectories to Matlab path
%
% Syntax: add_path_here(path_to_add)
%
% Inputs: 
%     path_to_add - name of path to add; if empty, assumes current directory
%
%
% Example:
%     add_path_here();
%

% John M. O' Toole, University College Cork
% Started: 16-02-2021
%
% last update: Time-stamp: <2021-02-19 11:58:31 (otoolej)>
%-------------------------------------------------------------------------------
function add_path_here(path_to_add)
if(nargin<1 || isempty(path_toadd)), path_toadd=pwd; end


% add all sub-dirs. apart from the following list:
addpath(genpath_selective(path_toadd, {'.git', 'old', 'OLD', 'data', 'codegen', 'pics'}));




function list_dirs=genpath_selective(dir_add, dirs_remove)
%-------------------------------------------------------------------------------
% Add subdirectories, excluding all those in 'dirs_remove'
%-------------------------------------------------------------------------------
if(nargin<2 || isempty(dirs_remove)), dirs_remove={'.git', 'old'}; end

if(~isempty(dirs_remove) && ~iscell(dirs_remove))
    dirs_remove={dirs_remove};
end

p=genpath(dir_add);

% extract paths from this:
subdirs=strsplit(p, ':');

if(length(subdirs{end})<2), subdirs(end)=[]; end

list_dirs='';
for n=1:length(subdirs)
    
    iexclude_dir=[];
    for p=1:length(dirs_remove)
        iexclude_dir{p}=strfind(subdirs{n}, [filesep dirs_remove{p}]);
    end

    if(all(cellfun(@isempty, iexclude_dir)))
        list_dirs=[list_dirs subdirs{n} ':'];
    end
end

DB = 1;
if(DB)
    fprintf('Adding paths:\n');
    str_ldirs=strsplit(list_dirs, ':');
    for n=1:length(str_ldirs)-1
        disp(['(' num2str(n) '): ' str_ldirs{n}]);
    end
end
