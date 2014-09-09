function [diff_vec,data_diff, data_A, data_B] = CompareLamemDebugOutput(folder1,folder2)
%
%   Compare LaMEM debug output that was created  with
%
%   ierr = DebugSave2Bin_GlobalVec();CHKERRQ(ierr);
%
%   or
%
%   ierr = DebugSave2Bin_GlobalMat();CHKERRQ(ierr);
%
%   call it like:
%   [diff_vec,data_diff, data_A, data_B] = CompareLamemDebugOutput('initialrun','restartrun')
%
%   Just give the names of the two directories that contain the data you want 
%   to compare with each other.
%   Contributed by Tobias Baumann(Mainz university, Sep 2012)
% -------------------------------------------------------------------------


current_path= pwd;
dir1_path   = [current_path '/' folder1];
dir2_path   = [current_path '/' folder2];


% --- first folder ---
data_A = PetscBin2struct('Debug_*.bin',dir1_path,current_path);


% --- second folder ---
data_B = PetscBin2struct('Debug_*.bin',dir2_path,current_path);


% --- get differences ---
[data_diff,diff_vec] = ComparePetscBinStructs(data_A,data_B);

end
% --- end of function -----------------------------------------------------


% -------------------------------------------------------------------------
function [data_struct] = PetscBin2struct(pattern,dir_path,current_path) 

cd(dir_path);

files = dir(pattern);
for k=1:length(files)
    varname = genvarname([files(k).name(7:end-4)]);
    eval([varname ' = PetscBinaryRead(files(k).name);']);
    data_struct.(genvarname(varname)) = eval(varname);
    clear(genvarname(varname));
end

cd(current_path);
end

% -------------------------------------------------------------------------
function [data_diff,diff_vec] = ComparePetscBinStructs(data_A,data_B)

fields_B = fieldnames(data_B);
fields_A = fieldnames(data_A);

n =1;
for k = 1:length(fields_B)
    test = find(strcmpi(fields_B(k),fields_A));   
    if ~isempty(test)
        samefields(n) = test;
        n = n+1;
    end
end

for k = 1:length(samefields)
    data_diff.(genvarname(fields_A{samefields(k)})) =  data_A.(genvarname(fields_A{samefields(k)})) - data_B.(genvarname(fields_B{k})) ;
    diff_vec(k) = sum(sum(data_diff.(genvarname(fields_A{samefields(k)}))));
end

fields_diff = fieldnames(data_diff);

if sum(diff_vec)~= 0
    for k = find(diff_vec ~= 0)
        disp(['There are differences in field ' fields_diff{k}]);
    end
else
    disp('There are no differences');
end
end
% -------------------------------------------------------------------------