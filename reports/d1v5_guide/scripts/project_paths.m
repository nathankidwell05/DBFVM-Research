function [solverFolder,solverFile] = project_paths
%PROJECT_PATHS finds the solver folder for either folder layout.
%   OneDrive research folder:  <root>/report_d1v5_math/scripts  and  <root>/d1v5_sodshock_muscl_MC_RK3.m
%   GitHub repository:         <repo>/reports/d1v5_mathematics/scripts  and  <repo>/matlab/d1v5/d1v5_sodshock_gminmod_rk3.m
% Results are read from and written to <solverFolder>/results in both layouts.

    here = fileparts(mfilename('fullpath'));
    oneDriveRoot = fileparts(fileparts(here));
    repoSolverFolder = fullfile(fileparts(fileparts(fileparts(here))),'matlab','d1v5');
    if isfile(fullfile(oneDriveRoot,'d1v5_sodshock_muscl_MC_RK3.m'))
        solverFolder = oneDriveRoot;
        solverFile = fullfile(oneDriveRoot,'d1v5_sodshock_muscl_MC_RK3.m');
    elseif isfile(fullfile(repoSolverFolder,'d1v5_sodshock_gminmod_rk3.m'))
        solverFolder = repoSolverFolder;
        solverFile = fullfile(repoSolverFolder,'d1v5_sodshock_gminmod_rk3.m');
    else
        error('Could not locate the KT-D1V5 solver from %s',here);
    end
end
