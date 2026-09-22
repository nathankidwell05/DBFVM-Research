function [solverFolder,solverFile] = project_paths
%PROJECT_PATHS finds the solver folder from a report scripts folder.
%   The OneDrive research folder and the GitHub repository now share one
%   layout, so a single rule locates the solver in both:
%       <root>/matlab/d1v5/<solver file>
%       <root>/reports/<report name>/scripts/   (this file)
%   Only the solver file name differs between them: the OneDrive copy is
%   d1v5_sodshock_muscl_MC_RK3.m and the repository copy is
%   d1v5_sodshock_gminmod_rk3.m. Both names are tried below.
%   Results are read from and written to <solverFolder>/results in both.

    here = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(fileparts(fileparts(here))); % up out of reports/<name>/scripts
    solverFolder = fullfile(projectRoot,'matlab','d1v5');

    candidateNames = {'d1v5_sodshock_muscl_MC_RK3.m', ... % OneDrive research folder
                      'd1v5_sodshock_gminmod_rk3.m'};     % GitHub repository
    for nameIndex = 1:numel(candidateNames)
        solverFile = fullfile(solverFolder,candidateNames{nameIndex});
        if isfile(solverFile); return; end
    end
    error('Could not locate the KT-D1V5 solver in: %s',solverFolder);
end
