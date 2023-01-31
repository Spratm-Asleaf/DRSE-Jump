%% Generate Main Contents of Latex Tables

fid = fopen('LatexTable.txt','w');

fprintf(fid, 'Scenario: %d\n', scenario_index);

if scenario_index <= 3
    fprintf(fid, 'KF*    & %.2f & %.2e & IMM-T* & %.2f & %.2e \\\\\n', RMSE_KF, time_KF, RMSE_IMM_true, time_IMM_true);
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'IMM-N  & %.2f & %.2e & IMM-B  & %.2f & %.2e \\\\\n', RMSE_IMM_nominal, time_IMM_nominal, RMSE_IMM_Bayesian, time_IMM_Bayesian);
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'IMM-M  & %.2f & %.2e & IMM-C  & %.2f & %.2e \\\\\n', RMSE_IMM_MaximumLikelihood, time_IMM_MaximumLikelihood, RMSE_IMM_compensation, time_IMM_compensation);
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'IMM-RS & %.2f & %.2e & DRIMM  & %.2f & %.2e \\\\\n', RMSE_IMM_RiskSensitive, time_IMM_RiskSensitive, RMSE_DRIMM, time_DRIMM);
elseif scenario_index == 4
    fprintf(fid, 'KF*    & %.2f & %.2e & IMM-T* & %.2f & %.2e \\\\\n', RMSE_KF, time_KF, RMSE_IMM_true, time_IMM_true);
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'IMM-N  & %.2f & %.2e & IMM-C  & %.2f & %.2e \\\\\n', RMSE_IMM_nominal, time_IMM_nominal, RMSE_IMM_compensation, time_IMM_compensation);
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'IMM-RS & %.2f & %.2e & DRIMM  & %.2f & %.2e \\\\\n', RMSE_IMM_RiskSensitive, time_IMM_RiskSensitive, RMSE_DRIMM, time_DRIMM);
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'IMM-R  & %.2f & %.2e &        &      &      \\\\\n', RMSE_IMM_RobustKF, time_IMM_RobustKF);
elseif scenario_index >= 5
    fprintf(fid, 'KF*    & %.2f & %.2e & IMM-T* & %.2f & %.2e \\\\\n', RMSE_KF, time_KF, RMSE_IMM_true, time_IMM_true);
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'IMM-N  & %.2f & %.2e & IMM-B  & %.2f & %.2e \\\\\n', RMSE_IMM_nominal, time_IMM_nominal, RMSE_IMM_Bayesian, time_IMM_Bayesian);
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'IMM-M  & %.2f & %.2e & IMM-C  & %.2f & %.2e \\\\\n', RMSE_IMM_MaximumLikelihood, time_IMM_MaximumLikelihood, RMSE_IMM_compensation, time_IMM_compensation);
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'IMM-RS & %.2f & %.2e & DRIMM  & %.2f & %.2e \\\\\n', RMSE_IMM_RiskSensitive, time_IMM_RiskSensitive, RMSE_DRIMM, time_DRIMM);
    fprintf(fid, '\\hline\n');
    fprintf(fid, 'IMM-R  & %.2f & %.2e &        &      &      \\\\\n', RMSE_IMM_RobustKF, time_IMM_RobustKF);
end

fclose(fid);