#!/bin/bash 
foldername=/home/steven/thesis/thesis_config/config/tester
ini_file=${foldername}/tester_oms8.ini
multinest_file=/home/steven/thesis/thesis_config/outputs/run_cluster10_SN_BAO_QSO_Marie/output_multinest.txt    # This can in principle be any file with the same columns as the multinest file (doesn't need the weights though). So it can be a list of parameter values with the same ordering that would be in a output multinest file produced by cosmosis usign the same inifile as ${ini_file}. Then code will look for the columns called post and like to find the row with the best fit values and starts from there. Alternatively if this is not given will look in the values.ini file that is named in the ${ini_file} for the start values. 
output_file_name=/home/steven/thesis/thesis_config/outputs/run_cluster10_SN_BAO_QSO_Marie/output_nelder_mead
bestfit_file_name=/home/steven/thesis/thesis_config/outputs/run_cluster10_SN_BAO_QSO_Marie/bestfit_nelder_mead
bestfit_prior_name=/home/steven/thesis/thesis_config/outputs/run_cluster10_SN_BAO_QSO_Marie/bestfit_priors # If  --best_fit_value set but this is not set it will still produce somethig with the same name but a "prior_" added to the start

# If -s is removed it will only save the best fit result and not the steps.

python maxlike_cosmosis.py -i ${ini_file} -m ${multinest_file} -o ${output_file_name} --max_post -s --best_fit_value ${bestfit_file_name} --best_fit_priors ${bestfit_prior_name} --maxiter 3000  

