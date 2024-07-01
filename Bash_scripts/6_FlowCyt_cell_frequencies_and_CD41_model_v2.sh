#!/bin/bash

eval "$(conda shell.bash hook)"



Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2

Master_path_analysis=$(echo "$MASTER_ROUTE""$analysis""/")

#rm -rf $Master_path_analysis
#mkdir -p $Master_path_analysis



Log_files=$(echo "$Master_path_analysis""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files



#### data_wrangling ####

module load R/4.1.0

type=$(echo "data_wrangling""_""$analysis")
outfile_data_wrangling=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_data_wrangling
echo -n "" > $outfile_data_wrangling
name_data_wrangling=$(echo "$type""_job")


Rscript_data_wrangling=$(echo "$Rscripts_path""22_FlowCyt_data_wrangling_cell_frquencies.R")


FlowCyt_results=$(echo "/group/soranzo/manuel.tardaguila/FlowCyt_parameters/""FlowCyt_results_for_model_added_all_quadrants_v2.csv")
selection_samples=$(echo "WT_A,WT_B,WT_C,clone_13,clone_27,clone_29,del_233,del_235,del_287")
selection_treatment=$(echo "5nM_PMA")
selection_time_point=$(echo "Basal,16_hrs,24_hrs,48_hrs,72_hrs")



myjobid_data_wrangling=$(sbatch --job-name=$name_data_wrangling --output=$outfile_data_wrangling --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=512M --parsable --wrap="Rscript $Rscript_data_wrangling --FlowCyt_results $FlowCyt_results --Master_path_analysis $Master_path_analysis --selection_samples $selection_samples --selection_treatment $selection_treatment --selection_time_point $selection_time_point --type $type --out $Master_path_analysis")
myjobid_seff_data_wrangling=$(sbatch --dependency=afterany:$myjobid_data_wrangling --open-mode=append --output=$outfile_data_wrangling --job-name=$(echo "seff""_""$name_data_wrangling") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_data_wrangling >> $outfile_data_wrangling")

#### clr transformation of cell frequencies ####

module load R/4.1.0

type=$(echo "ilr_transformation_and_princomp""_""$analysis")
outfile_ilr_transformation_and_princomp=$(echo "$Log_files""outfile_2_""$type"".log")
touch $outfile_ilr_transformation_and_princomp
echo -n "" > $outfile_ilr_transformation_and_princomp
name_ilr_transformation_and_princomp=$(echo "$type""_job")


Rscript_ilr_transformation_and_princomp=$(echo "$Rscripts_path""23_ilr_transformation_cell_frequencies_and_princomp_v3.R")


FlowCyt_results_subset=$(echo "$Master_path_analysis""FlowCyt_global_adapted.rds")

# --dependency=afterany:$myjobid_data_wrangling

myjobid_ilr_transformation_and_princomp=$(sbatch --dependency=afterany:$myjobid_data_wrangling --job-name=$name_ilr_transformation_and_princomp --output=$outfile_ilr_transformation_and_princomp --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=512M --parsable --wrap="Rscript $Rscript_ilr_transformation_and_princomp --FlowCyt_results_subset $FlowCyt_results_subset --Master_path_analysis $Master_path_analysis --type $type --out $Master_path_analysis")
myjobid_seff_ilr_transformation_and_princomp=$(sbatch --dependency=afterany:$myjobid_ilr_transformation_and_princomp --open-mode=append --output=$outfile_ilr_transformation_and_princomp --job-name=$(echo "seff""_""$name_ilr_transformation_and_princomp") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_ilr_transformation_and_princomp >> $outfile_ilr_transformation_and_princomp")



#### models ####

module load R/4.1.0

type=$(echo "models""_""$analysis")
outfile_models=$(echo "$Log_files""outfile_3_""$type"".log")
touch $outfile_models
echo -n "" > $outfile_models
name_models=$(echo "$type""_job")


Rscript_models=$(echo "$Rscripts_path""25_FlowCyt_Cell_frequencies_models_v3.R")

FlowCyt_results_subset=$(echo "$Master_path_analysis""FlowCyt_global_adapted_plus_ilr_and_princomp.rds")

# --dependency=afterany:$myjobid_ilr_transformation_and_princomp

myjobid_models=$(sbatch --dependency=afterany:$myjobid_ilr_transformation_and_princomp --job-name=$name_models --output=$outfile_models --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_models --FlowCyt_results_subset $FlowCyt_results_subset --Master_path_analysis $Master_path_analysis --type $type --out $Master_path_analysis")
myjobid_seff_models=$(sbatch --dependency=afterany:$myjobid_models --open-mode=append --output=$outfile_models --job-name=$(echo "seff""_""$name_models") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_models >> $outfile_models")


#### graphs ####

module load R/4.1.0

type=$(echo "graphs""_""$analysis")
outfile_graphs=$(echo "$Log_files""outfile_4_""$type"".log")
touch $outfile_graphs
echo -n "" > $outfile_graphs
name_graphs=$(echo "$type""_job")


Rscript_graphs=$(echo "$Rscripts_path""24_FlowCyt_Cell_frequencies_graphs_v3.R")

FlowCyt_results_subset=$(echo "$Master_path_analysis""FlowCyt_global_adapted_plus_ilr_and_princomp.rds")
model_results=$(echo "$Master_path_analysis""models/compositional_clr/compositional_ilr_model_result.tsv")


# --dependency=afterany:$myjobid_models

myjobid_graphs=$(sbatch --dependency=afterany:$myjobid_models --job-name=$name_graphs --output=$outfile_graphs --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_graphs --FlowCyt_results_subset $FlowCyt_results_subset --Master_path_analysis $Master_path_analysis --model_results $model_results  --type $type --out $Master_path_analysis")
myjobid_seff_graphs=$(sbatch --dependency=afterany:$myjobid_graphs --open-mode=append --output=$outfile_graphs --job-name=$(echo "seff""_""$name_graphs") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_graphs >> $outfile_graphs")


