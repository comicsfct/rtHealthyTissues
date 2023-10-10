
# select a folder containing bed files with the coordinates to use
# I wrote this script to apply to the bed_files_flanks folder

sample=$1/ # folder containing bed files

# this is the folder containing all the epigenome models - sample name will be used to search for the right model reference
models_samples=models_samples/
output_folder=$2 # OverlapEnrichment/


# define variables
name=$(basename $sample | rev |cut -d'_' -f2- | rev)
model=$(ls $models_samples/ | grep ${name})
output=${output_folder}/$(basename $sample)_OverlapEnrichment

echo "PROCESSING $name"
echo "MODEL $model"
echo "OUTPUT $output"

# run ChromoHMM on given sample
java -Xmx10G -jar ChromoHMM/ChromHMM.jar OverlapEnrichment -b 1 $models_samples/$model $sample $output
