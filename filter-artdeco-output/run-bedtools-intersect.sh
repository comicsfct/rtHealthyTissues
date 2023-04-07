# bash script

tissue_folder=$1 # standard artdeco output folder
python_script=$2 # costum python script

# The first input is the standard ARTDeco output folder, containing the dogs folder to filter
# The second input is the custom python script used to create a new dog-filtered folder containing the new filtered raw/fpkm files

# folders of interest are found directly from the input folder
dogs_folder=${tissue_folder}/dogs/
all_genes_bed=${tissue_folder}/preprocess_files/genes_condensed.bed

# create a filtered list of dogs for each .bed file

for file in $(ls $dogs_folder/*.bed); do

	echo "creating ${file}.filtered.intersect";

	bedtools intersect -a $file -b $all_genes_bed -v -S > ${file}.filtered.intersect
done

# -v Only report those entries in A that have no overlap in B. Restricted by -f and -r.
# -S Require different strandedness. That is, only report hits in B that overlap A on the _opposite_ strand.

# run custom python script

/mnt/c/Users/paulo/anaconda3//python.exe $python_script $tissue_folder

# temporary files with the suffix .filtered.intersect are created inside the original dog folder while the script is running
