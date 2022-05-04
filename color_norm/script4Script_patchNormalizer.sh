#!/bin/bash
source_dir=/proj/STOR/yaoli/data/melanoma/patches_all/patches/
target_dir=/pine/scr/y/a/yaoli/data/melanoma/patches_all_norm
jobs=./jobs_folder
mkdir $jobs
for patch_file in /proj/STOR/yaoli/data/melanoma/patches_all/patches/*.h5     # list directories in the form "/tmp/dirname/"
do
    # remove the trailing "/"
    echo $patch_file
    
    filename=$(basename $patch_file .h5)
    echo "processing $filename file"
    ext='.job'
    file=$filename$ext
    echo $file
    
    ext_h5='.h5'
    resNorm=$target_dir/patches/$filename$ext_h5

    echo '#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem 16g
#SBATCH -t 2:00:00

python normalizeStaining_patches.py -i '$patch_file' -o '$resNorm''> $file
        chmod +x $file
        mv $file jobs_folder/

done
