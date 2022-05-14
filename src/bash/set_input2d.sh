#!/bin/bash

WD=${1:?Provide a working directory}
Species=${2:?Provide a species name}

# Separate main parameter file into three files
cut -d , -f 1 "$WD"/data/parameters/parameters_runs.txt > "$WD"/data/input2d-files/number.txt
cut -d , -f 2 "$WD"/data/parameters/parameters_runs.txt > "$WD"/data/input2d-files/speed.txt

# Count number of lines in files
numlines=$(grep -c "^" "$WD"/data/parameters/parameters_runs.txt)
mkdir "$WD"/data/input2d-files
mkdir "$WD"/data/input2d-files/"$Species"

cd "$WD"/data/input2d-files

# initialize variables
number=0
speed=0
# For loop that will write files
for i in `seq 1 $numlines`;
do
# Sets Wo based on i
number=$(awk -v var="$i" 'NR==var' number.txt)
numspec=${Species}_${number}
# Writes file to replace line with vertex file number
awk -v var="$numspec" 'NR==118 {$0="   structure_names = \""'"var"'"\"    //"} 1' template-input2d > input2d_w1_${i}
awk -v var="$numspec" 'NR==119 {$0="   "'"var"'" {   // "} 1' input2d_w1_${i} > input2d_w2_${i}
awk -v var="$numspec" 'NR==221 {$0="   structure_names  = \""'"var"'"\"   //  "} 1' input2d_w2_${i} > input2d_w3_${i}
# Sets Freq based on i
speed=$(awk -v var="$i" 'NR==var' speed.txt)
# Writes file to replace line with speed
awk -v var="$speed" 'NR==5 {$0="FS =  "'"var"'"				  // Flow speed of tank (m/s)"} 1' input2d_w3_${i} > input2d_w4_${i}
awk -v var="$speed" 'NR==300 {$0="FS="'"var"'"    //"} 1' input2d_w4_${i} > input2d_w5_${i}
# Edits input2d to create different IBlog and visit files and folders
results_dir="../results/ibamr"
log_dump=${results_dir}/${Species}/"log-files/IB2d.log"${i}
visit_dir=${results_dir}/${Species}/"viz_IB2d"${i}
hier_dir=${results_dir}/${Species}/"hier_data_IB2d"${i}
awk -v var="$log_dump" 'NR==255 {$0="   log_file_name = \""'"var"'"\"         //"} 1' input2d_w5_${i} > input2d_w6_${i}
awk -v var="$visit_dir" 'NR==261 {$0="   viz_dump_dirname = \""'"var"'"\"                //"} 1' input2d_w6_${i} > input2d_w7_${i}
awk -v var="$hier_dir" 'NR==270 {$0="   data_dump_dirname = \""'"var"'"\"          //"} 1' input2d_w7_${i} > input2d_w8_${i}

# Cleans up folder
mv input2d_w8_${i} "$WD"/data/input2d-files/"$Species"/input2d${i}
rm input2d_w*_${i}


echo $i

done

#rm number.txt speed.txt

