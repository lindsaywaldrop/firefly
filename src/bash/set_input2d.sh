#!/bin/bash

WD=${1:?Provide a working directory}
Species=${2:?Provide a species name}

# Separate main parameter file into three files
cut -d , -f 1 "$WD"/data/parameters/parameters_runs.txt > "$WD"/data/input2d-files/number.txt
cut -d , -f 2 "$WD"/data/parameters/parameters_runs.txt > "$WD"/data/input2d-files/speed.txt

# Count number of lines in files
numlines=$(grep -c "^" "$WD"/data/parameters/parameters_runs.txt)
cd "$WD"/data/input2d-files

# initialize variables
number=0
speed=0
# For loop that will write files
for i in `seq 1 $numlines`;
do
# Sets Wo based on i
number=$(awk -v var="$i" 'NR==var' number.txt)
# Writes file to replace line with vertex file number
awk -v var="$number" 'NR==118 {$0="   structure_names = \"Lucidota_sp_"'"var"'"\"    //"} 1' template-input2d > input2d_w1_${i}
awk -v var="$number" 'NR==119 {$0="   Lucidota_sp_"'"var"'" {   // "} 1' input2d_w1_${i} > input2d_w2_${i}
awk -v var="$number" 'NR==221 {$0="   structure_names  = \"Lucidota_sp_"'"var"'"\"   //  "} 1' input2d_w2_${i} > input2d_w3_${i}
# Sets Freq based on i
speed=$(awk -v var="$i" 'NR==var' speed.txt)
# Writes file to replace line with speed
awk -v var="$speed" 'NR==5 {$0="FS =  "'"var"'"				  // Flow speed of tank (m/s)"} 1' input2d_w3_${i} > input2d_w4_${i}
awk -v var="$speed" 'NR==300 {$0="FS="'"var"'"    //"} 1' input2d_w4_${i} > input2d_w5_${i}
# Edits input2d to create different IBlog and visit files and folders
awk -v var="$i" 'NR==255 {$0="   log_file_name = \"../results/ibamr/log-files/IB2d.log"'"var"'"\"                //"} 1' input2d_w5_${i} > input2d_w6_${i}
awk -v var="$i" 'NR==261 {$0="   viz_dump_dirname = \"../results/ibamr/runs/viz_IB2d"'"var"'"\"                //"} 1' input2d_w6_${i} > input2d_w7_${i}
awk -v var="$i" 'NR==270 {$0="   data_dump_dirname = \"../results/ibamr/runs/hier_data_IB2d"'"var"'"\"          //"} 1' input2d_w7_${i} > input2d_w8_${i}

# Cleans up folder
mv input2d_w8_${i} "$WD"/data/input2d-files/input2d${i}
rm input2d_w*_${i}


echo $i

done

#rm number.txt speed.txt
