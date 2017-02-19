#!/bin/bash
#This script will take in the user's input for the Peptide chain, create a Polypeptide using a software created by the Pappu Lab known as CAMPARI and determine the number of simulations the user would like to perform
declare -a aarray
echo "Type the number of Amino Acids that your peptide chain will contain:"
read numAA
echo "Please Enter the Amino Acid sequence of the peptide you wish to simulate:"
read -a aarray

echo "File name?"
read Filename

declare -i num_pos=0
declare -i num_neg=0
declare -i sum=0
declare -i antisum=0

for (( i=0; i<$numAA; i++))
do
	if [ "${aarray[$i]}" == "ARG" ]|| [ "${aarray[$i]}" == "HIS" ] || [ "${aarray[$i]}" == "LYS" ]
	then
		num_pos=$(($num_pos + 1))
		echo ${aarray[$i]}
		echo $num_pos
	fi
	if [ "${aarray[$i]}" == "ASP" ]|| [ "${aarray[$i]}" == "GLU" ]
	then
		num_neg=$(($num_neg + 1))
		echo ${aarray[$i]}
		echo $num_neg
	fi
done

let sum=$((num_pos-num_neg))
echo $sum
		
mkdir /work/kbharat/Traj-Simulations/$Filename
cp /work/kbharat/template.key /work/kbharat/Traj-Simulations/$Filename/template.key
echo "ACE" >> /work/kbharat/Traj-Simulations/$Filename/seq.in

for (( k=0; k<=($numAA-1); k++))
do
	echo ${aarray["$k"]} >> /work/kbharat/Traj-Simulations/$Filename/seq.in
done

echo "NME" >> /work/kbharat/Traj-Simulations/$Filename/seq.in

if [ $sum < 0 ]
then
	let antisum=$((sum*-1))
	for  (( q=1; q<=($antisum); q++))
	do
		echo "K+" >> /work/kbharat/Traj-Simulations/$Filename/seq.in
	done
fi

if [ $sum > 0 ]
then
        for (( q=1; q<=($sum); q++))
        do
                echo "CL-" >> /work/kbharat/Traj-Simulations/$Filename/seq.in
        done
fi		



for (( j=1; j<=(267); j++)) # The number of Water Molecules (T3P) have been optimized based on calculations performed on previous simulations
do
	echo "t3p" >> /work/kbharat/Traj-Simulations/$Filename/seq.in
done 

echo "END" >> /work/kbharat/Traj-Simulations/$Filename/seq.in

echo "How many simulations would you like to perform?"
read Simnum

for (( l=1; l<=$Simnum; l++))
do	
	cd /work/kbharat/Traj-Simulations/$Filename/
	mkdir /work/kbharat/Traj-Simulations/$Filename/$Filename"$l"
	cp /work/kbharat/Traj-Simulations/$Filename/seq.in /work/kbharat/Traj-Simulations/$Filename/$Filename"$l"/seq.in
	cp /work/kbharat/Traj-Simulations/$Filename/template.key /work/kbharat/Traj-Simulations/$Filename/$Filename"$l"/template.key
	cd /work/kbharat/Traj-Simulations/$Filename/$Filename"$l"
	/packages/campari_new/bin/x86_64/campari -k template.key > log$Filename"$l"
        mv /work/kbharat/Traj-Simulations/$Filename/$Filename"$l"/POLYQ_START.pdb /work/kbharat/Traj-Simulations/$Filename/$Filename"$l"/"$Filename""$l".pdb
        rm /work/kbharat/Traj-Simulations/$Filename/$Filename"$l"/POLYQ_END.int
        rm /work/kbharat/Traj-Simulations/$Filename/$Filename"$l"/POLYQ_START.int
        rm /work/kbharat/Traj-Simulations/$Filename/$Filename"$l"/POLYQ_END.pdb
        cd ..
done

mkdir /work/kbharat/Traj-Simulations/$Filename/"$Filename"PDBs
for (( m=1; m<=$Simnum; m++))
do
	cp /work/kbharat/Traj-Simulations/$Filename/$Filename"$m"/"$Filename""$m".pdb /work/kbharat/Traj-Simulations/$Filename/"$Filename"PDBs/"$Filename""$m".pdb
	sed -i 's/T3P/HOH/g' /work/kbharat/Traj-Simulations/$Filename/"$Filename"PDBs/"$Filename""$m".pdb
        sed -i 's/1HW/HW1/g' /work/kbharat/Traj-Simulations/$Filename/"$Filename"PDBs/"$Filename""$m".pdb
        sed -i 's/2HW/HW2/g' /work/kbharat/Traj-Simulations/$Filename/"$Filename"PDBs/"$Filename""$m".pdb
done

scp -r /work/kbharat/Traj-Simulations/$Filename/"$Filename"PDBs/ krishbharat96@login01.chpc.wustl.edu:/scratch/krishbharat96/Pappu-Simulations/TrajSims/$Filename

	
