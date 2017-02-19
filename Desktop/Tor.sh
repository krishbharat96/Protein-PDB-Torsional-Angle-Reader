#!/bin/bash
declare -a carray
declare -a narray
declare -a calphaarray
declare -a cbarray
declare -a nbarray

declare -i ccount=0
declare -i ncount=0

declare -i cbvalue=0
declare -i cavalue=0
declare -i nbvalue=0
declare -i navalue=0

declare -i Ccounter=0
declare -i Ncounter=0

while read line 
do
	read -ra atom <<< "$line"

	if [ "${atom[2]}" == "C" ] 
	then
		carray[${#carray[@]}]=${atom[1]}
	
		cbvalue=$((${#carray[@]} - 2))
		cavalue=$((${#carray[@]} - 1))
		
		#echo "Number of Elements in C0Array:"
		#echo ${#carray[@]}
		#echo $cbvalue
		#echo $cavalue
		
		if (($Ccounter > 0)) 
		then
			let cbarray[$cbvalue]=${carray[$cavalue]}	
		fi

		#echo "Number of Elements in C1Array:"
                #echo ${#cbarray[@]}
		#echo 'End of one Loop!'

		Ccounter=$(($Ccounter + 1))
	fi

	if [ "${atom[2]}" == "N" ]
        then
		narray[${#narray[@]}]=${atom[1]}
                
                nbvalue=$((${#narray[@]} - 2))
                navalue=$((${#narray[@]} - 1))

                if (($Ncounter > 0))
                then
                        let nbarray[$nbvalue]=${narray[$navalue]}
		fi
		Ncounter=$(($Ncounter + 1))
	fi

	if [ "${atom[2]}" == "CA" ]
        then
		calphaarray[${#calphaarray[@]}]=${atom[1]}
	fi

done < 5ALA.pdb

declare -i count=0

count=$(($Ccounter - 1))

Linenum=$(wc -l < 5ALA.pdb)

for (( i=0; i<$count; i++))
do
	j=$(($i + 1))

	echo "Phi-Psi Angle Location $j :"
	echo "C0: ${carray[$i]}"
	echo "N1: ${narray[$i]}"
	echo "CA: ${calphaarray[$i]}"
	echo "C1: ${cbarray[$i]}"
	echo "N2: ${nbarray[$i]}"
	
	cp parse parse$j
	cp parse.c parse$j.c
	
	sed -i '6s/framenum/5000/i' parse$j.c
	sed -i '7s/Atomsnum/'$Linenum'/i' parse$j.c
	sed -i '8s/C0num/'${carray[$i]}'/i' parse$j.c
	sed -i '9s/N1num/'${narray[$i]}'/i' parse$j.c
	sed -i '10s/CAnum/'${calphaarray[$i]}'/i' parse$j.c
	sed -i '11s/C1num/'${cbarray[$i]}'/i' parse$j.c
	sed -i '12s/N2num/'${nbarray[$i]}'/i' parse$j.c  	 

	gcc parse$j.c -o parse$j -lm
	./parse$j 5ALA_traj.dcd > 5ALA_ParsedTraj$j.txt

done

