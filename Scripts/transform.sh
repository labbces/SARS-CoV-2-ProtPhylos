#! /bin/bash

sed -i 's/\//***/g' Oceania2.fasta
sed -i 's/|/^^^/g' Oceania2.fasta

filename="Oceania2.fasta"

while read -r line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile="{$line#>}.fasta"
        echo "$line" > "$outfile"
    else
        echo "$line" >> "$outfile"
    fi
done <"$filename"
