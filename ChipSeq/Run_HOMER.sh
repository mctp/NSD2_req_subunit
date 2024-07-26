#Activate HOMER enviroment and test path 

#Make additional bed files for top scoring peaks and/or center peaks 200bp width
for i in *.bed
do echo $i
sort -k5nr $i | head -n 1001 > Top1k_$i
awk '{print $1"\t"($3-$2)/2+$2-100"\t"($3-$2)/2+$2+100"\t"}' Top1k_$i > Top1k_200bp$i
done
 
#Launch HOMER runs. 
#Custom motif database: nohup perl findMotifsGenome.pl  $bed1 hg38 ./$out1/ -size given -preparse -mknown ../NewHybridsAdded.motifs >Log$out1.log 2>&1 &
for i in *.bed 
do echo $i
out1="${i%%.*}"
bed1=$i
nohup perl findMotifsGenome.pl  $bed1 hg38 ./$out1/ -size given -preparse >Log$out1.log 2>&1 &
done

#Move Output Around 
mkdir Known
mkdir Denovo
mkdir Table
find  -name known*html | xargs -I {} --verbose sh -c 'cp "{}" "$(basename "$(dirname "{}")")_$(basename "{}")"'
find  -name homer*html | xargs -I {} --verbose sh -c 'cp "{}" "$(basename "$(dirname "{}")")_$(basename "{}")"'
find  -name knownResults.txt| xargs -I {} --verbose sh -c 'cp "{}" "$(basename "$(dirname "{}")")_$(basename "{}")"'
mv *known*html ./Known
mv *homer*html ./Denovo
mv *known*txt ./Table
cd ./Table
for i in *.txt; do awk '{print FILENAME"\t"$0}' $i >> CombinedTable.txt; done

#End notes: Load table in R and plot using ggplot2 (Homer_script.R)
