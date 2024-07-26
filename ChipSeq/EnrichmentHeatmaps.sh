conda activate deeptools #or run in docker
set=Name
bw1=Sample.bw
path1=/pathto/bigwigs/
computeMatrix reference-point -p 15 -R bedfile1.bed  bedfile2.txt.bed -S $path1$bw1 --referencePoint center -a 2500 -b 2500 --skipZeros --blackListFileName /pathto/FixBlacklist_ENCFF356LFX.bed -o $set".matrix"
plotHeatmap -m ./$set".matrix" --missingDataColor white --whatToShow 'plot, heatmap and colorbar' --heatmapHeight=14 --samplesLabel 'bigwig 1'  --averageType mean --plotType se --colorMap Blues --dpi 300 -o $set".pdf"
