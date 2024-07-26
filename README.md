# NSD2_req_subunit
Code Repository for Paper : NSD2 is a requisite subunit of the AR/FOXA1 neo-enhanceosome in promoting prostate tumorigenesis

Paper Citation goes here. 


Repo Contents


RNA-Seq

  Run_RNAseq_kallisto.sh - basic pseudo alignment
  
  bulkRNA_Rscripts.R - R >3.6.0, additional analysis on bulk RNA data, GSEA, volcano outputs

  
  
ChIP-Seq 

  ChipAlignment_PE.sh - Basic alignment paired end pipeline (bwamem, samtools, picard wrapper)
  
  EnrichmentHeatmaps.sh - Uses Deeptools to create heatmaps
  
  OverlapCompare.R - Uses R>3.6.0 to implement a combination of ChipPeakAnno and ChipSeeker packages to compare output from MACS peak calling 
  
  Run_HOMER.sh - Runs HOMER motif analysis 
  
  TrackPlot.R - Uses R> 4.0.1 and R package Sushi to generate IGV style tracks with labeled motifs.
  

Copy of version information for software/tools needed to run these. Cites.  
