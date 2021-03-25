Analysis and data pertaining to Laumont et al. 2021 (Clin Cancer Res)

R scripts used for data analysis:
- 20210316_5GEXpreprocessing.R: used to process CD45+ 5'GEX data from quality control to cell type annotation.
- 20210316_5GEXtcell.R: used to process 5'GEX data, focusing on T cells.
- 20210316_TCRtcell.R: used to analyze TCR data

data folder contains:
- processed single-cell sequencing data were generated with cell ranger 3.1.0 (cellranger count,  cell ranger vdj).
  Each folder (cd39-01-s[0-9]) corresponds to an HGSC patient (simply referred to as s[0-9] in the manuscript).
  For each patient, 5'GEX data are stored in the 5gex folder.
  For each patient, TCR data are stored in the tcr folder. 
- annotation files used throughtout the workflow are all located in the folder annotationFiles.
  samples.txt: file path to data files
  riboGenes.txt: file used to annotate ribosomal genes
  mart_export_geneLevel_GRCh38.p12_Jun2020.txt: file used to annotate IG/TR variable genes    
  markers_ENSEMBL_GRCh38.93.txt: file containing the marker genes used to annotate cell types with garnett

results folder contains:
- 20210316_preprocessed.rds: R data file saved at the end of the 20210316_5GEXpreprocessing.R script 
- 20210316_tcells.rds: R data file saved at the end of the 20210316_5GEXtcell.R script
