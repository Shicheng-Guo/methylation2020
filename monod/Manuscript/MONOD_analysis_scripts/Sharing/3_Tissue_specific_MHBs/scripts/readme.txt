# Copy the code to Genome-miner to avoid code or input lost again
# data input 
WGBS_methHap_load_matrix_20Oct2015.txt   # GWBS MHL matrix phase-I
newsaminfo.txt                           # samplesheet
tissue2Layer.txt                         # germ-layer information
Heatmap.MHL.txt                          # other four cancer samples phase-II
# script
heatmap3.R                               # to label tissue type
heatmap4.R                               # to label germ-layer type
cd /home/shg047/oasis/monod/heatmap
# label the sample with tissue types
Rscript --vanilla heatmap3.R 1000  &
Rscript --vanilla heatmap3.R 3000 &
Rscript --vanilla heatmap3.R 5000 &
Rscript --vanilla heatmap3.R 8000 &
# label the sample with germ layer
Rscript --vanilla heatmap4.R 1000 &
Rscript --vanilla heatmap4.R 3000 &
Rscript --vanilla heatmap4.R 5000 &
Rscript --vanilla heatmap4.R 8000 &

# compare the similiary of the cluster
