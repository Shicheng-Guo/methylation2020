

http://smithlab.usc.edu/methbase/data/Thompson-Human-2015/Human_PancreaticCancer8/tracks_hg19/

# PBMC and Blood cells
wget http://smithlab.usc.edu/methbase/data/Hodges-Human-2011/Human_BCell/tracks_hg19/Human_BCell.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Hodges-Human-2011/Human_CD133HSC/tracks_hg19/Human_CD133HSC.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Hodges-Human-2011/Human_HSPC/tracks_hg19/Human_HSPC.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Hodges-Human-2011/Human_Neut/tracks_hg19/Human_Neut.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_PBMC/tracks_hg19/Human_PBMC.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_CD4T-100yr/tracks_hg19/Human_CD4T-100yr.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Heyn-Human-NewbornCentenarian-2012/Human_CD4T-Newborn/tracks_hg19/Human_CD4T-Newborn.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Li-PBMC-2010/Human_PBMC/tracks_hg19/Human_PBMC.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Heyn-Human-DNMT3BMut-2012/Human_BCell-Healthy/tracks_hg19/Human_BCell-Healthy.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Roadmap-Human-2015/Human_Macrophage/tracks_hg19/Human_Macrophage.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Roadmap-Human-2015/Human_Tcell/tracks_hg19/Human_Tcell.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Roadmap-Human-2015/Human_NK/tracks_hg19/Human_NK.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Roadmap-Human-2015/Human_HSC/tracks_hg19/Human_HSC.meth.bw & 
# Lung
http://smithlab.usc.edu/methbase/data/Xie-Human-2013/Human_IMR90/tracks_hg19/Human_IMR90.meth.bw
https://www.genboree.org/EdaccData/Current-Release/experiment-sample/Bisulfite-Seq/Lung/UCSD.Lung.Bisulfite-Seq.STL002.wig.gz
https://www.genboree.org/EdaccData/Current-Release/experiment-sample/Bisulfite-Seq/IMR90_Cell_Line/UCSD.IMR90.Bisulfite-Seq.combined.wig.gz

# Colon
wget http://smithlab.usc.edu/methbase/data/Berman-Human-2012/Human_ColonCancer/tracks_hg19/Human_ColonCancer.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Ziller-Human-2013/Human_Colon_Tumor_Primary/tracks_hg19/Human_Colon_Tumor_Primary.meth.bw &
wget http://smithlab.usc.edu/methbase/data/Hansen-Human-2011/Human_ColonCancer/tracks_hg19/Human_ColonCancer.meth.bw &

# Total 45 samples 
https://www.genboree.org/EdaccData/Current-Release/experiment-sample/Bisulfite-Seq/Adipose_Tissue/UCSD.Adipose_Tissue.Bisulfite-Seq.STL003.wig.gz

# bedgraph to bigwig
cd /home/shg047/oasis/Estellar2016/MF


# bigWigToBedGraph

for i in `ls *bw`
do
bigWigToBedGraph $i $i.bg &
done

# bigWigAverageOverBed

for i in `ls *bw`
do
bigWigAverageOverBed $i WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.renew.bed $i.aob &
done

bigWigAverageOverBed Human_CD4T-Newborn.meth.bw WGBS_pooled_mappable_bins.all_autosomes.mld_blocks_r2-0.5.bed Human_CD4T-Newborn.meth.bw.mhb.bg
