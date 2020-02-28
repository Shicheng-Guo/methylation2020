ls
R
ls
exit
R
ls
pwd
cd /mnt/
ls
cd bigdata/Genetic/Projects/
ls
cd shg047/
ls
cd rheumatology/RA/ASA/eqtl/GTEx_Analysis_v7_eQTL/
ls
ls -larth
pwd
data<-read.xlsx("/mnt/bigdata/Genetic/Projects/shg047/rheumatology/RA/ASA/eqtl/GTEx_Analysis_v7_eQTL/snpnexus_52995.xls",sheet=1)
R
ls
ll
ls -larth
qstat
exit
ll
ls
ll
exit
ls
exit
R
ls
rm DelayedMatrixStats_1.4.0.tar.gz 
ls
ln -s /mnt/bigdata/Genetic/Projects/shg047 hpc
ls
cd hpc/
ls
ll
cd rheumatology/
ls
cd RA/
ls
ll
cd NatureCommunication/
ls
ll
cd GSE112658/
ls
ll
R
top
ls -larth
top
ls -larth
cd hpc/rheumatology/RA/
ls
ll
ls
cd GS
l
ll
cd NatureCommunication/GSE112658/
ls
ll
ls
R
cd hpc/rheumatology/RA/NatureCommunication/GSE112658/
ls
cd WGBS/
cd ..
ll
R
SEQ<-seq(1,31121604,by=2000000)
DMRS<-c()
for(i in 1:(length(SEQ)-1)){
dmlTest.sm <- DMLtest(BSobj[SEQ[i]:SEQ[i+1],], group1=group1, group2=group2, smoothing=TRUE)
dmrs <- callDMR(dmlTest.sm, delta=0.1, p.threshold=0.05)
DMRS<-rbind(DMRS,dmrs)
print(i)
}
save.image("WGBS.DSS.RData")
ls
cd hpc/rheumatology/RA/NatureCommunication/
ls
ll
cd GSE112658/
ls
ll
R
ls
R
cd hpc/
cd rheumatology/RA/NatureCommunication/GSE112658/
ls
R
ls
ll
screen -s R
screen -S R
ls -larth
ll
screen
ls -larth
top
ls -larth
screen -ls
screen -R 28722.pts-0.birc10-lc
ls
ls -larth
cd
cd -
screen -R 28722.pts-0.birc10-lc
ls -larth
screen -R 28722.pts-0.birc10-lc
ll
q
ls -larth
cd /mnt/bigdata/Genetic/Projects/shg047/rheumatology/RA/NatureCommunication/GSE112658
ls -larth
less -S RA-Naturecommnunication.DMR.txt
ls -larth
less -S RA-Naturecommnunication.DMC.txt
ls -larth
less -S RA-Naturecommnunication.DMC.txt
ls -larth
wc -l RA-Naturecommnunication.DMR.txt
ls -larth
less -S RA-Naturecommnunication.DMR.txt
path
pwd
exit
ls
screen -R 28722.pts-0.birc10-lc
exit
screen -R 28722.pts-0.birc10-lc
ls
ll
ls
R
SEQ<-seq(1,31121604,by=20000)
DMRS<-c()
for(i in 1:(length(SEQ)-1)){
dmlTest.sm <- DMLtest(BSobj[SEQ[i]:SEQ[i+1],], group1=group1, group2=group2, smoothing=TRUE)
dmrs <- callDMR(dmlTest.sm, delta=0.1, p.threshold=0.05)
DMRS<-rbind(DMRS,dmrs)
print(i)
}
save.image("WGBS.DSS.RData")
R
ls
ls -larth
cd
exit
screen -R 28722.pts-0.birc10-lc
exit
R
ls
cd hpc/methylation/HCC/
ls
ll
R
cd hpc/methylation/HCC/
ls
ll
R
ll
pwd
R
ll
ls
cd
exit
R
ls
ll
exit
R
ls
ll
cd hpc/
ls
ll
ls
ll
vep
docker run -t -i ensemblorg/ensembl-vep ./vep
ls
cd
ll
cd /bin/
ls
cd
ls
ll
rm Tue_May_28_22_06_22_2019.tar.gz Tue_May_28_22_26_20_2019.tar.gz
ls
ll
cd GDCdata/
ls
ll
cd TCGA-READ
ls
ll
cd harmonized/
ls
cd Clinical/
ls
cd Clinical_Supplement/
ls
cd 388bc388-0c87-4000-817e-50bebc536696
ls
cd
ls
ll
cd ..
ls
ll
cd rfrahm/
ls
ll
cd ..
ll
cd ..
ls
cd MFLDCLIN/
ls
ll
wpd
pwd
ls
ls -larth
cd
cd hpc/tools/
ls
ll
ls -larth
wget http://org.gersteinlab.aloft.s3-website-us-east-1.amazonaws.com/aloft-1.0.zip
ls -larth
ls alo*
rm aloft-1.0.zip.1
cd aloft/
ls
ll
cd aloft-annotate/
ls
ll
./aloft 
pwd
vim ~/.bashrc 
source ~/.bashrc 
ls
cd hpc/
ls
ll
cd
ls
cd hpc/tools/
ls
ll
cd aloft/
ls
ll
cd aloft-annotate/
ls
ll
cd data/
ls
ll
unzip data.zip 
ls
ll
cd data/
ls
pwd
cd data_aloft_annotate/
ls
ll
less data.txt 
ll
cd hp
ll
exit
cd hpc/
cd db/
cd hg19/
ls
ls *tss
ls *tss*
ls
grep TSS *
ls 
dim(target)
docker
which docker
cd /usr/bin/docker
cd /usr/bin/
ls
cd ..
ll
cd bin/
ls
docker pull ensemblorg/ensembl-vep
ls
ll
ls
ll
exit
ls
cd hpc/
ls
cd alo
ls
ll
cd tools/
ls
cd aloft/
ls
cd aloft-annotate/
ls
./aloft 
python
ls
pwd
cd
cd hpc/
ls
cd autism/
ls
ll
cd da
cd data/
ls
ll
ls
cd aloft/
ls
ll
cd ..
ls
cd 2LOF/
ls
ll
cp All_samples_Exome_QC.DG.vcf ../aloft/
cd ../aloft/
ls
aloft
vim ~/.bashrc 
source ~/.bashrc 
cd hpc/
ls
ll
cd autism/
ls
cd data/
ls
cd aloft/
ls
aloft --vcf All_samples_Exome_QC.DG.vcf --output All_samples_Exome_QC.DG
ls
ll
aloft
pwd
ls
source ~/.bashrc 
R
aloft
pwd
ls
aloft --vcf All_samples_Exome_QC.DG.vcf --output All_samples_Exome_QC.DG 
cd
cd hpc/tools/
ls
ll
cd aloft/
ls
cd aloft-annotate/
ls
ll
cd ..
ls
cd aloft-predict/
ls
cd ..
ls
cd aloft-
ls
cd aloft-annotate/
ls
ll
mkdir data
cd data/
ls
ll
wget http://org.gersteinlab.aloft.s3-website-us-east-1.amazonaws.com/data.zip
ls
unzip data.zip 
ls
cd data/
ls
cd ..
ls
ll
cd data/
ls
cd data/
ls
cd ..
ls
mv data ../
ls
cd ..
ls
ll
cd data/
ls
cd data/
ls
mv * ../
cd ..
ls
ll
rm -rf data
ls
ll
pwd
ls
pwd
aloft
pwd
ls
cd data_aloft_
ls
cd data_aloft_annotate/
ls
less data.txt 
pwd
less data.txt 
ls
ll
cd
cd hpc/autism/
ls
cd data/
ls
cd loa
ls
cd aloft/
ls
ll
cd All_samples_Exome_QC.DG/
ls
ll
less snpinput_temp 
ls -larth
cd ..
ls
ls -larth
cd cache/
ls
ls -larth
cd ..
ls
ls -larth
cd All_samples_Exome_QC.DG/
ls
ll
ls
ls -larth
less bigwig.tab 
ls -larth
less All_samples_Exome_QC.DG.vcf.vat
less -S All_samples_Exome_QC.DG.vcf.vat
ls
ll
cd ..
ls
ll
wc -l All_samples_Exome_QC.DG.vcf
python
pwd
cd
cd hpc/tools/
ls
pwd
cd aloft/
cd aloft-annotate/
ls
pwd
cd
cd hpc/autism/data/aloft/
ls
wc -l All_samples_Exome_QC.DG.vcf 
cd All_samples_Exome_QC.DG/
ls
ll
less -S All_samples_Exome_QC.DG.vcf.vat.aloft.lof
wc -l All_samples_Exome_QC.DG.vcf.vat.aloft.vcf
wc -l *
ls
cd ..
ls
cd All_samples_Exome_QC.DG/
ls
ll
wc -l All_samples_Exome_QC.DG.vcf.vat
less -S All_samples_Exome_QC.DG.vcf.vat
ls
mkdir Gnomad
pwd
cd ~/hpc/db/Gnomad/vcf
ls
ll
ls
cd
cd hpc/tools/aloft/temp/
ls
aloft
aloft --vcf gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop.vcf --data 
aloft --vcf gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop.vcf --data ../aloft-annotate/data/data_aloft_annotate/data.txt --output gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop
cd ../aloft-annotate/data/data_aloft_annotate/
ls
pwd
cd -
ls
ll
aloft --vcf gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop.vcf --data /home/local/MFLDCLIN/guosa/hpc/tools/aloft/aloft-annotate/data/data_aloft_annotate/data.txt  --output gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop 
aloft --vcf gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop.vcf --data /home/local/MFLDCLIN/guosa/hpc/tools/aloft/aloft-annotate/data/data_aloft_annotate/  --output gnomad.exomes.r2.1.sites.chr21.rec.ExomeStop 
ls
ll
exit
aloft
ls
aloft --vcf All_samples_Exome_QC.DG.vcf --output All_samples_Exome_QC.DG
aloft --vcf All_samples_Exome_QC.DG.vcf --output All_samples_Exome_QC.DG --data /home/local/MFLDCLIN/guosa/hpc/tools/aloft/aloft-annotate/data/data_aloft_annotate
ls
ls -larth
cd All_samples_Exome_QC.DG/
ls
ll
cd ..
ls
ll
cd All_samples_Exome_QC.DG/
ls
ll
wc -l All_samples_Exome_QC.DG.vcf.vat.aloft.lof
less -S All_samples_Exome_QC.DG.vcf.vat.aloft.lof
ls
ll
cd ..
cd
cd hpc/tools/aloft/
cd aloft-predict/
ls
./aloft_predict.pl 
ls
cd /home/local/MFLDCLIN/guosa/hpc/autism/data/aloft/All_samples_Exome_QC.DG
cd
ls
cd hpc/tools/aloft/
ls
cd aloft-annotate/
ls
ll
cd data/
ls
wget http://org.gersteinlab.aloft.s3-website-us-east-1.amazonaws.com/hg19.zip &
ls -larth
unzip hg19.zip 
ls -larth
cd hg19/
ls
ll
gunzip *.gz
ls -larth
less -S chr1.vcf.vat.aloft.lof.predict
ls -larth
less -S chr15.vcf.vat.aloft.lof.predict
wc -l *predict
less aloft_scores.README 
ls -larth
less -S chr1.vcf.vat.aloft.lof.predict
R
ls
ll
less -S chr1.vcf.vat.aloft.lof.predict 
vim chr22.vcf.vat.aloft.lof.predict 
R
ls
ls -larth
perl -p -i -e 's/chr//i' aloft.hg19.txt 
ls -larth
less aloft.hg19.txt 
vim aloft.hg19.txt
head aloft.hg19.txt
perl -p -i -e 's/Recessive/Rec/i' aloft.hg19.txt 
ls -larth
less aloft.hg19.txt 
perl -p -i -e 's/Dominant/Dom/i' aloft.hg19.txt 
ls -larth
less aloft.hg19.txt 
gzip aloft.hg19.txt.gz aloft.hg19.txt
ls -larth
less aloft.hg19.txt.gz
ls -larth
less aloft.hg19.txt.gz
tail aloft.hg19.txt.gz
PuTTYPuTTYPuTTYPuTTYPuTTYPuTTYPuTTYPuTTYPuTTYPuTTYPuTTY
ls
ll
R
l
ll
cd
cd hpc/autism/
ls
cd data/
ls
ll
ls
cd aloft/
ls
ll
cd ..
ls
cd 2lof/
ls
ll
cd beagle/
ls
less -S All_samples_Exome_QC.temp.vcf.recode.clean.chr22.vcf.recode.vcf
ls -larth
less -S All_samples_Exome_QC.chr2.vcf.log
less -S All_samples_Exome_QC.chr2.vcf
less -S All_samples_Exome_QC.chr2.vcf.vcf.gz
R
ls
exit
vim ~/.bashrc 
ls
ll
ls
ll
./aloft
lsb_release -a
lsb_release
ls
ll
cd ..
aloft
cd
cd hpc/
ls
cd autism/
ls
cd data/
ls
ll
ls
cd ..
ls
ll
cd data/
ls
ls -larth
mkdir aloft
cd 2LOF/
ls
ll
ls -larth
rm All_samples_Exome_QC.DG.vcf.gz
aloft 
aloft --vcf All_samples_Exome_QC.DG.vcf --output All_samples_Exome_QC.DG 
cd ..
ll
ls
cd
cd hpc/tools/
ls
ll
cd aloft/
ls
cd aloft-
cd aloft-annotate/
ls
ll
./aloft
ls
ll
mkidr data
mkdir data
ls
cd data/
ls
ll
wget http://org.gersteinlab.aloft.s3-website-us-east-1.amazonaws.com/data.zip
R
exit
docker
docker --info
docker pull broadinstitute/gatk
exit
docker pull broadinstitute/gatk
docker
docker images
docker pull sevenbridges/maftools
docker pull ensemblorg/ensembl-vep
exit
R
ll
cd hpc/
ls
cd project/
cd LungBrainMetastasis/
ls
ll
ls -larth
cd vcf/
ls
ll
pwd
ls
cd annovar/
ls
ll
grep ERF *.txt
ll
ls -larth
wc -l *.txt
ls -larth
ls
ll
cd
cd hpc/
ls
cd methylation/
ls
ll
ls
cd ..
ls
ll
cd tcga/
ls
ll
cd ..
cd methylation/
cd Pancancer/
ls
ll
cd RNA-seq/
ls
cd ..
ls
cd ..
ls
ll
ls *BRCA*
ls -l *BRCA*
pwd
R
ls
ls -l *BRCA*
R
ll
cd
R
ll
docker ps
exit
R
ll
docker
docker ps
exit
docker
docker ps
cd
ls
cd vep_data/
ls
ll
cd ..
ll
rm Wu_Target.txt 
ll
pwd
exit
docker ps
docker ps STATUS
docker login
docker run broadinstitute/gatk
docker ps
docker run ensembl-vep
docker pull ensemblorg/ensembl-vep
docker run ensemblorg/ensembl-vep
vep
docker run ensemblorg/ensembl-vep
ls
./vep
docker run -it ensemblorg/ensembl-vep ./vep 
./vep
docker run -it ensemblorg/ensembl-vep ./vep 
ll
cd hpc/
ls
cd project/
cd LungBrainMetastasis/
ls
ll
ls
cd vcf/
ls
docker run -it ensemblorg/ensembl-vep ./vep -i F.T2.vcf 
less -S F.T2.vcf
mkdir $HOME/vep_data
chmod a+rwx $HOME/vep_data    
docker run -t -i -v $HOME/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep perl INSTALL.pl
docker run -t -i -v $HOME/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep perl INSTALL.pl -a cf -s homo_sapiens -y GRCh38
docker run -t -i -v $HOME/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep perl INSTALL.pl -a cf -s homo_sapiens -y GRCh19
ifconfig
exit
R
cd hpc/
cd methylation/
cd Ingrid/
ls
cd MCaldwell-Sept27-17-HuMethEPIC
ls
cd Raw_Data/
ls
ll
cd idat/
ls
ll
R
ls
ll
ls -larth
R
cd hpc/
cd methylation/Ingrid/
ls
cd MCaldwell-Sept27-17-HuMethEPIC/
ls
cd Raw_Data/
ls
cd idat/
ls
pwd
R
ll
less ingrid.csv 
R
ls
R
ls
cd hpc/methylation/Ingrid/
ls
cd MCaldwell-Sept27-17-HuMethEPIC/
ls
cd Raw_Data/
ls
ll
ls
cd idat/
ls
R
