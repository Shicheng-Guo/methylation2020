methyl2BED="/home/dinh/scripts/methylFreq2BED.pl"
for f in NC*methylFreq PC*methylFreq
do
	g=`echo $f | sed 's/.methylFreq//g'`
	$methyl2BED $g 10 < $f > $g.bed
done
