for f in *_1.methylFreq
do
	g=`echo $f | sed "s/_1/_2/g"`
	n=`echo $f | sed "s/_1//g"`
	cat $f $g > $n
done
