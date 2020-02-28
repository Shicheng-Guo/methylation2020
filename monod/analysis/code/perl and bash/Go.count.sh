for f in *bed
do
	echo $f 
	awk '{if($0 !~ /track/ && $0!~ /Lambda/){sum+=$5; cnt+=1;}} END {print sum"\t"cnt}' $f
	#awk '{if($0 ~ /Lambda/){sum+=$5; cnt+=$5*$4;}} END {print cnt"\t"sum}' $f
done
