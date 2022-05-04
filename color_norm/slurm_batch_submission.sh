FILES=./*.job
for f in $FILES;
do
	echo "processing" $f
	#sbatch $f;
done
