#/usr/bin
for file in *_circ.txt;
	do filename=`echo $file | cut -d "." -f 1`;
	cat file | cut -f 1,2,3,15,14,13 > ${filename}_circ_small.txt	
	done;
exit
