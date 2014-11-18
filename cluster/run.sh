for ((i=0;i<$1;i+=1)) do
	./generate $i
	./apcluster similarity.txt pre.txt data/cluster$i.txt
done
./prepare
