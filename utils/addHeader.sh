for i in *.[ch]
do
    cat head.txt $i >$i.new && mv $i.new $i
done
