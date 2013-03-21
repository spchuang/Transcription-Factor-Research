#script file to split the clusters results to gnuplot readable format

if [ "$1" ]
then
   FILE=$1
else
  echo "No input file name"
fi

if [ "$2" ]
then
   CSIZE=$2
else
  echo "No cluster size input"
fi

#remove debug message
cat $FILE | grep -v [DEBUG] | grep -v '^$' > cluster;
cat "testing_flip_consSig" | grep -v [DEBUG] | grep -v '^$' > cons;
cat "testing_filp_temp" | grep -v [DEBUG] | grep -v '^$' > t;

#split the files
awk "NR%30==1{x=\"F\"++i;}{print > x}" cluster;
awk "NR%30==1{y=\"T\"++j;}{print > y}" cons  ;
awk "NR%30==1{z=\"A\"++k;}{print > z}" t  ;
rm cluster;
rm cons;
rm t;

let CNUM=2166
let DISPLAY=20
let level=0
for (( min=1; min<=CNUM; min+=DISPLAY ))
   do

let max=$min+$DISPLAY-1

if (( max > CNUM ));
#if [$max>$CNUM];
then
  max=$CNUM
fi
#plot_script

gnuplot << EOF

reset
set output "data$min-$max.ps"
set terminal postscript eps monochrome enhanced dashed
set style data linespoints
set xlabel "base pair"
set ylabel "signal"
set format y "" 
set size 2,20
set multiplot

#set size 2,20


#number of cluster
n=$max
mm=$min
do for [i=mm:n]{
   j=i-mm+1

   set size 1, 1
   set origin 0, (20-j)
   set key top
   plot 'F'.i using 1:2 title "DnaseI level".i with lines lc rgb "red" 
   set key bottom
   plot 'T'.i using 1:2 title "conservation level ".i with lines lc rgb "blue"
   
   #plot 'F'.i lc rgb "red" title "cluster ".i
   
   set size 1, 1
   set key top
   set origin 1, (20-j)
   plot 'A'.i lc rgb "red" title "anchor ".i
}
EOF
let level+=1

done


#gnuplot < plot_script;

#rm F*;
#rm T*;
#rm A*;
