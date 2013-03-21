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
cat "consSig" | grep -v [DEBUG] | grep -v '^$' > cons;

#split the files
awk "NR%30==1{x=\"F\"++i;}{print > x}" cluster;
awk "NR%30==1{y=\"T\"++j;}{print > y}" cons  ;
rm cluster;
rm cons;

let CNUM=80
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
set size 1,20
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
   set offsets 0,0,0,0
   
   set key bottom
    set offsets 0,0,0.5,0.5 
   plot 'T'.i using 1:2 title "conservation level ".i with lines lc rgb "blue"
    
   #plot 'F'.i lc rgb "red" title "cluster ".i
   
}
EOF
let level+=1

done


#gnuplot < plot_script;

rm F*;
rm T*;
