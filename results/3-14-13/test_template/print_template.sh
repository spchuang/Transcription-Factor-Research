#script file to split the clusters results to gnuplot readable format

CSIZE=30
#remove debug message
cat "testing_filp_temp" | grep -v [DEBUG] | grep -v '^$' > t;

#split the files
awk "NR%30==1{z=\"A\"++k;}{print > z}" t  ;

rm t;


#plot_script

gnuplot << EOF

reset
set output "data_template.ps"
set terminal postscript eps monochrome enhanced dashed
set style data linespoints
set xlabel "base pair"
set ylabel "signal"
#set format y "" 
set size 10,10
set multiplot


#number of cluster

do for [i=1:99]{

   set size 1, 1
   set origin (i/10), (i%10)
   set key top
   plot 'A'.i lc rgb "red" title "anchor ".i
   MAX=GPVAL_Y_MAX
    MIN=GPVAL_Y_MIN
    set yrange [0:19]
}
EOF



#gnuplot < plot_script;
rm A*;
