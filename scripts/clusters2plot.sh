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
cat $FILE | grep -v [DEBUG] | grep -v '^$' > c;

#split the files
awk "NR%$CSIZE==1{x=\"F\"++i;}{print > x}"  c;
rm c;

#plot_script
gnuplot < plot_script;

rm F*;
