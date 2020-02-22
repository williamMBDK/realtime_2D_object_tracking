./commands/universal_autonome_segmentation/compile
COUNT=0
FILES=./tests/natural_images/*/*
for file in $FILES
do
  ((COUNT += 1))
  ppm=${file%.jpg}".ppm"
  echo "$COUNT: testing $ppm ..."
  ./commands/getppm $file $ppm
  ./compiled/universal_autonome_segmentation/main 3 $ppm no_output_needed noprint >> out/universal_autonome_segmentation/results.txt
  rm $ppm
done
