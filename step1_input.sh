year=$(date +%Y)
month=$(date +%m)
echo "The input data comes from: "$year"-"$month
cat /ldfssz1/ST_BIGDATA/USER/st_bigdata/outdumpinfo/$year"-"$month > out.txt
