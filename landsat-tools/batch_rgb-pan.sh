#!/bin/bash

# Taryn Black, June 2020
# to run: ./batch_rgb-pan.sh filename.ext
# filename.ext is a file containing a list of Landsat Product IDs (scene identifiers)

scenelist="$1"

cat $scenelist | while read id; do
    echo -e "\nProcessing $id..."
    python3 composite_pansharpen.py --id=$id --path='/mnt/d/KEFJ_CoastalGlaciers/imagery/Landsat-L1/kefj/september/'
done

echo "Image processing complete ($scenelist)."
