#!/bin/bash

# Taryn Black, October 2018
# to run: ./downloadScenes.sh filename.ext
# filename.ext is a file containing a list of Landsat Product IDs (scene identifiers)

scenelist="$1"

cat $scenelist | while read id; do
    echo -e "\n===== DOWNLOADING $id ====="
    python3 getLandsat.py --productID=$id --basedir='/mnt/d/KEFJ_CoastalGlaciers/imagery/Landsat-L1/kefj' --epsg=6393
    echo -e "=============== DONE WITH SCENE ===============\n"
done

echo "Scene list download complete ($scenelist)."
#for id in 'cat $scenelist'; do
#    echo $id
#done
