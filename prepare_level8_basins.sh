
unzip -j shapes.zip
# how many rows?
ogrinfo -al  -geom=no hybas_lev00_all_8.shp | grep PFAF_8 | wc -l  
# row number is new ID
ogrinfo -al  -geom=no hybas_lev00_all_8.shp | grep PFAF_8  | awk '{ if (NR>1) print $NF , NR-1  }'  > HYBAS_ID_NROW.txt
# attach new ID to shapefile
oft-addattr-new.py  hybas_lev00_all_8.shp  PFAF_8   NROW_Int  Int   HYBAS_ID_NROW.txt
# rasterize based on new ID
gdal_rasterize -ot UInt32 -tap   -ot UInt32 -a_srs EPSG:4326 -l hybas_lev00_all_8  -a NROW_Int  -a_nodata 0  -tr  0.008333333333333 0.008333333333333 -co COMPRESS=LZW    hybas_lev00_all_8.shp    hybas_lev00_all_8.tif


