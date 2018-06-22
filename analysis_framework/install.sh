########################################################################
# File Name: install.sh
# Author: Mohammad Hattawy
# mail: mohammad.hattawy@gmail.com
# Created Time: July 3rd, 2017
#########################################################################
#!/bin/bash
cd convert/convert_caen/
cp Bin2TxtConv_CAEN.cc $ROOTDEV/bin
cd ../../analysis
make clean
make
make install
cp SplineFitDB.rdb $ROOTDEV/bin
cd ../paf
cp paf $ROOTDEV/bin
cp paf_agilent $ROOTDEV/bin
cd ../

echo " -> lappd is copied to $ROOTDEV/bin"
echo " -> SplineFitDB.rdb is copied to $ROOTDEV/bin"
echo " -> Bin2TxtConv_CAEN.cc is copied to $ROOTDEV/bin"
echo " -> paf is copied to $ROOTDEV/bin"
echo " -> paf_agilent is copied to $ROOTDEV/bin"
echo " ---> Installation done!"

