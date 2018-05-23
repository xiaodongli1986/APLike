

for catname in CMASS LOWZ
do
for ibin in 1 2 3
do
file1=/home/xiaodongli/SparseFilaments/data/input/boss2pcf/data/DR12v4-${catname}/xyzw.binsplitted/data.xyzw.${ibin}of3.cosmo-converted.om0.2600_w-1.0000.rmax150.750rbins.600mubins.2pcf
file2=../2pcfs/DR12_${catname}_data.xyzw.${ibin}of3.cosmo-converted.om0.2600_w-1.0000.rmax150.750rbins.600mubins.2pcf
cp $file1 $file2
ls $file1 -alh
ls $file2 -alh
for nummock in 0 1 2 3 
do
file1=/home/xiaodongli/SparseFilaments/data/input/boss2pcf/data/DR12v4-$catname/xyzw.binsplitted/J08.RSD.00${nummock}.xyzw.${ibin}of3.rmax150.150rbins.120mubins.2pcf
file2=../2pcfs/DR12_${catname}_J08.RSD.00$nummock.xyzw.${ibin}of3.rmax150.150rbins.120mubins.2pcf
cp $file1 $file2
done
done
done


#cp base*.txt /home/xiaodongli/software/cosmomcs/12Oct_generic/chains/wcdm/
#mv base*.txt nuisance
#mv base*info nuisance
