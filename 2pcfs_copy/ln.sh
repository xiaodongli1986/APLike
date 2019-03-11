
for cat in LOWZ CMASS
do
for ibin in 1 2 3 
do
#for cos in om0.2600_w-1.0000 
for cos in om0.3100_w-1.0000 om0.3100_w-1.4000 om0.3100_w-0.6000 om0.4100_w-1.0000 om0.2600_w-1.5000 om0.2600_w-1.4000  om0.2600_w-0.6000 om0.2600_w-0.5000 om0.1100_w-1.0000 om0.1100_w-2.0000
do
	file1=/home/xiaodongli/LSS/boss2pcf/data/DR12v4-$cat/xyzw.binsplitted/data.xyzw.${ibin}of3.cosmo-converted.$cos.rmax150.750rbins.600mubins.2pcf
	file2=DR12_${cat}_data.xyzw.${ibin}of3.cosmo-converted.$cos.rmax150.750rbins.600mubins.2pcf
	#file1=/home/xiaodongli/SparseFilaments/data/input/boss2pcf/data/DR12v4-$cat/xyzw.binsplitted/data.xyzw.${ibin}of3.cosmo-converted.$cos.rmax150.150rbins.120mubins.2pcf
	#file2=DR12_${cat}_data.xyzw.${ibin}of3.cosmo-converted.$cos.rmax150.150rbins.120mubins.2pcf
        #ls $file1
        #echo $file2
        #rm $file2
	ln -s $file1 $file2
done
done
done
