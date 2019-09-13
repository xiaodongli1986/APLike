for nmu in 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
do
	for iz in 2 3 4 5 6
	do
		file1=${nmu}mubins.mumax0.97.iz$iz.CovMock_1900.s6.0to40.0.NBComp.covmat
		cp /home/xiaodongli/software/APLike/covmat_files/$file1 ../covmat_files
	done
done

for nowfile in /home/xiaodongli/software/APLike//2pcfs/DR12_LOWZ_J08.RSD.000.xyzw.1of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_LOWZ_J08.RSD.001.xyzw.1of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_LOWZ_J08.RSD.002.xyzw.1of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_LOWZ_J08.RSD.003.xyzw.1of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_LOWZ_J08.RSD.000.xyzw.2of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_LOWZ_J08.RSD.001.xyzw.2of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_LOWZ_J08.RSD.002.xyzw.2of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_LOWZ_J08.RSD.003.xyzw.2of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_LOWZ_J08.RSD.000.xyzw.3of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_LOWZ_J08.RSD.001.xyzw.3of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_LOWZ_J08.RSD.002.xyzw.3of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_J08.RSD.000.xyzw.1of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_J08.RSD.001.xyzw.1of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_J08.RSD.002.xyzw.1of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_J08.RSD.003.xyzw.1of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_J08.RSD.000.xyzw.2of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_J08.RSD.001.xyzw.2of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_J08.RSD.002.xyzw.2of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_J08.RSD.003.xyzw.2of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_J08.RSD.000.xyzw.3of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_J08.RSD.001.xyzw.3of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_J08.RSD.002.xyzw.3of3.rmax150.150rbins.120mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_J08.RSD.003.xyzw.3of3.rmax150.150rbins.120mubins.2pcf
do
cp $nowfile ../2pcfs
done
exit

for nowfile in /home/xiaodongli/software/APLike//2pcfs/DR12_LOWZ_data.xyzw.1of3.cosmo-converted.om0.2600_w-1.0000.rmax150.750rbins.600mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_LOWZ_data.xyzw.2of3.cosmo-converted.om0.2600_w-1.0000.rmax150.750rbins.600mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_LOWZ_data.xyzw.3of3.cosmo-converted.om0.2600_w-1.0000.rmax150.750rbins.600mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_data.xyzw.1of3.cosmo-converted.om0.2600_w-1.0000.rmax150.750rbins.600mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_data.xyzw.2of3.cosmo-converted.om0.2600_w-1.0000.rmax150.750rbins.600mubins.2pcf /home/xiaodongli/software/APLike//2pcfs/DR12_CMASS_data.xyzw.3of3.cosmo-converted.om0.2600_w-1.0000.rmax150.750rbins.600mubins.2pcf 
do
cp $nowfile ../2pcfs
done
 
 
