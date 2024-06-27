cd ../work_dir/dataset/raw/0
for i in {0..80}; do

cd ../$i
mpirun -np 64 /home/t.hsu/openmx/source/openmx openmx_in.dat > openmx.std
cat openmx.out >> openmx.scfout
rm -r openmx_rst *.cube

done
