## Parallel Matrix Condensation for Calculating Log-Determinant of Large Matrix  
This repository contains the source code, sample input, output and result for following parallel algorithms that calculate the determinant of large dense matrix.  All algorithms are implemented by Fortran and are parallelized by MPI.

* Matrix Condensation (mc_p.f95)
* Gaussian Elimination (ge_p.f95)
* [Gaussian Elimination by LU decoposition of ScaLAPACK (ge_scalapack.f)](https://www.ibm.com/support/knowledgecenter/en/SSNR5K_4.2.0/com.ibm.cluster.pessl.v4r2.pssl100.doc/am6gr_lgetrf.htm)


### Prerequisites
* gfortran
* OpenMPI
* ScaLAPACK


### Instruction of Execution

#### Compile the source code:
```
make mc_p
make ge_p
make ge_scalapack
```

#### Run the executable:
Calculate log(abs(det)) of random 1000x1000 dense matrix using 4 processors by matrix condensation:
```
mpirun -n 4 mc_p 1000 
```
Calculate log(abs(det)) of a specific 1000x1000 dense matrix using 4 processors by matrix condensation:
```
mpirun -n 4 mc_p 1000 ./input/m1000x1000.bin 
```


### Output
The ouput folder is the raw experimental results from [OSCER cluster(OU Supercomputing Center for Education&Research)](http://www.ou.edu/oscer/resources/hpc). If you want to run the experiment in OSCER too, you can use the batchrun.sh and comp.sbatch provided in this repository.


### Results
The results folder contains the following:
* Runtimeresults.xlsx: spreadsheet that stores the average run time, speedup, communication time and data distribution time for all algorithms.
* Extract_results_plot: jupyter notebook that extracts the raw output, analyze data and generate the plots.

### Showcase
#### Execution time for three algorithms:
![Matrix Condensation Execution Time](https://github.com/vbvg2008/MatrixCondensation/blob/master/images/MC_average_time.png)

![Gaussian Elimination Execution Time](https://github.com/vbvg2008/MatrixCondensation/blob/master/images/GE_average_time.png)

![Gaussian Elimination Scalapack Execution Time](https://github.com/vbvg2008/MatrixCondensation/blob/master/images/GE_scalapack_average_time.png)

#### Speed-up for all problem sizes:
![Speed-up for all sizes](https://github.com/vbvg2008/MatrixCondensation/blob/master/images/Speedup_allsizes.png)

#### Average speed-up of all problem sizes:
![Average speed-up](https://github.com/vbvg2008/MatrixCondensation/blob/master/images/Average_speedup.png)

#### Average CPU communication time of all problem sizes(for mc_p and ge_p):
![Average communication](https://github.com/vbvg2008/MatrixCondensation/blob/master/images/Average_communication.png)

#### Average data distribution time of all problem sizes(for mc_p and ge_p):
![Average distribution](https://github.com/vbvg2008/MatrixCondensation/blob/master/images/Average_distribution.png)


### Contact:
Xiaomeng.dong-1@ou.edu
