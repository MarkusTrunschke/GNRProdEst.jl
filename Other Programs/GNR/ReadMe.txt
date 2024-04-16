These files contain the programs and data for replicating the tables and figures in: 
Amit Gandhi, Salvador Navarro, and David A. Rivers, "On the Identification of Gross Output Production Functions", Forthcoming, Journal of Political Economy.
The files are described here, organized by the associated table/figure.
Figure 1:
	Monte_Carlo_Code: This folder contains the Fortran code files used to generate our Monte Carlo simulated data. As documented in globvar.f90 and main.f90, one needs to change the value of "var_pr" depending on the degree of time series variation used (in globvar.f90), and rename the corresponding output file accordingly (in main.f90). While we provide all of the auxiliary routines needed, the code assumes that both MPI (since the code is implemented in parallel) and the Intel Math Kernel library (and associated Lapack and Blas packages) are installed. 
	Monte_Carlo_Data: This folder contains Stata (.dta) versions of the Monte Carlo simulated data, produced using the Monte Carlo code described above, and used to produce Figure 1.
	Estimation_Code: This folder contains the Stata code needed to produce the numbers used in Figure 1. Running "compile_all.do" runs all of the necessary routines.

Table 1:
	Monte_Carlo_Code: This folder contains the Fortran code files used to generate our Monte Carlo simulated data. There are three separate sub-folders, one corresponding to each parametric functional form specification (CD--Cobb-Douglas; CES; and TL--Translog). While we provide all of the auxiliary routines needed, the code assumes that both MPI (since the code is implemented in parallel) and the Intel Math Kernel library (and associated Lapack and Blas packages) are installed. 
	Monte_Carlo_Data: This folder contains text file (.out) versions of the Monte Carlo simulated data, produced using the Monte Carlo code described above, and used to produce Table 1.
	Estimation_Code: This folder contains the Stata code needed to produce the numbers used in Table 1. Running "our_estimator.do" runs all of the necessary routines.

Tables 2 and 3: This folder contains Stata code and Chilean and Colombian data used to produce Tables 2 and 3. There is one folder for each country. These two folders are further subdivided by industry. Within each industry sub-folder, there is the dataset in Stata (.dta) format, as well as folders containing the Stata code used for estimation (one for OLS and one for our procedure, labeled GNR). Running boot_compile.do (for OLS) and call_boot.do (for GNR) runs all the necessary routines for computing point estimates and standard errors.


