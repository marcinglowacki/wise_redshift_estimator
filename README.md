# wise_redshift_estimator
Code which uses the mid-infrared WISE information to estimate the redshift of radio sources

Created by Marcin Glowacki, University of Sydney.
Requirements: numpy, scipy, pylab, pandas

If used for your own publications, please cite the paper! "WISE data as a photometric redshift indicator for radio AGN", Glowacki et al. 2017, found at https://arxiv.org/abs/1709.08634 (submitted to MNRAS).

Simply copy files to your local directory and run as a python script. Example test data files/variables are provided, as are the outputs.

WISE_redshift_estimator.py is for a single radio source whose WISE colour information you input at the top of the file (W1_in, W1_er, etc for W1, W2 and W3 magnitudes).

WISE_redshift_estimator_mult.py is for multiple radio sources who WISE information is given in a .csv file; see test.csv file. Plots and statistics are given in the output/ directory and results.csv respectively. To adjust file type, location, etc., changes should be made in the relevant sections of the code. 
