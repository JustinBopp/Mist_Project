# Mist_Project

### To run:
Download the Mist File folder, Mist_Project.py and read_mist_models.py. Place both .py files into the folder. Run Mist_Project.py from the command line.

### Instructions:
##### 1st input: 
This the name of the files. Only put in the number in front of the M that you would like to open. This is the initial mass of the star. '00100' means an initial mass of 1.00 solar masses
##### 2nd input: 
Enter the process that you would like to run on the data. Options:\
Interpolation: random selects 50 datapoints from the file and interpolates 100 data points\
Integration: intergrates the Luminosity over time of both the main sequence and post-main sequence stars and prints the total energy for each\
Visualize: creates plots based on the style
##### 3rd input:
Enter the style you would like to perform. Options:\
Interpolation : 'linear', 'gpr',‘nearest’, ‘nearest-up’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’, ‘next’, ‘zero’, ‘slinear’, ‘quadratic’\
Integration : 'trapezoid' , 'simpson', 'midpoint'\
Visualize :  'compare','hr','core', 'radius', 'surface'
##### 4th input:
(default = star_age) The X value
##### 5th input:
(default = log_L) The Y value

