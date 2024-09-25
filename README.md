It is developed for Improvement of Gaussian Basis Set Accuracy
Utilizing Linear Regression Technique for parameter optimization.
Starting with the Universal Basis (UB) of Gaussian-type orbitals (GTO).
Alignment with GRASP92 results.
The method, implemented in FORTRAN, enhances basis set accuracy, facilitating precise Correlation Calculations.
Additionally enables accurate calculation of Specific Mass Shift (SMS), Nuclear Mass Shift (NMS), and Field Shift (FS) components in Isotopic Shift (IS).
Algorithm Steps
-	Assign initial guess values to the parameters of the Gaussian basis function and calculate three key quantities 〈E〉, 〈r〉 and 〈1⁄r〉 for each orbital. 
-	Calculate the error for key quantities and fit the values to a linear regression curve.
-	Determine the slope of the linear regression curve and iterate until the global minima (Optimized) for the slope are found.
