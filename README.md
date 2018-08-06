# ThesisCppSolution
This is the full solution of the code developed for writing my Master's thesis in Mathematical Finance. 06/08/2018

The code includes:
  - A MRG32k3a random number generator.
  - A general purpose matrix class.
  - A Automatic Differentiation object for reverse mode.
  - Vasicek and Hull-White short rate models. Hull-White is fitted with a Nelder-Mead routine.
  - Hermitic spline interpolation class.
  - Code examples for the simulation of CVA of interest rate swaps and thus also all necessary finance objects to do so (derivative objects, CDS curve etc.)
 
Note: This library not built by a profesional and as such it is very specific in areas where being general is possible. When the code is oddly specific it should be pointed out with a comment.
