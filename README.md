Summary

This set of Matlab codes reproduce the figures and experimental results published in our paper:
Holographic Phase Retrieval and Optimal Reference Design.
David A. Barmherzig, Ju Sun, Emmanuel J. Candes, T.J. Lane, and Po-Nan Li. 

•	Folder ref_scaling to reproduce Figure 7.

•	Folder alg_compar to reproduce Figure 8.

•	Folder mimi_compar to reproduce Figure 9.

•	Folder flat_mi_compar to reproduce Figure 10.

•	Folder ref_deconv contains code to run the Referenced Deconvolution algorithm and compute the expected error for a given specimen and reference. The code is contains a special, optimized implementation if the reference is specified to be a pinhole, slit, or block.
  
  o	Ref_Deconv_fcn is a function implementing the Referenced Deconvolution algorithm.
  
  o	Ref_Deconv_tester is a sample file implementing Ref_Deconv_fcn.
  
  o	ref2mtrx is a function which generates the matrix M_R for any given reference R
  
  o	img_recov is a function which recovers a given specimen from the top-left quadrant of the cross-correlation of the specimen and a given reference.

•	Folder HIO contains code to run the HIO algorithm
  
  o	HIO_fcn is a function implementing the HIO algorithm.
  
  o	HIO_tester is a sample file implementing HIO_fcn.

•	Folder ref_HIO contains code to run the HIO algorithm, with the reference pixels enforced as an additional constraint (i.e. an additional projection)
  
  o	Ref_HIO_fcn is a function implementing HIO with the reference enforced.
  
  o	Ref_HIO_tester is a sample file implementing Ref_HIO_fcn.

Codes written by David Barmherzig and Ju Sun. For questions or bug reports please send email to David Barmherzig, davidbar@stanford.edu
Thanks to bug reporters:

