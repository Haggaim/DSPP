# DS++: A flexible, scalable and provably tight relaxation for matching problems
Nadav Dym, Haggai Maron, Yaron Lipman
SIGGRAPH Asia 2017

Abstract
------------
Correspondence problems are often modelled as quadratic optimization problems over permutations. Common scalable methods for approximating solutions of these NP-hard problems are the spectral relaxation for non-convex energies and the doubly stochastic (DS) relaxation for convex energies. Lately, it has been demonstrated that semidefinite programming relaxations can have considerably improved accuracy at the price of a much higher computational cost. We present a convex quadratic programming relaxation which is provably stronger than both DS and spectral relaxations, with the same scalability as the DS relaxation. The derivation of the relaxation also naturally suggests a projection method for achieving meaningful integer solutions which improves upon the standard closest-permutation projection. Our method can be easily extended to optimization over doubly stochastic matrices, partial or injective matching, and problems with additional linear constraints. We employ recent advances in optimization of linear-assignment type problems to achieve an efficient algorithm for solving the convex relaxation. 
We present experiments indicating that our method is more accurate than local minimization or competing relaxations for non-convex problems. We successfully apply our algorithm to shape matching and to the problem of ordering images in a grid, obtaining results which compare favorably with state of the art methods. We believe our results indicate that our method should be considered the method of choice for quadratic optimization over permutations.

--------
This code supports matching n points on the first shape to k points on the second shape (k<=n).
Run demo.m for a quick demostration
questions regarding the code can be sent to nadavdym@gmail.com

Disclaimer:
------------
The code is provided as-is for academic use only and without any guarantees. Please contact the authors to report any bugs.


Code written by Haggai Maron and Nadav Dym. 
