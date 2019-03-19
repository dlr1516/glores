# glores - Globally Optimal Registration based on Fast Branch and Bound
#### Copyright (C) 2019 Luca Consolini, Mattia Laurini, Marco Locatelli, Dario Lodi Rizzini.

OVERVIEW
-------------------------------------------------

Library **glores** implements the globally optimal registraion method 
using two lower bound estimations. 
It has been kept to a minimal design. 

If you use this library, please cite the following paper: 

L. Consolini, M. Laurini, M. Locatelli, D. Lodi Rizzini. 
Globally Optimal Registration based on Fast Branch and Bound. 
under review, 
arXiv [arXiv:1901.09641](https://arxiv.org/abs/1901.09641).

````
@misc{consolini2019glores,
  author={Consolini, L. and Laurini, M. and Locatelli, M. and Lodi Rizzini, D.},
  title = {Globally Optimal Registration based on Fast Branch and Bound},
  Year = {2019},
  Eprint = {arXiv:1901.09641},
}
````

DEPENDENCIES
-------------------------------------------------

The software depends on the following external libraries

- Boost 
- Eigen 3.0 

If Matlab (we tested version R2018a=9.4) is available, it also compiles the glores_mex 
extension.
Other dependencies are placed in directories named thirdparty/. 

