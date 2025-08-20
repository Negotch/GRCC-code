# GRCC-code

## About

This is a small code example for the GRCC method in our paper [Uncertainty Quantification of Set-Membership Estimation in Control and Perception: Revisiting the Minimum Enclosing Ellipsoid](https://arxiv.org/abs/2311.15962)


## Preparation
1. Install Matlab
2. Clone this git repo (make sure to pull the submodules as well)
3. cd into spotless in Matlab, and run `spot_install.m` (we use spotless to define polynomial optimization problems)
4. Install [MOSEK](https://www.mosek.com/downloads/) (we use MOSEK as the SDP solver, make sure to get a free academic license)

## Get started
We provide two examples that are also included in the paper.

### Example 1:
In this example. We find the minimum enclosing ball of the TV Screen example.

Run `example_minimum_enclosing_ball.m'.

### Example 2:
In this example we show how we give error bounds for translation and rotation estimations in the perception example. For more details, please refer to the paper.

Run `example_purse.m'.

### Main function

- GRCC (See Example 1)
    - Input: 
        - problem: includes the defining variable and constraints (both equalities and inequalities)
        - kappa: the relaxation order of the Lasserre's Hierarchy
        - P: the projection matrix. The shape should be d by n, where d is the subspace dimension, n is the dimension of x
        - Q: the shape matrix of the ellipsoid. The shape should be d by d.
        - d: the subspace dimension that you care.
    - Output:
        - sol: a struct variable.
            - sol.upper_bound: the size of the ellipsoid
            - sol.a: the center of the ellipsoid
            - sol.Xopt: the moment matrix of the Lasserre's Hierarchy

