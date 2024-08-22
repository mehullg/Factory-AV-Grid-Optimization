# Factory-AV-Grid-Optimization

This repository contains a Julia implementation of a traffic intersection optimization model. The model schedules vehicles at intersections in a grid-based topology to minimize delays while adhering to constraints like speed limits, acceleration limits, and safety gaps.

## Table of Contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Running the Code](#running-the-code)
- [Explanation of the Code](#explanation-of-the-code)
  - [Parameters](#parameters)
  - [Functions](#functions)
- [Optimization Model](#optimization-model)
- [Results](#results)
- [Authors](#authors)

## Requirements
- Julia 1.6 or later
- The following Julia packages:
  - JuMP
  - HiGHS
  - DataFrames

## Installation
1. **Install Julia:** Download and install Julia from the official website.

2. **Install required packages:** Open Julia REPL and run the following commands:
    ```julia
    import Pkg
    Pkg.add("JuMP")
    Pkg.add("HiGHS")
    Pkg.add("DataFrames")
    ```

## Running the Code
1. **Clone the repository:**
    ```bash
    git clone <repository-url>
    cd <repository-folder>
    ```

2. **Open the Julia REPL:** Run `julia` in your terminal.

3. **Load the code:** Run the following commands in the Julia REPL to load and execute the code:
    ```julia
    include("factory_multi_intersection.jl")
    main(Φ, d, w1, w2, t0, v_avg, v_max, avg_a, max_a, grid_dim, L)
    ```
   Replace the parameters in the `main` function with your desired values.

4. **View Results:** The results will be printed in the REPL, and you can analyze them further by integrating them into a DataFrame or exporting them to a file.

## Explanation of the Code

### Parameters
The `main` function takes the following parameters:
- **Φ:** A list of directions for each vehicle in the grid. Possible directions include NORTH, EAST, SOUTH, and WEST.
- **d:** A list of initial distances of vehicles from their first intersection.
- **w1, w2:** Weights for different parts of the objective function.
- **t0:** Initial time.
- **v_avg:** Average speed of the vehicles.
- **v_max:** Maximum speed of the vehicles.
- **avg_a:** Average acceleration.
- **max_a:** Maximum acceleration.
- **grid_dim:** Dimensions of the grid (e.g., a 2x2 grid).
- **L:** Distance between intersections.

### Functions
- **custom_mod(x, n):** Helper function to handle modulo operations.
- **main():** The main function that sets up and solves the optimization problem.

## Optimization Model
The model is built using the JuMP package and solved using the HiGHS solver. The optimization problem is formulated to:
- **Minimize Delays:** The objective function minimizes the delay of each vehicle while crossing the intersections.
- **Enforce Speed and Acceleration Constraints:** The model ensures that vehicles do not exceed the speed and acceleration limits.
- **Ensure Safety Gaps:** The model includes constraints to maintain safety gaps between vehicles on the same movement direction.

## Results
The model outputs the optimal access times (`t_access`) for each vehicle at each intersection in the grid. These results can be used to schedule vehicle movements in a way that minimizes delays and maintains safety.

## Authors
- Mehull Girdhar [@mehullg](https://github.com/mehullg)
- Hans Van Rooij
- Benoît Legat [@blegat](https://github.com/blegat)
