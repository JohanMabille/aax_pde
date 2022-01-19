The project contains one main file and 6 cpp + header files: <br> 

Payoff: allows to take a vector of underlying values and compute a payoff vector. By default, payoffs of call and put are implemented. <br>
Boundary Condition: creates the boundary conditions associated with the PDE. <br>
Coef Equation: compute the parameters /alpha, /beta, /gamma, /delta associated with the PDE. <br>
Matrix System: Build the matrix system X_t = AX_{t+1} + b associated with the PDE. <br>
Mesh: Solve the above matrix system at each iteration until price is found. <br>
Closed form: Compute the price of a call option according to the Black-Scholes closed form formula. <br> 

Main: compute the price of a call option and its greeks via finite difference (PDE) and via closed form with the test_mesh function. The user is more than welcome to change the inputs of the function (in the main()). It takes approximately 20 seconds with 3 months maturity to get the price and the greeks. The user can uncomment rows 48-49 to see the mesh. 
Remark: The gamma computation seems to have an issue as we end up with twice the value of the closed form gamma. We were not able to find the source of the issue. 
