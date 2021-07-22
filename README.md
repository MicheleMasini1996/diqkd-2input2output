# diqkd-2input2output

The code contained in the two folders is based on the article "Simple and practical DIQKD security analysis via BB84-type uncertainty relations and Pauli correlation constraints" (https://arxiv.org/abs/2107.08894). The program language used is Julia.

The file contained in the folder two-basis allows to compute the key rate as a function of the channel error rate or the conditional entropy H(A|E) as a function of the CHSH statistic. To run the code, besides a recent Julia distribution, one needs to install the packages Optim, ForwardDiff, Roots, Plots, LinearAlgebra, DynamicPolynomials, SumOfSquares, MosekTools.

One can use the program, e.g., creating a Jupyter Notebook in the same folder where the file two-basis.jl is contained and typing

    include("two-basis.jl")
    
    rate(0.08,0.5,0.1) # computes the key rate when the channel error rate is 8%, the probability to use each of Alice's basis is 50% and q=10%
    
    plot_rate(0.6,0.3) # plots the key rate as a function of the channel error rate with 60% of probability of using A_1 and q=30%
    
    fzero(delta->rate(delta,0.5,0.49),0.0911) # finds the threshold channel error rate for 50% of probability of using A_1 and q=49%. It starts looking for the threshold close to delta=9.11%
