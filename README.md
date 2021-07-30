# diqkd-2input2output

The code contained in the two folders is based on the article "Simple and practical DIQKD security analysis via BB84-type uncertainty relations and Pauli correlation constraints" (https://arxiv.org/abs/2107.08894). The program language used is Julia.

The file contained in the folder two-bases allows to compute the key rate as a function of the channel error rate or the conditional entropy H(A|E) as a function of the CHSH statistic when both of Alice's bases are exploited to generate the secret key. To run the code, besides a recent Julia distribution, one needs to install the packages Optim, ForwardDiff, Roots, Plots, LinearAlgebra, DynamicPolynomials, SumOfSquares, MosekTools. Notice that, to solve problems for which the probability of Alice using her basis is different from 0.5 or 1, you will further need a Mosek license in order to solve the Lasserre hierarchy.

The installation of the packages can be done by launching julia and typing

    using Pkg
    Pkg.add(["Optim", "ForwardDiff", "Roots", "Plots", "LinearAlgebra", "DynamicPolynomials", "SumOfSquares", "MosekTools"])

One can use the program, e.g., creating a Jupyter Notebook in the same folder where the file two-bases.jl is contained and typing

    include("two-bases.jl")
    
    rate(0.08,0.5,0.1) # computes the key rate when the channel error rate is 8%, the probability to use each of Alice's bases is 50% and q=10%
    
    plot_rate(0.6,0.3) # plots the key rate as a function of the channel error rate with 60% of probability of using A_1 and q=30%
    
    fzero(delta->rate(delta,0.5,0.49),0.0911) # finds the threshold channel error rate for 50% of probability of using A_1 and q=49%. It starts looking for the threshold close to delta=9.11%

The files contained in the folder Biases allow to compute the key rate in the case in which we use a single bases of Alice to generate the key. The key rate is computed as a function of the detection efficiency and of the visibility (that can be computed as a function of the channel error rate delta as v=1-2delta).

Again, one can use the program creating a Jupyter Notebook in the same folder where the Julia files are contained and typing

    include("entropies.jl")
    include("dibounds.jl")
    include("quantumcorrelations.jl")
    include("keyrate.jl");
    
    (eta,k,hae,hab,a1,s,x)=find_best_eta_n([0.4 0 pi/2 pi/8],1,0.0,1e-6,convexbound) # finds the threshold detection efficiency assuming v=1, q=0. The program will take the threshold as the value of eta at which the key rate becomes smaller than 10^(-6)
    
    certifykeyrate(x...,1,eta,0.,1e-10,0.) # with our method it is important to certify the key rate that we found running this command. If the output is "true", we can consider the values obtained proven
