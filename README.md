# SineFit.jl

![Example plot](example.svg)

This package provides functions for estimating and accurately fitting a sine
function to a given signal.

This is basically an improvement in robustness to only using `LsqFit`: clever
initial parameter guessing helps the least-squares fit converge quickly and
accurately.
