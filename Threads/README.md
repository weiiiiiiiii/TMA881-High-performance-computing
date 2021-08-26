We will use Newton’s method to practice programming with the C11 thread library. Give a function f, say f(x) = x^3 - 1, where x^3 denotes the third power of x, Newton’s method iteratively strives to find a root of f. Starting with a value x_0, one sets x_{k+1} = x_k - f(x_k) / f’(x_k). Convergence of this is unclear, and even if it does converge, it is unclear what root of f is approached.

Here is an example run of Newton’s method for three different (complex) values of x_0.
