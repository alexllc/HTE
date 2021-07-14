
library("survivalsvm")

survivalsvm(formula = NULL, 
            data = NULL, 
            subset = NULL,
            type = "regression", 
            diff.meth = NULL, 
            gamma.mu = NULL,
            opt.meth = "quadprog", 
            kernel = "lin_kernel", 
            kernel.pars = NULL,
            time.variable.name = NULL, 
            status.variable.name = NULL, 
            sgf.sv = 5,
            sigf = 7, 
            maxiter = 20, 
            margin = 0.05, 
            bound = 10, 
            eig.tol = 1e-06,
            conv.tol = 1e-07, 
            posd.tol = 1e-08
            )