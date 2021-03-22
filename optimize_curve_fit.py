from scipy import optimize 
def fit_curve(x,y):
    var_fun = lambda mu, phi: mu + phi * mu ** 2
    phi_hat, _ = optimize.curve_fit(var_fun, x, y)
    return(phi_hat[0])
