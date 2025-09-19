# g and h functions which create the closed form solution for the EI ODE
# E(t) = g1(gamma, nu, alpha, t)e0 + g2(gamma, nu, alpha, t)i0
function power(a,b)
    return a^b
end
function g1(params)
    gamma, nu, alpha, t = params
    sqrt_term = sqrt(4 * alpha * gamma + (gamma - nu)^2)
    exp_term_t = exp(sqrt_term * t)  # exp( sqrt(4*alpha*gamma + (gamma - nu)^2) * t )  
    return (4 * alpha * (1 + exp_term_t) * gamma + 
              (gamma - nu) * ((1 + exp_term_t) * gamma + sqrt_term - nu - 
              exp_term_t * (sqrt_term + nu))) / 
              (2 * exp((gamma + sqrt_term + nu) * t / 2) * (4 * alpha * gamma + (gamma - nu)^2))  
end
function g2(params)
 gamma, nu, alpha, t = params
 sqrt_term = sqrt(4 * alpha * gamma + (gamma - nu)^2)
 exp_term_t = exp(sqrt_term * t)  # exp( sqrt(4*alpha*gamma + (gamma - nu)^2) * t )
 return (alpha * (-1 + exp_term_t)) / 
         (exp(((gamma + sqrt_term + nu) * t) / 2) * sqrt_term)
end 
# I(t) = h1(gamma, nu, alpha, t)e0 + h2(gamma, nu, alpha, t)i0
function h1(params)
    gamma, nu, alpha, t = params
    # Compute the common terms
sqrt_term = sqrt(4 * alpha * gamma + (gamma - nu)^2)
exp_term_t = exp(sqrt_term * t)  # This calculates exp(sqrt(4 * alpha * gamma + (gamma - nu)^2) * t)

# Now rewrite the expression in Julia
return((-1 + exp_term_t) * gamma) / 
         (exp(((gamma + sqrt_term + nu) * t) / 2) * sqrt_term)
end 
function h2(params)
    gamma, nu, alpha, t = params
    # Compute the common terms
sqrt_term = sqrt(4 * alpha * gamma + (gamma - nu)^2)
exp_term_t = exp(sqrt_term * t)  # This calculates exp(sqrt(4 * alpha * gamma + (gamma - nu)^2) * t)
return (4 * alpha * (1 + exp_term_t) * gamma + 
          (gamma - nu) * (gamma - sqrt_term + 
          exp_term_t * sqrt_term + 
          exp_term_t * (gamma - nu) - nu)) / 
         (2 * exp(((gamma + sqrt_term + nu) * t) / 2) * 
         (4 * alpha * gamma + (gamma - nu)^2))
end
