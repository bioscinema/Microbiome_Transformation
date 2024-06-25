import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import skew
import warnings
warnings.filterwarnings("ignore") 

def dual_group_boxcox_transformation2(p1, p2):
    
    # Combining both likelihood functions into one
    def combined_negative_log_likelihood(params):
        lambda_ = params[0]
        total_nll = 0
        
        for p_i in [p1, p2]:
            # Box-Cox transformation
            if lambda_ != 0:
                y_i = (p_i ** lambda_ - 1) / lambda_
            else:
                y_i = np.log(p_i)
                
            mu_hat = np.mean(y_i)
            sigma_hat_squared = np.sum((y_i - mu_hat) ** 2) / len(p_i)
            n = len(p_i)
            
            nll = -(n / 2) * np.log(2 * np.pi * sigma_hat_squared) - (n / 2) + (lambda_ - 1) * np.sum(np.log(p_i))
            total_nll += -nll 
        
        return total_nll
    
    bounds = [(-10, 10)]
    initial_guess = [0.1]  # A reasonable initial guess for lambda
    result = minimize(combined_negative_log_likelihood, initial_guess, method='L-BFGS-B', bounds=bounds)
    optimized_lambda = result.x[0]
    optimized_lambda = round(optimized_lambda, 2)
    
    # Transforming the samples using the optimized parameters
    def transform(p, lambda_):
        if lambda_ != 0:
            return (p ** lambda_ - 1) / lambda_
        else:
            return np.log(p)
    
    transformed_p1 = transform(p1, optimized_lambda)
    transformed_p2 = transform(p2, optimized_lambda)
    
    return transformed_p1, transformed_p2



def new_logit_transformation_s(p1, p2):
    
    def combined_negative_log_likelihood(params):
        phi, varphi = params
        total_nll = 0
        for p_i in [p1, p2]:
            p_i = np.array(p_i)
            log_args_1 = (p_i + phi) / (1 - p_i + varphi)
            log_args_2 = (phi + varphi + 1)
            log_args_3 = (p_i + phi) * (varphi + 1 - p_i)
            y_i = np.log(log_args_1)
            mu_hat = np.mean(y_i)
            sigma_hat_squared = np.sum((y_i - mu_hat)**2) / (len(p_i))
            n = len(p_i)
            nll = -(n/2)*np.log(2*np.pi*sigma_hat_squared) - 1/(2*sigma_hat_squared) * np.sum((y_i - mu_hat)**2) + n*np.log(log_args_2) - np.sum(np.log(log_args_3))
            total_nll += -nll 
        return total_nll
    bounds = [(0, None), (0, None)]
    initial_guess = [0.001, 0]
    result = minimize(combined_negative_log_likelihood, initial_guess, method='L-BFGS-B', bounds=bounds)
    optimized_phi, optimized_varphi = result.x
    
    def transform(p, phi, varphi):
        return np.log((p + phi) / (1 - p + varphi))
    p1 = np.array(p1)
    p2 = np.array(p2)
    transformed_p1 = transform(p1, optimized_phi, optimized_varphi)
    transformed_p2 = transform(p2, optimized_phi, optimized_varphi)
    
    return transformed_p1, transformed_p2