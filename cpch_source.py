import scipy
from scipy.stats import norm, chi2, truncnorm, bernoulli, t, uniform
import numpy as np
from numpy.random import seed, normal
from timeit import default_timer as timer
import pandas as pd
from functools import partial
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from itertools import product, combinations, chain
import copy
import argparse
import statsmodels.stats.multitest as mt
import rpy2
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
#make plots not look bad
sns.set(font='DejaVu Sans',
        rc={
 'axes.axisbelow': False,
 'axes.edgecolor': 'lightgrey',
 'axes.facecolor': 'None',
 'axes.grid': False,
 'axes.labelcolor': 'dimgrey',
 'axes.spines.right': False,
 'axes.spines.top': False,
 'figure.facecolor': 'white',
 'lines.solid_capstyle': 'round',
 'patch.edgecolor': 'w',
 'patch.force_edgecolor': True,
 'text.color': 'dimgrey',
 'xtick.bottom': False,
 'xtick.color': 'dimgrey',
 'xtick.direction': 'out',
 'xtick.top': False,
 'ytick.color': 'dimgrey',
 'ytick.direction': 'out',
 'ytick.left': False,
 'ytick.right': False})
sns.despine(left = True, bottom = True)
sns.set_context("notebook", rc={"font.size":16,
                                "axes.titlesize":20,
                                "axes.labelsize":18})


lookup_table = pd.read_csv('lookup_table.csv')

###################
## Utility Fxns ###
###################

def generate_S(m, r):
    '''
    @ param m: int, m>= 2, number of base hypotheses
    @ param r: int, 2 <= r <= m

    Generates the set S as presented in the Computing cPCH p-values section of the paper
    '''
    #gets all combinations of length n-r+1 from range(1, 2, ..., n)
    S = list(itertools.combinations(list(np.arange(1, m+1, 1)), m-r+1))
    S = [list(ele) for ele in S]
    remaining = []
    #get indices that are not present in S
    for ls in S:
        remaining.append(list(list(set(np.arange(1, m+1, 1))-set(ls))))
    test = []
    for i in np.arange(len(S)):
        small = S[i]
        rest = remaining[i]
        all_perms = list(itertools.permutations(rest))
        final = list(itertools.product([small], all_perms))
        test.append(final)
    return(list(chain.from_iterable(test)))

####################################
# Generate Data for MT Experiments #
####################################

def hierarchical(nature, p, pi1):
    if bernoulli.rvs(p) == 1:
        return nature
    else:
        return bernoulli.rvs(pi1)

vhierarchical = np.vectorize(hierarchical)

def generate_indicator_matrix(pi1, p, mus, m, r, M):
    '''
    generates matrix, 0 for null, 1 for nonnull
    '''
    nature = bernoulli.rvs(pi1, size = M)
    null_nonnull = []
    for i in np.arange(1, m+1, 1):
         null_nonnull.append(vhierarchical(nature, p, pi1))
    null_nonnull = np.array(null_nonnull).transpose()
    return(null_nonnull)

def generate_XX_mt(pi1, p, mus, m, r, M, power_iters = 1):
    '''
    generates data based on null-nonnull indicator mu_matrix (0 for null, 1 for nonnull)
    '''
    null_nonnull = generate_indicator_matrix(pi1, p, mus, m, r, M)
    mu_matrix = np.random.choice(mus, size = (M, m)) * null_nonnull
    #sample T from normal with above true mu
    T = norm.rvs(mu_matrix, scale = 1, size = (power_iters, M, m))
    return T, null_nonnull


######################
# Truncated Sampling #
######################
def trunc_t(truncs, dof, loc, scale, N):
    '''
    @ param truncs: a list of 2 float values containing the lower and upper bound of the truncation
    @ param dof: int, the DOF to use for the t-distribution
    @ param loc: float, the location to use for the t-distribution
    @ param scale: float, the scale of the t-distribution
    @ param N: int, the number of samples from the specified t-distribution to take

    Outputs a N length vector of samples from a t-distribution with a particular location and dof (scale = 1)
    '''
    if min(truncs) == max(truncs):
        return np.repeat(min(truncs), N)
    F_a =  t.cdf(min(truncs), dof, loc, scale)
    F_b = t.cdf(max(truncs), dof, loc, scale)
    #uniform on [F_a, F_b]
    u = uniform.rvs(loc = F_a, scale = F_b - F_a, size = N)
    #ppf is quantile function
    return (t.ppf(u, dof, loc))


#######################
# Combining Functions #
#######################
def f_bonferroni(samples, axis):
    '''
    @ param samples: a numpy array of numpy arrays, each element representing
    the samples from a particular truncnorm representing the m-r+1 smallest test statistics, e.g. np.array([t1, t2])
    @ param axis:  int for computing the Bonferroni combining function down an axis of a numpy array
    '''
    return np.max(np.abs(samples), axis)

def f_fisher(pvals, axis):
    '''
    @ param pvals: a numpy array of np arrays containing the m-r+1 largest p-values (aka those associated with the m-r+1 smallest test statistics)
    @ param axis: int for computing the Simes combining function down an axis of a numpy array
    '''
    return (-2 * np.sum(np.log(pvals), axis))

def f_simes(pvals, axis):
    '''
    @ param pvals: a numpy array of np arrays containing the m-r+1 largest p-values (aka those associated with the m-r+1 smallest test statistics)
    @ param axis:  int for computing the Fisher combining function down an axis of a numpy array
    '''
    if len(pvals.shape) == 1:
        pvals = pvals.reshape(1, pvals.shape[0], 1)
    if len(pvals.shape) == 2:
    #need to shape into a 3-d array so that we can handle XXs and the null samples
        pvals = pvals.reshape(pvals.shape[0], pvals.shape[1], 1)
    mr1 = pvals.shape[axis]
    sorted_pvals = np.sort(pvals, axis)
    i = np.arange(1, mr1+1, 1)
    sorted_pvals
    const = np.tile(mr1/i.reshape(mr1, 1), pvals.shape[2])
    #calculates n-r+1/i * each column of sorted_pvals
    simes = np.multiply(const, sorted_pvals)
    #gets simes test statistic
    return np.min(simes, axis)


######################
# Storey's Procedure #
######################
#As detailed in https://candes.su.domains/teaching/stats300c/Lectures/Lecture08.pdf

def storey(p, fdr_level):
    '''
    @ param p: np.array of p-values on which to perform Storey's procedure
    @ fdr_level: float, nominal level to perform Storey's procedure

    Outputs p-length np.array of rejections (1 for rejection, 0 for no rejection)
    '''
    qvals = p
    rm_na = np.logical_not(np.isnan(p))
    p = p[rm_na]
    pi0hat = min(1, (1+len(p) - sum(p <= 0.5))/(len(p)/2))
    m = len(p)
    i = np.arange(m, 0, -1)
    o = np.argsort(-p)
    ro = np.argsort(o)
    for j in o:
        p[j] * m/i
    qvals = pi0hat * np.minimum(np.minimum.accumulate(p[o] * m /i), 1)[ro]
    return (qvals <= fdr_level)

###############################################
# Vectorized Fxns for Computing cPCH p-values #
###############################################

def get_thetahat_T_small(T, m, r):
    '''
    @param T: numpy array of base test statistics
    @param m: The number of base hypotheses
    @param r: The replicability constant
    @param df_mat: a matrix the same dimension as T containing the associated DOFs

    Outputs T_{(1:m-r+1)}
    '''
    if T.shape[1] != m:
        raise ValueError("Each row must have length m")
    Vthetahat = T.copy()
    repl = T.shape[0]
    ignore = (np.abs(T).argsort(axis=1).argsort(axis = 1) < m-r+1)
    Vthetahat[ignore] = 0
    Ts = T[ignore].reshape(repl, m-r+1)
    return Vthetahat, Ts

def get_df_t_dist(T, m, r, df_mat):
    '''
    @param T: numpy array of base test statistics (num columns = m)
    @param m: The number of base hypotheses
    @param r: The replicability constant
    @param df_mat: a matrix the same dimension as T containing the associated DOFs

    Outputs associated DOFs for T_{(1:m-r+1)} for t-distibution base test statistics
    '''
    if T.shape[1] != m:
        raise ValueError("Each row must have length m")
    repl = T.shape[0]
    ignore = (np.abs(T).argsort(axis=1).argsort(axis = 1) < m-r+1)
    df_Xs = df_mat[ignore].reshape(repl, m-r+1)
    return df_Xs

def get_T_big(T, m, r):
    '''
    @param T: numpy array of base test statistics (num columns = m)
    @param m: The number of base hypotheses
    @param r: The replicability constant

    Outputs the value of the observed T_{(m-r+2:m)}
    '''
    M = T.shape[0]
    idx = (np.abs(T).argsort(axis=1).argsort(axis = 1) >= m-r+1)
    T_big = T[idx].reshape(M, r-1)
    arr_inds = np.argsort(np.abs(T_big), axis = 1)
    final_t_big = np.take_along_axis(np.sign(T_big),arr_inds,axis=1)* np.sort(np.abs(T_big), axis = 1)
    return final_t_big

def get_oracle_theta(true_theta, T, m, r):
    '''
    @param true_theta: a M x m array of true means associated with each element of X, nonzero elements should always
                    appear first in descending order based on magnitude e.g. np.array([-3, 2, 1, 0, 0])
    @param T: numpy array of base test statistics, (num columns = m)
    @param m: The number of base hypotheses
    @param r: The replicability constant
    Outputs M x m matrix where every row is the oracle theta under the PC null, a vector
    identical to theta except with only the top r-1 entries in magnitude remaining the same and the rest set to 0
    '''
    if true_theta.shape[1] != m:
        raise ValueError("number of columns of true_theta must equal m")
    M = T.shape[0]
    #true_theta_sorted = np.sort(true_theta, axis = 1)
    zero_pads = np.repeat(np.repeat(0, M), m-r+1).reshape(M, m-r+1)
    tm = np.concatenate((true_theta[:, :r-1], zero_pads), axis = 1)
    return(tm)

def cpch_m2_r2(T, m, r):
    '''
    @param T: numpy array of base test statistics (num columns = 2)
    @param m: The number of base hypotheses
    @param r: The replicability constant

    Outputs cPCH p-values for the r=m=2 setting
    '''
    Vmin_val = np.abs(T).min(axis = 1)
    Vmax_val = np.abs(T).max(axis = 1)
    Vthetahat = get_thetahat_T_small(T, m, r)[0]

    idx = np.argmax(abs(T), axis = 1)
    M = T.shape[0]
    T_2 = T[np.arange(M), idx]

    mixture_comp_1 = (norm.cdf(Vmax_val,Vthetahat[:,0],1)-norm.cdf(Vmin_val,Vthetahat[:,0],1) +
                           norm.cdf(-Vmin_val,Vthetahat[:,0],1)-norm.cdf(-Vmax_val,Vthetahat[:,0],1))/(norm.cdf(Vmax_val,Vthetahat[:,0],1)-norm.cdf(-Vmax_val,Vthetahat[:,0],1))

    mix_weight_1 = (norm.pdf(T_2,Vthetahat[:,1], 1))*(norm.cdf(Vmax_val-Vthetahat[:,0]) - norm.cdf(-Vmax_val-Vthetahat[:,0]))

    mixture_comp_2 = (-norm.cdf(Vmin_val,Vthetahat[:,1],1)+norm.cdf(Vmax_val,Vthetahat[:,1],1) +
                           norm.cdf(-Vmin_val,Vthetahat[:,1],1)-norm.cdf(-Vmax_val,Vthetahat[:,1],1))/(norm.cdf(Vmax_val,Vthetahat[:,1],1)-norm.cdf(-Vmax_val,Vthetahat[:,1],1))
    mix_weight_2 = (norm.pdf(T_2,Vthetahat[:,0], 1))*(norm.cdf(Vmax_val-Vthetahat[:,1]) - norm.cdf(-Vmax_val-Vthetahat[:,1]))

    P = mix_weight_1 + mix_weight_2

    pvals = mixture_comp_1 * (mix_weight_1/P) + mixture_comp_2 * (mix_weight_2/P)

    return(pvals)

def cpch_m2_r2_orac(T, m, r, true_loc):
    '''
    @ param T: numpy array of base test statistics (num columns = 2)
    @param m: The number of base hypotheses
    @param r: The replicability constant

    Outputs cPCH oracle p-values for the r=m=2 setting
    '''
    Vmin_val = np.abs(T).min(axis = 1)
    Vmax_val = np.abs(T).max(axis = 1)
    Vthetahat = get_oracle_theta(true_loc, T, m, r)

    idx = np.argmax(abs(T), axis = 1)
    M = T.shape[0]
    T_2 = T[np.arange(M), idx]

    mixture_comp_1 = (norm.cdf(Vmax_val,Vthetahat[:,0],1)-norm.cdf(Vmin_val,Vthetahat[:,0],1) +
                           norm.cdf(-Vmin_val,Vthetahat[:,0],1)-norm.cdf(-Vmax_val,Vthetahat[:,0],1))/(norm.cdf(Vmax_val,Vthetahat[:,0],1)-norm.cdf(-Vmax_val,Vthetahat[:,0],1))

    mix_weight_1 = (norm.pdf(T_2,Vthetahat[:,1], 1))*(norm.cdf(Vmax_val-Vthetahat[:,0]) - norm.cdf(-Vmax_val-Vthetahat[:,0]))

    mixture_comp_2 = (-norm.cdf(Vmin_val,Vthetahat[:,1],1)+norm.cdf(Vmax_val,Vthetahat[:,1],1) +
                           norm.cdf(-Vmin_val,Vthetahat[:,1],1)-norm.cdf(-Vmax_val,Vthetahat[:,1],1))/(norm.cdf(Vmax_val,Vthetahat[:,1],1)-norm.cdf(-Vmax_val,Vthetahat[:,1],1))
    mix_weight_2 = (norm.pdf(T_2,Vthetahat[:,0], 1))*(norm.cdf(Vmax_val-Vthetahat[:,1]) - norm.cdf(-Vmax_val-Vthetahat[:,1]))

    P = mix_weight_1 + mix_weight_2

    pvals = mixture_comp_1 * (mix_weight_1/P) + mixture_comp_2 * (mix_weight_2/P)

    return(pvals)


def cpch_mc_norm(Ts, S, f, null_theta, T_big, N):
    '''
    Performs the MC sampling scheme outlined in the Computation section and outputs the cpch p-value
    Ts: A repl x 2 matrix of the m-r+1 smallest test statistics we observed
    S: A list of lists containing all the ways to some combination T1, ..., Tm corresponding to T_(1:m-r+1), and some ordered sets
    corresponding to the remaining order stats
    f: the function to use on the samples, takes in a N x n-r+1 dimensional matrix (N = 1, 2, ...)
    null_theta: repl x n matrix of theta_hats for each observation
    T_big: repl x r-1 vector of the observed T_(m-r+2:m), sets truncation for truncated normal sampling
    N: The number of Monte Carlo simulations to run for the sampling
    '''
    Ts_pvals = 2 * (1-norm.cdf(np.abs(Ts), loc = 0, scale = 1))
    if f == f_simes or f == f_fisher:
        test_stat = f(np.array(Ts_pvals), axis = 1)
    elif f == f_bonferroni:
        test_stat = f(np.array(Ts), axis = 1)
    else:
        test_stat = f(np.array(Ts), axis = 1)
    repl = Ts.shape[0]
    #cond_thres is repl x 1 vector of the observed T_(m-r+2), sets truncation for truncated normal sampling
    cond_thres = T_big[:,0]
    joint_probs, null_ts_probs = [], []
    for i in np.arange(len(S)):
        #each S[i][0] represents the set of T_i corresponding to the smallest order statistics
        small = S[i][0]
        #indices: 1D array of indices corresponding to each T_i that correspond to the smallest order statistics
        indices = np.array(small) - 1
        #so we take these indices, and subset each row of theta_hat with this indices vector
        #this gives us a repl x m-r+1 matrix of the theta for the T_i corresponding to the smallest test statistics
        small_theta = null_theta[:,indices]
        #get T_(m-r+2)
        theta_thres = null_theta[:, list(S[i][1])[0]-1]


        if len(S[i][1]) == 1:
            large_theta = []
        else:
        #get list of mu corresponding to T_i larger than T_{(m-r+2)}, as long as it is not empty
            l_indices = np.array(list(S[i][1])[1:]) - 1
            #so we take these indices, and subset each row of mu_hat with this indices vector
            #this gives us a repl x r-2 matrix of the predicted mu for the T_i corresponding to the largest test statistics
            large_theta = null_theta[:,l_indices]

        # # calculate P(B_i|T_{(m-r+2:m)}) terms
        joint_prob_i = norm.pdf(cond_thres, loc = theta_thres, scale = 1)
        s_probs = norm.cdf(np.abs(cond_thres).reshape(repl, 1), loc = small_theta, scale = 1) - norm.cdf(-np.abs(cond_thres).reshape(repl, 1), loc = small_theta, scale = 1)
        #multiply together the columns so you're multiplying across the small_theta
        joint_prob_i *= np.prod(s_probs, axis = 1)
        if len(S[i][1]) > 1:
            l_probs = norm.pdf(T_big[:,1:], loc = large_theta, scale = 1)
            joint_prob_i *= np.prod(l_probs, axis = 1)
        #joint probs will be a len(S) x repl matrix
        joint_probs.append(joint_prob_i)

        #calculate P(f(T_(1:m-r+1)) > f^obs|T_{(m-r+2:m)}) terms
        # #sample from truncnorm with corresponding null_theta
        null_samples = []
        #loop over replications
        for i in range(small_theta.shape[0]):
            samples = []
            #look over the values in small_theta
            for j in range(small_theta.shape[1]):
                a = -np.abs(cond_thres[i]) - small_theta[i,j]
                b =  np.abs(cond_thres[i]) - small_theta[i,j]
                if a == b:
                    samp = np.repeat(a, N)
                else:
                    samp = truncnorm.rvs(a = -np.abs(cond_thres[i]) - small_theta[i,j], b = np.abs(cond_thres[i]) - small_theta[i,j], loc = small_theta[i,j], scale = 1, size = N)
                samples.append(samp)
            null_samples.append(samples)
        samples_pvals = 2 * (1-norm.cdf(np.abs(np.array(null_samples)), loc = 0, scale = 1))
        if f == f_bonferroni:
            #axis = 1 means apply this down the rows of each matrix in an array of matrices, which is what we want since each row of each matrix in null_samples represents a single sample
            func_evals = f(np.array(null_samples), axis = 1)
            null_ts_probs.append((1 + np.sum(func_evals >= test_stat.reshape(repl, 1), axis = 1))/(N+1))
        elif f == f_fisher:
            #need to transform null_samples into pvalues to plug into fisher and simes function
            func_evals = f(np.array(samples_pvals), axis = 1)
            null_ts_probs.append((1 + np.sum(func_evals >= test_stat.reshape(repl, 1), axis = 1))/(N+1))
        elif f == f_simes:
            #need to transform null_samples into pvalues to plug into fisher and simes function
            func_evals = f(np.array(samples_pvals), axis = 1)
            null_ts_probs.append((1 + np.sum(func_evals <= test_stat.reshape(repl, 1), axis = 1))/(N+1))
        else:
            func_evals = f(np.array(null_samples), axis = 1)
            null_ts_probs.append((1 + np.sum(func_evals >= test_stat.reshape(repl, 1), axis = 1))/(N+1))
    cond_probs = np.array(joint_probs)/np.sum(np.array(joint_probs), axis = 0)
    #perform element-wise thetaltiplication and sum the columns, gives us a 1 x repl matrix of cpch-pvals as desired
    cpch_pvals = np.array(null_ts_probs) * cond_probs
    return np.sum(cpch_pvals, axis = 0)

def cpch_mc_t(Ts, S, f, df, null_theta, T_big, N, df_Ts):
    '''
    Performs the MC sampling scheme outlined in the Computation section and outputs the cpch p-value
    Ts: A repl x 2 matrix of the m-r+1 smallest test statistics we observed
    S: A list of lists containing all the ways to some combination T1, ..., Tm corresponding to T_(1:m-r+1), and some ordered sets
    corresponding to the remaining order stats
    f: the function to use on the samples, takes in a N x n-r+1 dimensional matrix (N = 1, 2, ...)
    df: a repl x n matrix of degrees of freedom for each observation
    null_theta: repl x n matrix of theta_hats for each observation
    T_big: repl x r-1 vector of the observed T_(m-r+2:m), sets truncation for truncated normal sampling
    N: The number of Monte Carlo simulations to run for the sampling
    df_Ts: degrees of freedom corresponding to Ts
    '''
    Ts_pvals = 2 * (1-t.cdf(np.abs(Ts), df_Ts, loc = 0, scale = 1))
    if f == f_simes or f == f_fisher:
        test_stat = f(np.array(Ts_pvals), axis = 1)
    elif f == f_bonferroni:
        test_stat = f(np.array(Ts), axis = 1)
    else:
        test_stat = f(np.array(Ts), axis = 1)
    repl = Ts.shape[0]
    cond_thres = T_big[:,0]
    joint_probs, null_ts_probs = [], []
    for i in np.arange(len(S)):
        #each S[i][0] represents the set of T_i corresponding to the smallest order statistics
        small = S[i][0]
        #indices: 1D array of indices corresponding to each T_i that correspond to the smallest order statistics
        indices = np.array(small) - 1
        #so we take these indices, and subset each row of theta_hat with this indices vector
        #this gives us a repl x m-r+1 matrix of the predicted theta for the T_i corresponding to the smallest test statistics
        small_theta = null_theta[:,indices]
        small_df = df[:, indices]
        #get the m-r+2 order stat
        theta_thres = null_theta[:, list(S[i][1])[0]-1]
        df_thres = df[:,  list(S[i][1])[0]-1]
        #gives a repl x 1 array (picks out the S[i][1]-1th column of null_theta)
        if len(S[i][1]) == 1:
            large_theta = []
            large_df = []
        else:
        #get list of mu corresponding to T_i larger than T_{(m-r+2)}, as long as it is not empty
            l_indices = np.array(list(S[i][1])[1:]) - 1
            #so we take these indices, and subset each row of mu_hat with this indices vector
            #this gives us a repl x r-2 matrix of the predicted mu for the T_i corresponding to the largest test statistics
            large_theta = null_theta[:,l_indices]
            large_df = df[:,l_indices]
        # # calculate P(B_i|T_{(m-r+2:m)}) terms
        joint_prob_i = norm.pdf(cond_thres, loc = theta_thres, scale = 1)
        s_probs = norm.cdf(np.abs(cond_thres).reshape(repl, 1), loc = small_theta, scale = 1) - norm.cdf(-np.abs(cond_thres).reshape(repl, 1), loc = small_theta, scale = 1)
        #multiply together the columns so you're multiplying across the small_theta
        joint_prob_i *= np.prod(s_probs, axis = 1)
        if len(S[i][1]) > 1:
            l_probs = norm.pdf(T_big[:,1:], loc = large_theta, scale = 1)
            joint_prob_i *= np.prod(l_probs, axis = 1)
        #joint probs will be a len(S) x repl matrix
        joint_probs.append(joint_prob_i)

        #calculate P(f(T_(1:m-r+1)) > f^obs|T_{(m-r+2:m)}) terms
        # #sample from truncnorm with corresponding null_theta
        null_samples, null_pvals = [], []

        #loop over replications
        for i in range(small_theta.shape[0]):
            samples, sample_pvals = [], []
            #look over the values in small_theta
            for j in range(small_theta.shape[1]):
                a = -np.abs(cond_thres[i])
                b =  np.abs(cond_thres[i])
                samp = trunc_t([a,b], dof = small_df[i, j], loc = small_theta[i, j], scale = 1, N = N)
                samp_pval = 2 * (1-t.cdf(np.abs(samp), small_df[i, j], loc = 0, scale = 1))
                samples.append(samp)
                sample_pvals.append(samp_pval)
            null_samples.append(samples)
            null_pvals.append(sample_pvals)
        if f == f_bonferroni:
            func_evals = f(np.array(null_samples), axis = 1)
            null_ts_probs.append((1 + np.sum(func_evals >= test_stat.reshape(repl, 1), axis = 1))/(N+1))
        elif f == f_fisher:
            #need to transform null_samples into pvalues to plug into fisher and simes function
            func_evals = f(np.array(null_pvals), axis = 1)
            null_ts_probs.append((1 + np.sum(func_evals >= test_stat.reshape(repl, 1), axis = 1))/(N+1))
        elif f == f_simes:
            #need to transform null_samples into pvalues to plug into fisher and simes function
            func_evals = f(np.array(null_pvals), axis = 1)
            null_ts_probs.append((1 + np.sum(func_evals <= test_stat.reshape(repl, 1), axis = 1))/(N+1))
        else:
            func_evals = f(np.array(null_samples), axis = 1)
            null_ts_probs.append((1 + np.sum(func_evals >= test_stat.reshape(repl, 1), axis = 1))/(N+1))
    cond_probs = np.array(joint_probs)/np.sum(np.array(joint_probs), axis = 0)
    #perform element-wise multiplication and sum the columns, gives us a 1 x repl matrix of cpch pvals as desired
    cpch_pvals = np.array(null_ts_probs) * cond_probs
    return np.sum(cpch_pvals, axis = 0)

#sampling for general m and r
def cpch_mc(Ts, S, f, pdf, cdf, trunc_rvs, null_theta, T_big, N):
    '''
    Performs the MC sampling scheme outlined in the Computation section and outputs the cpch p-value
    Ts: A repl x 2 matrix of the m-r+1 smallest test statistics we observed
    S: A list of lists containing all the ways to some combination T1, ..., Tm corresponding to T_(1:m-r+1), and some ordered sets
    corresponding to the remaining order stats
    N: The number of Monte Carlo simulations to run for the sampling
    f: the function to use on the samples, takes in a N x n-r+1 dimensional matrix (N = 1, 2, ...)
    null_theta: repl x n matrix of theta_hats for each observation
    '''
    lower_tail = cdf(Ts, loc = 0)
    Ts_pvals = 2*np.abs((lower_tail >= 0.5)-lower_tail)
    if f == f_simes or f == f_fisher:
        test_stat = f(np.array(Ts_pvals), axis = 1)
    elif f == f_bonferroni:
        test_stat = f(np.array(Ts), axis = 1)
    else:
        test_stat = f(np.array(Ts), axis = 1)
    #generate the test stat to compare our samples to, just apply f along the rows of our given Ts matrix
    repl = Ts.shape[0]
    cond_thres = T_big[:,0]
    joint_probs, null_ts_probs = [], []
    for i in np.arange(len(S)):
        small = S[i][0]
        indices = np.array(small) - 1
        small_theta = null_theta[:,indices]
        theta_thres = null_theta[:, list(S[i][1])[0]-1]
        if len(S[i][1]) == 1:
            large_theta = []
        else:
            l_indices = np.array(list(S[i][1])[1:]) - 1
            large_theta = null_theta[:,l_indices]

        # # calculate P(B_i|T_{(m-r+2:m)}) terms
        joint_prob_i = norm.pdf(cond_thres, loc = theta_thres, scale = 1)
        s_probs = norm.cdf(np.abs(cond_thres).reshape(repl, 1), loc = small_theta, scale = 1) - norm.cdf(-np.abs(cond_thres).reshape(repl, 1), loc = small_theta, scale = 1)
        #multiply together the columns so you're multiplying across the small_theta
        joint_prob_i *= np.prod(s_probs, axis = 1)
        if len(S[i][1]) > 1:
            l_probs = norm.pdf(T_big[:,1:], loc = large_theta, scale = 1)
            joint_prob_i *= np.prod(l_probs, axis = 1)
        #joint probs will be a len(S) x repl matrix
        joint_probs.append(joint_prob_i)

        #calculate P(f(T_(1:m-r+1)) > f^obs|T_{(m-r+2:m)}) terms
        # #sample from truncnorm with corresponding null_theta
        null_samples, null_pvals = [], []

        #loop over replications
        for i in range(small_theta.shape[0]):
            samples, sample_pvals = [], []
            #look over the values in small_theta
            for j in range(small_theta.shape[1]):
                a = -np.abs(cond_thres[i])
                b =  np.abs(cond_thres[i])
                samp = trunc_rvs(truncs = [a,b], loc = small_theta[i,j], N = N)
                samp_pval = 2 * (1 - cdf(np.abs(samp),  loc = 0))
                samples.append(samp)
                sample_pvals.append(samp_pval)
            null_samples.append(samples)
            null_pvals.append(sample_pvals)
        if f == f_bonferroni:
            func_evals = f(np.array(null_samples), axis = 1)
            null_ts_probs.append((1 + np.sum(func_evals >= test_stat.reshape(repl, 1), axis = 1))/(N+1))
        elif f == f_fisher:
            #need to transform null_samples into pvalues to plug into fisher and simes function
            func_evals = f(np.array(null_pvals), axis = 1)
            null_ts_probs.append((1 + np.sum(func_evals >= test_stat.reshape(repl, 1), axis = 1))/(N+1))
        elif f == f_simes:
            #need to transform null_samples into pvalues to plug into fisher and simes function
            func_evals = f(np.array(null_pvals), axis = 1)
            null_ts_probs.append((1 + np.sum(func_evals <= test_stat.reshape(repl, 1), axis = 1))/(N+1))
        else:
            func_evals = f(np.array(null_samples), axis = 1)
            null_ts_probs.append((1 + np.sum(func_evals <= test_stat.reshape(repl, 1), axis = 1))/(N+1))
    cond_probs = np.array(joint_probs)/np.sum(np.array(joint_probs), axis = 0)
    #perform element-wise multiplication and sum the columns, gives us a 1 x repl matrix of cpch pvals as desired
    cpch_pvals = np.array(null_ts_probs) * cond_probs
    return np.sum(cpch_pvals, axis = 0)


def cpch_unadjusted(T, m, r, f, pdf = norm.pdf, cdf = None, trunc_rvs = None, N = 10000, dof = []):
    '''
    @ param T: a matrix of base test statistics
    @ param m: int, m>= 2, number of base hypotheses
    @ param r: int, 2 <= r <= m
    @ param f: a combining function, see f_fisher as an example
    @ param pdf: a function for computing the pdf of a distribution (e.g. norm.pdf)
    @ param cdf: a function for computing the cdf of a distribution (e.g. norm.cdf)
    @ param trunc_rvs: a function for computing random draws from a truncated distribution (e.g. truncnorm.rvs)
    @ param N: The number of samples to use for the Monte Carlo procedure for computing cPCH p-values

    Outputs a cPCH p-value
    '''
    if pdf not in [norm.pdf, t.pdf] and (trunc_rvs == None or cdf == None):
        raise ValueError("pdf not norm.pdf or t.pdf but no cdf and/or truncated sampling function provided")
    if m == 2 and r == 2 and pdf == norm.pdf:
        return (cpch_m2_r2(T, m, r))
    S = generate_S(m, r)
    T_big = get_T_big(T, m, r)
    thetahat, Ts = get_thetahat_T_small(T, m, r)
    if pdf == norm.pdf:
        cpch_pvals = cpch_mc_norm(Ts, S, f, thetahat, T_big, N)
    elif pdf == t.pdf:
        df_Ts = get_df_t_dist(T, m, r, dof)
        cpch_pvals = cpch_mc_t(Ts, S, f, dof, thetahat, T_big, N, df_Ts)
    else:
        cpch_pvals = cpch_mc(Ts, S, f, pdf, cdf, trunc_rvs, thetahat, cond_thres, N)
    return (cpch_pvals)


def cpch_oracle_unadjusted(T, m, r, f, true_loc, pdf = norm.pdf, cdf = None, trunc_rvs = None, N =10000, dof = []):
    '''
    @ param T: a matrix of base test statistics
    @ param m: int, m>= 2, number of base hypotheses
    @ param r: int, 2 <= r <= m
    @ param f: a combining function, see f_fisher as an example
    @ param true_loc: a M x m array of true means associated with each element of T, nonzero elements should always
                    appear first in descending order based on magnitude e.g. np.array([-3, 2, 1, 0, 0])
    @ param pdf: a function for computing the pdf of a distribution (e.g. norm.pdf)
    @ param cdf: a function for computing the cdf of a distribution (e.g. norm.cdf)
    @ param trunc_rvs: a function for computing random draws from a truncated distribution (e.g. truncnorm.rvs)
    @ param N: The number of samples to use for the Monte Carlo procedure for computing cPCH p-values

    Outputs an unadjusted cPCH Oracle p-value
    '''
    if pdf not in [norm.pdf, t.pdf] and (trunc_rvs == None or cdf == None):
        raise ValueError("pdf not norm.pdf or t.pdf but no cdf and/or truncated sampling function provided")
    if m == 2 and r == 2 and pdf == norm.pdf:
        return (cpch_m2_r2_orac(T, m, r, true_loc))
    S = generate_S(m, r)
    repl = T.shape[0]
    T_big = get_T_big(T, m, r)
    thetahat, Ts = get_thetahat_T_small(T, m, r)
    oracle_theta = get_oracle_theta(true_loc, T, m, r)
    if pdf == norm.pdf:
        cpch_pvals =  cpch_mc_norm(Ts, S, f, oracle_theta, T_big, N)
    elif pdf == t.pdf:
        df_Ts = get_df_t_dist(T, m, r, dof)
        cpch_pvals = cpch_mc_t(Ts, S, f, dof, oracle_theta, T_big, N, df_Ts)
    else:
        cpch_pvals = cpch_mc(Ts, S, f, pdf, cdf, trunc_rvs, oracle_theta, cond_thres, N)
    return (cpch_pvals)


def cpch(T, m, r, f, pdf = norm.pdf, cdf = None, trunc_rvs = None, N = 10000, dof = []):
    #subset lookup table to right m, r, and method
    if f == f_fisher:
        method = 'Fisher'
    elif f==f_simes:
        method='Simes'
    else:
        raise Error('undefined combining function')
    unadjusted_cpch_pvals = cpch_unadjusted(T, m, r, f, pdf, cdf, trunc_rvs, N, dof)
    lookup_subset = lookup_table[(lookup_table['methods'] == method)*(lookup_table['m'] == m)*(lookup_table['r'] == r)]
    #get alpha-prime values from this subsetted table and calcuate alpha' closest to each cpch p-value
    alpha_primes = np.append(0, np.unique(lookup_subset[['alpha_prime']]))
    def get_closest_alpha_prime(pval):
        return(min(alpha_primes, key=lambda x:abs(x-pval)))
    
    v_get_closest_alpha_prime = np.vectorize(get_closest_alpha_prime)
    
    closest_alpha_prime = v_get_closest_alpha_prime(unadjusted_cpch_pvals)

    def get_next_up_index(alpha):
        if alpha == max(alpha_primes):
            return(list(alpha_primes).index(alpha))
        else:
            return(list(alpha_primes).index(alpha)+1)
    
    def calc_adjusted_pval(i):   
        if closest_alpha_prime[i] == max(alpha_primes):
            y2 = 1
            x2 = 1
            y1 = lookup_subset[lookup_subset['alpha_prime'] == alpha_primes[len(alpha_primes)-2]]['alpha'].values[0]
            x1 = alpha_primes[len(alpha_primes)-2]
            adjusted_pval = (y2-y1)/(x2-x1) * (unadjusted_cpch_pvals[i] - x1) + y1 
        elif closest_alpha_prime[i] == 0:
            #interpolate between (0, 0) and first point 
            y2 = lookup_subset[lookup_subset['alpha_prime'] == min(np.unique(lookup_subset[['alpha_prime']]))]['alpha'].values[0]
            y1 = 0
            x1 = 0
            x2 = min(np.unique(lookup_subset[['alpha_prime']]))
            adjusted_pval = max(0, (y2-y1)/(x2-x1) * (unadjusted_cpch_pvals[i] - x1) + y1)
        else:
            #interpolate between closest alpha' and next up
            y1 = lookup_subset[lookup_subset['alpha_prime'] == closest_alpha_prime[i]]['alpha'].values[0]
            y2 = lookup_subset[lookup_subset['alpha_prime'] == next_up_alpha_primes[i]]['alpha'].values[0]
            x1 = closest_alpha_prime[i]
            x2 = next_up_alpha_primes[i]
            adjusted_pval = (y2-y1)/(x2-x1) * (unadjusted_cpch_pvals[i] - x1) + y1
        return(adjusted_pval)

    v_calc_adjusted_pval = np.vectorize(calc_adjusted_pval)
    
    next_up_index = np.vectorize(get_next_up_index)

    next_up_alpha_primes = alpha_primes[next_up_index(closest_alpha_prime)]

    return(v_calc_adjusted_pval(np.arange(0, len(unadjusted_cpch_pvals))))


def cpch_oracle(T, m, r, f, true_loc, pdf = norm.pdf, cdf = None, trunc_rvs = None, N =10000, dof = []):
    #subset lookup table to right m, r, and method
    if f == f_fisher:
        method = 'Fisher'
    elif f==f_simes:
        method='Simes'
    else:
        raise ValueError('undefined combining function')
    unadjusted_oracle_cpch_pvals = cpch_oracle_unadjusted(T, m, r, f, true_loc, pdf, cdf, trunc_rvs, N, dof)
    lookup_subset = lookup_table[(lookup_table['methods'] == method)*(lookup_table['m'] == m)*(lookup_table['r'] == r)]
    #get alpha-prime values from this subsetted table and calcuate alpha' closest to each cpch p-value
    alpha_primes = np.append(0, np.unique(lookup_subset[['alpha_prime']]))
    def get_closest_alpha_prime(pval):
        return(min(alpha_primes, key=lambda x:abs(x-pval)))
    
    v_get_closest_alpha_prime = np.vectorize(get_closest_alpha_prime)
    
    closest_alpha_prime = v_get_closest_alpha_prime(unadjusted_oracle_cpch_pvals)

    def get_next_up_index(alpha):
        if alpha == max(alpha_primes):
            return(list(alpha_primes).index(alpha))
        else:
            return(list(alpha_primes).index(alpha)+1)
    
    def calc_adjusted_pval(i):   
        if closest_alpha_prime[i] == max(alpha_primes):
            y2 = 1
            x2 = 1
            y1 = lookup_subset[lookup_subset['alpha_prime'] == alpha_primes[len(alpha_primes)-2]]['alpha'].values[0]
            x1 = alpha_primes[len(alpha_primes)-2]
            adjusted_pval = (y2-y1)/(x2-x1) * (unadjusted_oracle_cpch_pvals[i] - x1) + y1 
        elif closest_alpha_prime[i] == 0:
            #interpolate between (0, 0) and first point 
            y2 = lookup_subset[lookup_subset['alpha_prime'] == min(np.unique(lookup_subset[['alpha_prime']]))]['alpha'].values[0]
            y1 = 0
            x1 = 0
            x2 = min(np.unique(lookup_subset[['alpha_prime']]))
            adjusted_pval = max(0, (y2-y1)/(x2-x1) * (unadjusted_oracle_cpch_pvals[i] - x1) + y1)
        else:
            #interpolate between closest alpha' and next up
            y1 = lookup_subset[lookup_subset['alpha_prime'] == closest_alpha_prime[i]]['alpha'].values[0]
            y2 = lookup_subset[lookup_subset['alpha_prime'] == next_up_alpha_primes[i]]['alpha'].values[0]
            x1 = closest_alpha_prime[i]
            x2 = next_up_alpha_primes[i]
            adjusted_pval = (y2-y1)/(x2-x1) * (unadjusted_oracle_cpch_pvals[i] - x1) + y1
        return(adjusted_pval)

    v_calc_adjusted_pval = np.vectorize(calc_adjusted_pval)
    
    next_up_index = np.vectorize(get_next_up_index)

    next_up_alpha_primes = alpha_primes[next_up_index(closest_alpha_prime)]

    return(v_calc_adjusted_pval(np.arange(0, len(unadjusted_oracle_cpch_pvals))))
