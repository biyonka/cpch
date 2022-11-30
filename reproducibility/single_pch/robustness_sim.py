#replace path with location of cpch_source.py locally
exec(open("../cpch_source.py").read())
###################
#### Set Params ###
###################
parser = argparse.ArgumentParser(description='Generate grid plot given some n')
parser.add_argument('m', metavar='m', type=int,
                    help='the number of studies for the PC test')
parser.add_argument('r', metavar='r', type=int,
                            help='replicability number')
parser.add_argument('k', metavar='k', type=int,
                    help='number of non nulls')
parser.add_argument('sig_strength', metavar='sig_strength', type=int,
                    help='the value of mu for non nulls')
parser.add_argument('alpha', metavar='alpha', type=float,
                    help='the alpha level of the PC test')
parser.add_argument('replicates', metavar='replicates', type=int,
                    help='the number of iterations to run for the t1 error and power estimation')
parser.add_argument('N', metavar='N', type=int,
                    help='the number of iterations to use for the MC cpch sampling')
a = parser.parse_args()
print(a)

m = a.m
r = a.r
k = a.k
sig_strength = a.sig_strength
alpha = a.alpha
replicates = a.replicates
N = a.N

dof_vec = np.array([1, 3, 5, 7, 9])

#set a unique seed per configuration, making it so each r is evaluated on
#same dataset with config n, k and ss
np.random.seed(int(20000*m + 1000*k + 10*sig_strength))

def get_pvals_t_robustness(m, r, k, sig_strength, replicates, df, alpha, N = 10000):
    '''
    df: a n by replicates numpy array (element-wise matching up w true_mu)
    '''
    true_mu = np.tile(np.append(np.repeat(sig_strength, k), np.repeat(0, m-k)), replicates).reshape(replicates, m)
    #generate the t-distributed data according to true mu
    XX = t.rvs(df, loc = true_mu, scale = 1,  size = (replicates, true_mu.shape[1]))
    pvals_cpch_f_t = cpch(XX, m, r, f_fisher, t.pdf, t.cdf, trunc_t, N, df)
    pvals_cpch_b_t = cpch(XX, m, r, f_bonferroni, t.pdf, t.cdf, trunc_t, N, df)
    pvals_cpch_s_t = cpch(XX, m, r, f_simes, t.pdf, t.cdf, trunc_t, N, df)
    #transform XX into p-values, then into normal test statistics
    t_pvals =  2 * (1-t.cdf(np.abs(XX), df, loc = 0, scale = 1))
    norm_ts = -norm.ppf(t_pvals/2)*np.sign(XX)
    pvals_cpch_f_norm = cpch(norm_ts, m, r, f_fisher, norm.pdf, norm.cdf, truncnorm.rvs, N)
    pvals_cpch_b_norm = cpch(norm_ts, m, r, f_bonferroni, norm.pdf, norm.cdf, truncnorm.rvs, N)
    pvals_cpch_s_norm = cpch(norm_ts, m, r, f_simes, norm.pdf, norm.cdf, truncnorm.rvs, N)
    return (pvals_cpch_f_t, pvals_cpch_s_t, pvals_cpch_b_t,
             pvals_cpch_f_norm, pvals_cpch_s_norm, pvals_cpch_b_norm)


def estimate_power_robustness_t(m, r, k, sig_strength, dof_vec, alpha, N, replicates):
    power_cpch_f_t, power_cpch_s_t, power_cpch_b_t  = [],[],[]
    power_cpch_f_norm, power_cpch_s_norm, power_cpch_b_norm  = [],[],[]
    se_cpch_f_t, se_cpch_s_t, se_cpch_b_t  = [],[],[]
    se_cpch_f_norm, se_cpch_s_norm, se_cpch_b_norm  = [],[],[]
    for elem in dof_vec:
        df = np.tile(np.repeat(elem, m), replicates).reshape(replicates, m)
        pvals_cpch_f_t, pvals_cpch_s_t, pvals_cpch_b_t, pvals_cpch_f_norm, pvals_cpch_s_norm, pvals_cpch_b_norm = get_pvals_t_robustness(m, r, k, sig_strength, replicates, df, alpha, N)
        power_cpch_f_t.append(np.mean(pvals_cpch_f_t <= alpha))
        power_cpch_s_t.append(np.mean(pvals_cpch_s_t <= alpha))
        power_cpch_b_t.append(np.mean( pvals_cpch_b_t <= alpha))
        power_cpch_f_norm.append(np.mean(pvals_cpch_f_norm <= alpha))
        power_cpch_s_norm.append(np.mean(pvals_cpch_s_norm <= alpha))
        power_cpch_b_norm.append(np.mean( pvals_cpch_b_norm <= alpha))
        se_cpch_f_t.append(np.std((pvals_cpch_f_t <= alpha), ddof = 1)/(replicates**0.5))
        se_cpch_s_t.append(np.std((pvals_cpch_s_t <= alpha), ddof = 1)/(replicates**0.5))
        se_cpch_b_t.append(np.std((pvals_cpch_b_t <= alpha), ddof = 1)/(replicates**0.5))
        se_cpch_f_norm.append(np.std((pvals_cpch_f_norm <= alpha), ddof = 1)/(replicates**0.5))
        se_cpch_s_norm.append(np.std((pvals_cpch_s_norm <= alpha), ddof = 1)/(replicates**0.5))
        se_cpch_b_norm.append(np.std((pvals_cpch_b_norm <= alpha), ddof = 1)/(replicates**0.5))
    return ([power_cpch_f_t, power_cpch_s_t, power_cpch_b_t, power_cpch_f_norm, power_cpch_s_norm, power_cpch_b_norm],
           [se_cpch_f_t, se_cpch_s_t, se_cpch_b_t, se_cpch_f_norm, se_cpch_s_norm, se_cpch_b_norm])

power, se = estimate_power_robustness_t(m, r, k, sig_strength, dof_vec, alpha, N, replicates)
methods = np.array(['cPCH-Fisher-t', 'cPCH-Simes-t', 'cPCH-Bonferroni-t', 'cPCH-Fisher-norm', 'cPCH-Simes-norm', 'cPCH-Bonferroni-norm'])
dof_len = len(dof_vec)

data_df = pd.DataFrame({'m': np.repeat(m, dof_len*6), 'r': np.repeat(r, dof_len*6), 'k': np.repeat(k, dof_len*6),
                        'sig_strength': np.repeat(sig_strength, dof_len*6), 'dof': np.tile(dof_vec, 6),
                    'methods': np.repeat(methods, dof_len),
                        'p_reject': np.array(power).reshape(len(dof_vec)*6),
                        'se': np.array(se).reshape(len(dof_vec)*6)})

#change path to where and how you want the data to be saved
csv_name = "Data/robustness/m_" + str(m) + "/m_" + str(m) + "_r_" + str(r) + '_k_' + str(k) + '_ss_' + str(sig_strength) + '.csv'
data_df.to_csv(csv_name, index = False)
