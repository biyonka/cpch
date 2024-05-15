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
                    help='the number of iterations to use for the MC cPCH sampling')
a = parser.parse_args()
print(a)

m = a.m
r = a.r
k = a.k
sig_strength = a.sig_strength
alpha = a.alpha
replicates = a.replicates
N = a.N
mus = np.array([sig_strength])


#set a unique seed per configuration, making it so each r is evaluated on
#same dataset with config n, k and ss
np.random.seed(int(20000*m + 1000*k + 10*sig_strength))


adaFilter = importr('adaFilter')
af = adaFilter.adaFilter
classical_methods = adaFilter.ClassicalMTP

rpy2.robjects.numpy2ri.activate()

def single_pc_power(m, r, k, alpha, sig_strength, replicates, N, tau = 0.1, save_XX = False):
    #generate 10000 replicates of normal vector with k having mean mu and n-k having mean 0
    mu_vec = np.append(np.repeat(sig_strength, k), np.repeat(0, m-k))
    XX = norm.rvs(mu_vec, scale = 1, size = (replicates, m))
    #get individual p-values from it
    XX_pvals = 2 * (1-norm.cdf(np.abs(XX)))

    cpch_f, cpch_f_oracle, fisher = [], [], []
    cpch_s, cpch_s_oracle, simes = [], [], []
    cpch_b, cpch_b_oracle, bonferroni= [], [], []
    # # #run each of the classical methods
    nr, nc = XX_pvals.shape
    Br = ro.r.matrix(XX_pvals, nrow=nr, ncol=nc)
    fisher_results = classical_methods(Br, r = r, alpha = alpha, method = "Fisher")
    simes_results = classical_methods(Br, r = r, alpha = alpha, method = "Simes")
    bon_results = classical_methods(Br, r = r, alpha = alpha, method = "Bonferroni")
    fisher_pvals = np.array(list(pd.DataFrame(fisher_results).to_numpy()[1, :]))
    simes_pvals = np.array(list(pd.DataFrame(simes_results).to_numpy()[1, :]))
    bon_pvals = np.array(list(pd.DataFrame(bon_results).to_numpy()[1, :]))

    cpch_pvals_f = cpch(XX, m, r, f_fisher, norm.pdf, norm.cdf, truncnorm.rvs, N)
    cpch_pvals_s = cpch(XX, m, r, f_simes,  norm.pdf, norm.cdf, truncnorm.rvs, N)
    cpch_pvals_b = cpch(XX, m, r, f_bonferroni, norm.pdf, norm.cdf, truncnorm.rvs, N)

    cpch_pvals_f_oracle = cpch_oracle(XX, m, r, f_fisher, np.tile(mu_vec, replicates).reshape(replicates, m), norm.pdf, norm.cdf, truncnorm.rvs, N)
    cpch_pvals_s_oracle = cpch_oracle(XX, m, r, f_simes,  np.tile(mu_vec, replicates).reshape(replicates, m), norm.pdf, norm.cdf, truncnorm.rvs, N)
    cpch_pvals_b_oracle = cpch_oracle(XX, m, r, f_bonferroni,  np.tile(mu_vec, replicates).reshape(replicates, m), norm.pdf, norm.cdf, truncnorm.rvs, N)

    p_reject = np.array([np.mean(fisher_pvals <= alpha), np.mean(simes_pvals <= alpha), np.mean(bon_pvals <= alpha),
          np.mean(cpch_pvals_f <= alpha), np.mean(cpch_pvals_s <= alpha), np.mean(cpch_pvals_b <= alpha),
    np.mean(cpch_pvals_f_oracle <= alpha), np.mean(cpch_pvals_s_oracle <= alpha), np.mean(cpch_pvals_b_oracle <= alpha)])

    ses = np.array([np.std(fisher_pvals <= alpha, ddof = 1), np.std(simes_pvals <= alpha, ddof = 1), np.std(bon_pvals <= alpha, ddof=1),
          np.std(cpch_pvals_f <= alpha, ddof=1), np.std(cpch_pvals_s <= alpha, ddof = 1), np.std(cpch_pvals_b <= alpha, ddof = 1),
    np.std(cpch_pvals_f_oracle <= alpha, ddof = 1), np.std(cpch_pvals_s_oracle <= alpha, ddof=1), np.std(cpch_pvals_b_oracle <= alpha, ddof = 1)])/(replicates**0.5)

    methods = np.array(['Fisher', 'Simes', 'Bonferroni', 'cPCH-Fisher', 'cPCH-Simes', 'cPCH-Bonferroni','cPCH-Oracle-Fisher', 'cPCH-Oracle-Simes', 'cPCH-Oracle-Bonferroni'])

    data_df = pd.DataFrame({'m': np.repeat(m, 9), 'r': np.repeat(r, 9), 'k': np.repeat(k, 9), 'sig_strength': np.repeat(sig_strength, 9),
                        'methods': methods, 'p_reject': p_reject, 'se': ses})
    return data_df

data_df = single_pc_power(m, r, k, alpha, sig_strength, replicates, N, tau = 0.1, save_XX = False)
#change path to where and how you want the data to be saved
csv_name = "Data/m_" + str(m) + "/m_" + str(m) + "_r_" + str(r) + '_k_' + str(k) + '_ss_' + str(sig_strength) + '.csv'
data_df.to_csv(csv_name, index = False)
