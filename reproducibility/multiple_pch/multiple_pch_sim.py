exec(open("../cpch_source.py").read())
###################
#### Set Params ###
###################
parser = argparse.ArgumentParser(description='Generate grid plot given some m')
parser.add_argument('m', metavar='m', type=int,
                    help='the number of studies for the PC test')
parser.add_argument('r', metavar='r', type=int,
                            help='replicability number')
parser.add_argument('M', metavar='M', type=int,
                            help='number of PC hypotheses')
parser.add_argument('pi1', metavar='pi1', type=float,
                    help='true proportion of alternative')
parser.add_argument('p', metavar='p', type=float,
                    help='probability we identify a hypothesis correctly, or resample from bernoulli, the folly of man')
parser.add_argument('sig_strength', metavar='sig_strength', type=int,
                    help='the multiplier on the mu vector')
parser.add_argument('alpha', metavar='alpha', type=float,
                    help='the alpha level of the PC test')
parser.add_argument('iters', metavar='iters', type=int,
                    help='the number of iterations to run a sample of power_iters hypotheses power estimation')
parser.add_argument('N', metavar='N', type=int,
                    help='the number of iterations to use for the MC cPCH sampling')
a = parser.parse_args()
print(a)

m = a.m
r = a.r
p = a.p
pi1 = a.pi1
M = a.M
sig_strength = a.sig_strength
alpha = a.alpha
iteration = a.iters
N = a.N

mus = np.array([sig_strength])
power_iters = 20

#set a unique seed per configuration
np.random.seed(int(100000000*pi1 +100000*p + 1000*sig_strength + 100*r + iteration))
adaFilter = importr('adaFilter')
af = adaFilter.adaFilter
classical_methods = adaFilter.ClassicalMTP

rpy2.robjects.numpy2ri.activate()


def sim_study_general_n_r(pi1, mus, p, m, r, M, alpha, power_iters, N, tau = 0.1, save_XX = False):
    XX, null_nonnull = generate_XX_mt(pi1, p, mus, m, r, M, power_iters)
    # 1 if null, 0 alt
    null_bool = np.sum(null_nonnull, axis = 1) < r
    #need this chunk to prevent division by zero errors caused by every hypothesis being null
    if sum(1-null_bool) == 0:
        while sum(1-null_bool) == 0:
            XX, null_nonnull = generate_XX_mt(pi1, p, mus, m, r, M, power_iters)
            null_bool = np.sum(null_nonnull, axis = 1) < r

    XX_reshaped = XX.reshape(XX.shape[0], -1)
    if save_XX:
        filename = "Data/XX_m_" + str(m) + '_r_' + str(r) + '_M_' + str(M) + '.txt'
        #save data in case we want to access it later
        np.savetxt(filename, XX_reshaped)

    alts = np.logical_not(null_bool)
    cpch_f_bh_fdr, cpch_f_storey_fdr, f_bh_fdr, f_storey_fdr, f_dhh_fdr, f_sdhh_fdr, ada_fdr = [], [], [], [], [],[],[]
    cpch_s_bh_fdr, cpch_s_storey_fdr, s_bh_fdr, s_storey_fdr, s_dhh_fdr, s_sdhh_fdr = [], [], [], [], [],[]
    cpch_f_bh_power, cpch_f_storey_power, f_bh_power, f_storey_power, f_dhh_power, f_sdhh_power, ada_power = [],[],[],[],[],[],[]
    cpch_s_bh_power, cpch_s_storey_power, s_bh_power, s_storey_power, s_dhh_power, s_sdhh_power = [],[],[],[],[],[]
    # # #run each of the classical methods
    XX_pvals = 2 * (1-norm.cdf(np.abs(XX)))
    for i in np.arange(power_iters):
        XXi = XX[i]
        XX_pvals_i = XX_pvals[i]
        cpch_pvals_f = cpch(XXi, m, r, f_fisher, norm.pdf, norm.cdf, truncnorm.rvs, N)
        cpch_pvals_s = cpch(XXi, m, r, f_simes, norm.pdf, norm.cdf, truncnorm.rvs, N)
        #cpch_pvals_b = cpch(XXi, m, r, f_bonferroni, norm.pdf, norm.cdf, truncnorm.rvs, N)
        #adafilter
        nr, nc = XX_pvals_i.shape
        Br = ro.r.matrix(XX_pvals_i, nrow=nr, ncol=nc)
        ada_results = af(Br, r = r, alpha = alpha)
        ada_decisions = pd.DataFrame(ada_results).to_numpy()[0,:]

        #cpch-Fisher-BH
        cpch_f_bh_decisions = mt.multipletests(cpch_pvals_f, alpha, method = 'fdr_bh')[0]
        #cpch-Fisher-Storey
        cpch_f_storey_decisions = storey(cpch_pvals_f, alpha)

        #cpch-Simes-BH
        cpch_s_bh_decisions = mt.multipletests(cpch_pvals_s, alpha, method = 'fdr_bh')[0]
        #cpch-Simes-Storey
        cpch_s_storey_decisions = storey(cpch_pvals_s, alpha)

        # #cpch-Bonferroni-BH
        # cpch_b_bh_decisions = mt.multipletests(cpch_pvals_b, alpha, method = 'fdr_bh')[0]
        # #cpch-Bonferroni-Storey
        # cpch_b_storey_decisions = storey(cpch_pvals_b, alpha)

        ##Fisher-BH
        fisher_results = classical_methods(Br, r = r, alpha = alpha, method = "Fisher")
        f_bh_decisions = pd.DataFrame(fisher_results).to_numpy()[0, :]

        ##Simes-BH
        simes_results = classical_methods(Br, r = r, alpha = alpha, method = "Simes")
        s_bh_decisions = pd.DataFrame(simes_results).to_numpy()[0, :]

        # ##Bonferroni-BH
        # bon_results = classical_methods(Br, r = r, alpha = alpha, method = "Bonferroni")
        # b_bh_decisions = pd.DataFrame(bon_results).to_numpy()[0, :]

        ##Fisher-Storey
        fisher_pvals = np.array(list(pd.DataFrame(fisher_results).to_numpy()[1, :]))
        f_storey_decisions = storey(fisher_pvals, alpha)

      ##Simes-Storey
        simes_pvals = np.array(list(pd.DataFrame(simes_results).to_numpy()[1, :]))
        s_storey_decisions = storey(simes_pvals, alpha)

        shortened_fisher = fisher_pvals[fisher_pvals <= tau]/tau
        shortened_simes = simes_pvals[simes_pvals <= tau]/tau
       # shortened_bon = bon_pvals[bon_pvals <= tau]/tau
        if len(shortened_fisher) >= 1:
            ##Fisher-DHH
            f_dhh_decisions = mt.multipletests(shortened_fisher, alpha, method = 'fdr_bh')[0]
            shortened_null_bool_f = null_bool[fisher_pvals <= tau]
            ##Fisher-sDHH
            f_sdhh_decisions = storey(shortened_fisher, alpha)

            f_dhh_fdr.append(sum(f_dhh_decisions[shortened_null_bool_f])/max(1, sum(f_dhh_decisions)))
            f_sdhh_fdr.append(sum(f_sdhh_decisions[shortened_null_bool_f])/max(1, sum(f_sdhh_decisions)))
            f_dhh_power.append(sum(f_dhh_decisions[np.logical_not(shortened_null_bool_f)])/sum(1-null_bool))
            f_sdhh_power.append(sum(f_sdhh_decisions[np.logical_not(shortened_null_bool_f)])/sum(1-null_bool))
        else: #none of the p-values were less than tau, so dont reject anything
            f_dhh_fdr.append(0)
            f_sdhh_fdr.append(0)
            f_dhh_power.append(0)
            f_sdhh_power.append(0)

        if len(shortened_simes) >= 1:
            ##Simes-DHH
            s_dhh_decisions = mt.multipletests(shortened_simes, alpha, method = 'fdr_bh')[0]
            shortened_null_bool_s = null_bool[simes_pvals <= tau]
            ##Simes-sDHH
            s_sdhh_decisions = storey(shortened_simes, alpha)

            s_dhh_fdr.append(sum(s_dhh_decisions[shortened_null_bool_s])/max(1, sum(s_dhh_decisions)))
            s_sdhh_fdr.append(sum(s_sdhh_decisions[shortened_null_bool_s])/max(1, sum(s_sdhh_decisions)))
            s_dhh_power.append(sum(s_dhh_decisions[np.logical_not(shortened_null_bool_s)])/sum(1-null_bool))
            s_sdhh_power.append(sum(s_sdhh_decisions[np.logical_not(shortened_null_bool_s)])/sum(1-null_bool))
        else: #none of the p-values were less than tau, so dont reject anything
            s_dhh_fdr.append(0)
            s_sdhh_fdr.append(0)
            s_dhh_power.append(0)
            s_sdhh_power.append(0)

        cpch_f_bh_fdr.append(np.sum(cpch_f_bh_decisions[null_bool])/max(1, np.sum(cpch_f_bh_decisions)))
        cpch_f_storey_fdr.append(np.sum(cpch_f_storey_decisions[null_bool])/max(1, np.sum(cpch_f_storey_decisions)))
        f_bh_fdr.append(np.sum(f_bh_decisions[null_bool])/max(1, np.sum(f_bh_decisions)))
        f_storey_fdr.append(np.sum(f_storey_decisions[null_bool])/max(1, np.sum(f_storey_decisions)))

        cpch_s_bh_fdr.append(np.sum(cpch_s_bh_decisions[null_bool])/max(1, np.sum(cpch_s_bh_decisions)))
        cpch_s_storey_fdr.append(np.sum(cpch_s_storey_decisions[null_bool])/max(1, np.sum(cpch_s_storey_decisions)))
        s_bh_fdr.append(np.sum(s_bh_decisions[null_bool])/max(1, np.sum(s_bh_decisions)))
        s_storey_fdr.append(np.sum(s_storey_decisions[null_bool])/max(1, np.sum(s_storey_decisions)))

        ada_fdr.append(np.sum(ada_decisions[null_bool])/max(1, np.sum(ada_decisions)))

        #if num_alt > 0:
        cpch_f_bh_power.append(np.mean(cpch_f_bh_decisions[alts]))
        cpch_f_storey_power.append(np.mean(cpch_f_storey_decisions[alts]))
        f_bh_power.append(np.mean(f_bh_decisions[alts]))
        f_storey_power.append(np.mean(f_storey_decisions[alts]))

        cpch_s_bh_power.append(np.mean(cpch_s_bh_decisions[alts]))
        cpch_s_storey_power.append(np.mean(cpch_s_storey_decisions[alts]))
        s_bh_power.append(np.mean(s_bh_decisions[alts]))
        s_storey_power.append(np.mean(s_storey_decisions[alts]))


        ada_power.append(np.mean(ada_decisions[alts]))
    #if num_alt > 0:
    return np.array([cpch_f_bh_fdr, cpch_s_bh_fdr, 
                     cpch_f_storey_fdr, cpch_s_storey_fdr,
                     f_bh_fdr, s_bh_fdr,
                     f_storey_fdr, s_storey_fdr,
                     f_dhh_fdr, s_dhh_fdr,
                     f_sdhh_fdr, s_sdhh_fdr,
                     ada_fdr,
                    cpch_f_bh_power, cpch_s_bh_power,
                    cpch_f_storey_power, cpch_s_storey_power,
                    f_bh_power, s_bh_power,
                    f_storey_power, s_storey_power, 
                    f_dhh_power, s_dhh_power,
                    f_sdhh_power, s_sdhh_power,
                     ada_power])#


data = sim_study_general_n_r(pi1, mus, p, m, r, M, alpha, power_iters, N)

data_df = pd.DataFrame(data.T)
data_df.insert(0,'m', m)
data_df.insert(1,'M', M)
data_df.insert(2,'r', r)
data_df.insert(3,'p', p)
data_df.insert(4,'sig_strength', sig_strength)
data_df.insert(5,'pi1', pi1)
data_df.columns = [ 'm', 'M', 'r', 'p', 'sig_strength', 'pi1',
'cpch_f_bh_fdr', 'cpch_s_bh_fdr',
                     'cpch_f_storey_fdr', 'cpch_s_storey_fdr',
                     'f_bh_fdr', 's_bh_fdr',
                     'f_storey_fdr', 's_storey_fdr', 
                     'f_dhh_fdr', 's_dhh_fdr',
                     'f_sdhh_fdr', 's_sdhh_fdr',
                     'ada_fdr',
                    'cpch_f_bh_power', 'cpch_s_bh_power',
                    'cpch_f_storey_power', 'cpch_s_storey_power',
                    'f_bh_power', 's_bh_power',
                    'f_storey_power', 's_storey_power', 
                    'f_dhh_power', 's_dhh_power', 
                    'f_sdhh_power', 's_sdhh_power', 
                     'ada_power']

csv_name = "Data/M_"+str(M)+"_m_"+str(m)+"/r_"+str(r)+'_pi1_'+str(pi1)+'_p_'+str(p)+'_ss_'+str(sig_strength)+"/it_"+str(iteration)+'_r_'+str(r)+'_pi1_'+str(pi1)+'_p_'+str(p)+'_ss_'+str(sig_strength)+"_m_"+str(m)+'.csv'
data_df.to_csv(csv_name, index = False)

print("python script complete")
