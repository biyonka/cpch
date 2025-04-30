#replace path with location of cpch_source.py locally
exec(open("../cpch_source.py").read())
###################
#### Set Params ###
###################
parser = argparse.ArgumentParser(description='Calculate Type I error for each theta1 theta2 pair when m=3, r=2, and r*=2')
parser.add_argument('theta1', metavar='theta1', type=float,
                    help='theta1 of test stat')
parser.add_argument('theta2', metavar='theta2', type=float,
                    help='theta2 of test stat')
parser.add_argument('alpha', metavar='alpha', type=float,
                    help='the alpha level of the PC test')
parser.add_argument('replicates', metavar='replicates', type=int,
                    help='the number of iterations to run for the t1 error and power estimation')
parser.add_argument('N', metavar='N', type=int,
                    help='the number of iterations to use for the MC cPCH sampling')
a = parser.parse_args()
print(a)

m = 3
r = 3
theta1=a.theta1
theta2=a.theta2
#k = a.k
#sig_strength = a.sig_strength
alpha = a.alpha
#rho = a.rho
replicates = a.replicates
N = a.N
#mus = np.array([sig_strength])

#set a unique seed per configuration, making it so each r is evaluated on
#same dataset with config n, k and ss
#np.random.seed(int(20000*theta1 + 1000*theta2))

np.random.seed(int(90000*theta1 + 1000*theta2))
#for m=3, r=3, r*=2, for each theta1, theta2 value, generate many test statistics sampled from those values
#calculate cPCH p-value across the many test statistics and get estimate of Type I error
#for m=3, r=2, r*=1, for each theta1 value, do the same thing?
def single_pc_power(theta1, theta2, alpha, replicates, N=10000):
    #generate 10000 replicates of normal vector with k having mean mu and n-k having mean 0
    mu_vec = np.array([theta1, theta2, 0])
    if theta1==theta2:
        XX = norm.rvs(mu_vec, scale = 1, size = (2*replicates, m)) #double replicates for diagonal since off diagonals have symmetric points
    else:
        XX = norm.rvs(mu_vec, scale = 1, size = (replicates, m))

    cpch_pvals_f = cpch_unadjusted(XX, m, r, f_fisher, norm.pdf, norm.cdf, truncnorm.rvs, N)
    cpch_pvals_s = cpch_unadjusted(XX, m, r, f_simes,  norm.pdf, norm.cdf, truncnorm.rvs, N)

    total_reject = np.array([#np.mean(fisher_pvals <= alpha), np.mean(simes_pvals <= alpha), #np.mean(bon_pvals <= alpha),
          np.sum(cpch_pvals_f <= alpha), np.sum(cpch_pvals_s <= alpha)#, np.mean(cpch_pvals_b <= alpha),
   # np.mean(cpch_pvals_f_oracle <= alpha), np.mean(cpch_pvals_s_oracle <= alpha)# np.mean(cpch_pvals_b_oracle <= alpha)
    ]) #in the end, combine total_reject across reflection to double sample size

    ses = np.array([#np.std(fisher_pvals <= alpha, ddof = 1), np.std(simes_pvals <= alpha, ddof = 1),# np.std(bon_pvals <= alpha, ddof=1),
          np.std(cpch_pvals_f <= alpha, ddof=1), np.std(cpch_pvals_s <= alpha, ddof = 1)#, np.std(cpch_pvals_b <= alpha, ddof = 1),
   # np.std(cpch_pvals_f_oracle <= alpha, ddof = 1), np.std(cpch_pvals_s_oracle <= alpha, ddof=1), np.std(cpch_pvals_b_oracle <= alpha, ddof = 1)
    ])/((replicates)**0.5) #se of each pixel, can just divide by sqrt{2} to get estimate for full sample size later on

    methods = np.array(['cPCH-Fisher', 'cPCH-Simes'])

    data_df = pd.DataFrame({'theta1': theta1, 'theta2': theta2, 
                        'methods': methods, 'total_reject': total_reject, 'se': ses})
    return data_df

data_df = single_pc_power(theta1, theta2, alpha, replicates, N)
#change path to where and how you want the data to be saved
csv_name = "Data_t1err_grid_2/theta1_" + str(theta1) + "_theta2_" + str(theta2) + '.csv'
data_df.to_csv(csv_name, index = False)
