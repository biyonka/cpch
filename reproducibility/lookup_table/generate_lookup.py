exec(open("../cpch_source.py").read())
###################
#### Set Params ###
###################
parser = argparse.ArgumentParser(description='Generate grid plot given some n')
parser.add_argument('m', metavar='m', type=int,
                    help='the number of studies for the PC test')
parser.add_argument('r', metavar='r', type=int,
                            help='replicability number')
parser.add_argument('alpha', metavar='alpha', type=float,
                    help='the alpha level of the PC test')
parser.add_argument('N', metavar='N', type=int,
                    help='the number of iterations to use for the MC cPCH sampling')
#parser.add_argument('inputs', metavar='a', type=int, nargs='+',
                    #help='n, power_iters, N')
a = parser.parse_args()
print(a)

m = a.m
r = a.r
#k = a.k
#sig_strength = a.sig_strength
alpha = a.alpha
#threshold = a.threshold
#replicates = a.replicates
N = a.N
#mus = np.array([sig_strength])


#set a unique seed per configuration, making it so each r is evaluated on
#same dataset with config n, k and ss
np.random.seed(int(round(1000*alpha)))#int(20000*m + 1000*k + 10*sig_strength))

#packnames = ('devtools', 'BiocManager')
#devtools = importr('devtools')
#biohub = importr('BiocManager')
#adaFilter = importr('adaFilter')
#af = adaFilter.adaFilter
#classical_methods = adaFilter.ClassicalMTP

#rpy2.robjects.numpy2ri.activate()

#get appropriate alpha'
# %run cpch_new_method_utils.py
# %run src.py
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from tqdm.notebook import tqdm_notebook
print(cpu_count())

cpus = cpu_count()
#alpha = 0.05

from numpy.linalg import norm as npnorm    
from scipy.optimize import lsq_linear
from scipy.linalg import toeplitz
from numpy.linalg import inv


def lsq(y):
    y = np.array(y)
    d = y.shape[0]
    ub = np.array([np.inf]*d)
    lb = np.zeros(d)
    A = toeplitz([1] + [0]*(d-1), [1,-1] +[0]*(d-2))
    Ainv = inv(A)
    beta = lsq_linear(Ainv, y, bounds=(lb, ub), lsmr_tol='auto', verbose=0)
    return Ainv.dot(beta.x)    

def proj(v):
    v = np.array(v)
    if len(v) == 1:
        return np.maximum(v, 0)
    elif len(v) == 2:
        x, y = v
        if y <= 0:
            return np.array([max(x, 0), 0])
        elif y <= -x:
            return np.zeros(2)
        elif y <= x:
            return v
        else:
            e0 = np.ones(2)/np.sqrt(2)
            return v.dot(e0)*e0
    elif len(v) > 2:
        return lsq(v)
def gen_init(m, r, scale=5):
    arr = np.zeros(m)
    for j in range(r-1):
        if j == 0:
            arr[j] = np.random.rand(1)*scale
        else:
            arr[j] = np.random.rand(1)*arr[j-1]
    return arr

def get_pvals_core(Z_hold, mu):
    XX = Z_hold + mu
    return cpch(XX, m, r, f, pdf=norm.pdf, cdf=norm.cdf, N=N)
def get_grad_core(Z_hold, pvals_hold):
    return  -(Z_hold)[:,:(r-1)]*(pvals_hold <= alpha)[:, np.newaxis]

def sgd(m, r, f, alpha, N=500):
    mu_cur = gen_init(m=m, r=r)
    #print(f'init at mu={mu_cur}')
    mu_seq = []
    err_seq = []
    gradnorm_seq = [] 
    Tmax = 200
    lrate0 = 30
    lrate_seq = [max(lrate0/(1 + 0.02*tt), 2) for tt in range(Tmax)]
    K_seq = []
    # lrate_seq = [lrate0*0.95**tt for tt in range()]
    for tt in range(Tmax):
    ## exponential decay
        lrate = lrate_seq[tt]
        if tt % 10 == 0:
            print(f'epoch: {tt}, learn rate: {lrate}')
        if tt < 100:
            K = 1e3
            # K = 96*1280
        elif tt < 150:
            K = 1e4
        else:
            K = 3e4#1e5
        K = int(K)
        K_seq.append(K)
        ## shuffle for mini batch
        # Z_mini = normal(scale=1, size=[K, n])
        # X_mini = Z_mini + mu_cur
        # pvals = get_pvals_core(Z_mini, mu_cur)
        # grads = get_grad(n, r, X_mini, mu_cur, pvals, alpha)

        Z_list = [normal(scale=1, size=[K//cpus, m]) for _ in range(cpus)]
        Z_all = np.vstack(Z_list)
        star_args = [[Z_list[o], mu_cur] for o in range(cpus)]
        with Pool(processes=cpus) as pool:
            pvals= np.hstack(pool.starmap(get_pvals_core, star_args))
        grads = get_grad_core(Z_all, pvals)

        diff = lrate*grads.mean(axis=0)
        mu_cur[:(r-1)] = proj(mu_cur[:(r-1)] - diff)
        mu_seq.append(mu_cur[:(r-1)].copy())
        err_seq.append((pvals <= alpha).mean())
        gradnorm_seq.append(npnorm(grads.mean(axis=0), np.inf))
    mu_seq = np.vstack(mu_seq)
    print('final epoch error:', err_seq[-1])
    print('se cal:', np.sqrt(1/20/K))
    print('max type I error:', err_seq[-1], 'at mu =', mu_seq[-1,:])
    return(err_seq[-1])

def binary_search(m, r, f, nominal_alpha, threshold, low_alpha, high_alpha):
    # Check base case
    if high_alpha >= low_alpha:
 
        mid_alpha = (high_alpha + low_alpha) / 2
        
        t1_err_mid = sgd(m, r, f, mid_alpha, N=500)
        # If the max t1 error at the midpoint is within some epsilon ball of the desired nominal level
        if np.abs(t1_err_mid - nominal_alpha) <= threshold:
            return [mid_alpha, t1_err_mid]
 
        # If the max t1 error at the midpoint is too far from nominal on the high end
        elif t1_err_mid - nominal_alpha > threshold:
            print('new alpha range: '+str(low_alpha) +' , '+ str(mid_alpha))
            return binary_search(m, r, f, nominal_alpha, threshold, low_alpha, mid_alpha)
 
        # if max t1 error at midpoint is too far from nominal on the lower end
        elif nominal_alpha - t1_err_mid > threshold:
            print('new alpha range: '+str(mid_alpha) +' , '+ str(high_alpha))
            return binary_search(m, r, f, nominal_alpha, threshold, mid_alpha, high_alpha)
    else:
        return (-1)

alpha_primes = []
t1_err_at_alpha_prime = []
if alpha >= 0.04:
    f=f_fisher
    alpha_prime_fisher = binary_search(m=m, r=r, f=f_fisher, nominal_alpha=alpha, threshold=0.004, low_alpha=alpha-0.03, high_alpha=alpha)
    print("Fisher alpha'" + str(alpha_prime_fisher[0]))
    alpha_primes.append(alpha_prime_fisher[0])
    t1_err_at_alpha_prime.append(alpha_prime_fisher[1])
    f=f_simes
    alpha_prime_simes = binary_search(m=m, r=r, f=f_simes, nominal_alpha=alpha, threshold=0.004, low_alpha=alpha-0.03, high_alpha=alpha)
    print("Simes alpha'" + str(alpha_prime_simes[0]))
    alpha_primes.append(alpha_prime_simes[0])
    t1_err_at_alpha_prime.append(alpha_prime_simes[1])
elif alpha < 0.04 & alpha >=0.005:
    f=f_fisher
    alpha_prime_fisher = binary_search(m=m, r=r, f=f_fisher, nominal_alpha=alpha, threshold=alpha/3, low_alpha=alpha/2, high_alpha=alpha)
    print("Fisher alpha'" + str(alpha_prime_fisher[0]))
    alpha_primes.append(alpha_prime_fisher[0])
    t1_err_at_alpha_prime.append(alpha_prime_fisher[1])
    f=f_simes
    alpha_prime_simes = binary_search(m=m, r=r, f=f_simes, nominal_alpha=alpha, threshold=alpha/3, low_alpha=alpha/2, high_alpha=alpha)
    print("Simes alpha'" + str(alpha_prime_simes[0]))
    alpha_primes.append(alpha_prime_simes[0])  
    t1_err_at_alpha_prime.append(alpha_prime_simes[1])
elif alpha < 0.005:
    f=f_fisher
    N=2000 #need to set N larger to make sure p-values can attain values below alpha
    alpha_prime_fisher = binary_search(m=m, r=r, f=f_fisher, nominal_alpha=alpha, threshold=alpha/3, low_alpha=alpha/2, high_alpha=alpha)
    print("Fisher alpha'" + str(alpha_prime_fisher[0]))
    alpha_primes.append(alpha_prime_fisher[0])
    t1_err_at_alpha_prime.append(alpha_prime_fisher[1])
    f=f_simes
    alpha_prime_simes = binary_search(m=m, r=r, f=f_simes, nominal_alpha=alpha, threshold=alpha/3, low_alpha=alpha/2, high_alpha=alpha)
    print("Simes alpha'" + str(alpha_prime_simes[0]))
    alpha_primes.append(alpha_prime_simes[0])  
    t1_err_at_alpha_prime.append(alpha_prime_simes[1])

data_df = pd.DataFrame({'m': np.repeat(m, 2), 'r': np.repeat(r, 2), 'alpha': np.repeat(alpha, 2),
                        'methods': np.array(['Fisher', 'Simes']), 'alpha_prime': alpha_primes, 't1_err_alpha_prime': t1_err_at_alpha_prime})


csv_name = "Lookup_table/m_" + str(m) + "/m_" + str(m) + "_r_" + str(r) + '_alpha_' + str(alpha) + '.csv'
data_df.to_csv(csv_name, index = False)
