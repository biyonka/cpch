for m in 3; do
  mkdir -p Data_dependence/m_${m}
for r in $(seq 2 $m); do
for k in $(seq 0 $m); do
for sig_strength in 0 1 2 3 4; do
for alpha in 0.05; do
for rho in 0.1 0.25 0.5 0.75 0.9; do
for replicates in 5000; do
for N in 10000; do

 echo "${m} ${r} ${k} ${sig_strength} ${alpha} ${rho} ${replicates} ${N}"
export m r k sig_strength alpha rho replicates N


sbatch -o output/out_m_${m}_r_${r}_k_${k}_ss_${sig_strength}_rho_${rho}.stdout.txt \
-e err/err_m_${m}_r_${r}_k_${k}_ss_${sig_strength}_rho_${rho}.stdout.txt \
--job-name="cPCH_single_m_${m}_r_${r}_k_${k}_ss_${sig_strength}_rho_${rho}" \
batch_dependence.sh
#
#
sleep 1 # pause to be kind to the scheduler
done
done
done
done
done
done
done
done