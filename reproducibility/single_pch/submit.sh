for m in 3 4 5; do
  mkdir Data/m_${m}
for r in $(seq 2 $m); do
for k in $(seq 0 $m); do
for sig_strength in 0 1 2 3 4; do
for alpha in 0.05; do
for replicates in 5000; do
for N in 10000; do

 echo "${m} ${r} ${k} ${sig_strength} ${alpha} ${replicates} ${N}"
export m r k sig_strength alpha replicates N


sbatch -o output/out_m_${m}_r_${r}_k_${k}_ss_${sig_strength}.stdout.txt \
-e err/err_m_${m}_r_${r}_k_${k}_ss_${sig_strength}.stdout.txt \
--job-name="cPCH_single_m_${m}_r_${r}_k_${k}_ss_${sig_strength}" \
batch.sh
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
