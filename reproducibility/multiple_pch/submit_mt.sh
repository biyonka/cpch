for m in 4; do
for M in 2000; do
  #change to where you want data to be saved
  mkdir Data/M_${M}_m_${m}
  mkdir Data/M_${M}_m_${m}/final_data
for r in $(seq 2 $m); do
for pi1 in 0.05 0.1 0.4 0.8; do
for p in 0.0 0.5 1.0; do
for sig_strength in 2 3; do
    #change to where you want data to be saved
    mkdir Data/M_${M}_m_${m}/r_${r}_pi1_${pi1}_p_${p}_ss_${sig_strength}
for alpha in 0.1; do
for N in 10000; do
for iters in $(seq 1 50); do

 echo "${m} ${M} ${r} ${pi1} ${p} ${sig_strength} ${alpha} ${iters} ${N}"
export m r M pi1 p sig_strength alpha iters N

#change directory to where and how you want output and error to be written
sbatch -o output/cpchmt__output_iter_${iters}_r_${r}_pi1_${pi1}_p_${p}_ss_${sig_strength}_M_${M}_m_${m}_alpha_${alpha}.stdout.txt \
-e err/cpchmt__errors_iter_${iters}_r_${r}_pi1_${pi1}_p_${p}_ss_${sig_strength}_M_${M}_m_${m}_alpha_${alpha}.stdout.txt \
--job-name="cpchmt_iter_${iters}_r_${r}_pi1_${pi1}_p_${p}_ss_${sig_strength}_M_${M}_m_${m}" \
batch_mt.sh
#
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
done
