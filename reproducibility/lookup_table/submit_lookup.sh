for m in 2 3 4 5; do
  mkdir -p Lookup_table/m_${m}
for r in $(seq 2 $m); do
for alpha in 0.001 0.005 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.12 0.14 0.16 0.18 0.20; do
for N in 100000; do

 echo "${m} ${r} ${alpha} ${N}"
export m r k alpha N

sbatch -o output/out_m_${m}_r_${r}_alpha_${alpha}.stdout.txt \
-e err/err_m_${m}_r_${r}_alpha_${alpha}.stdout.txt \
--job-name="cPCH_single_m_${m}_r_${r}_alpha_${alpha}" \
batch_lookup.sh
#
#
sleep 1 # pause to be kind to the scheduler
done
done
done
done



# for m in 3 4 5; do
#   mkdir -p Lookup_table/m_${m}
# for r in $(seq 2 $m); do
# for alpha in 0.001; do
# for N in 2000; do

#  echo "${m} ${r} ${alpha} ${N}"
# export m r k alpha N

# sbatch -o output/out_m_${m}_r_${r}_alpha_${alpha}.stdout.txt \
# -e err/err_m_${m}_r_${r}_alpha_${alpha}.stdout.txt \
# --job-name="cPCH_single_m_${m}_r_${r}_alpha_${alpha}" \
# batch_lookup.sh
# #
# #
# sleep 1 # pause to be kind to the scheduler
# done
# done
# done
# done


# for m in 2; do
#   mkdir -p Lookup_table/m_${m}
# for r in 2; do
# for alpha in 0.02 0.04 0.06 0.08; do # 0.001 0.01 0.03 0.05 0.07 0.09 0.10 0.12 0.14 0.16 0.18 0.20; do
# for N in 500; do

#  echo "${m} ${r} ${alpha} ${N}"
# export m r alpha N


# sbatch -o output/out_m_${m}_r_${r}_alpha_${alpha}.stdout.txt \
# -e err/err_m_${m}_r_${r}_alpha_${alpha}.stdout.txt \
# --job-name="cPCH_single_m_${m}_r_${r}_alpha_${alpha}" \
# batch_lookup.sh
# #
# #
# sleep 1 # pause to be kind to the scheduler
# done
# done
# done
# done