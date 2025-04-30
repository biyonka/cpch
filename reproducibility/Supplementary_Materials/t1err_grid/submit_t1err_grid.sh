
mkdir -p Data_t1err_grid_2
for theta1 in $(seq 0 0.2 5); do
for theta2 in $(seq 0 0.2 5); do
for alpha in 0.05; do
for replicates in 10000; do
for N in 20000; do

 echo "${theta1} ${theta2} ${alpha} ${replicates} ${N}"
export theta1 theta2 alpha replicates N


sbatch -o output/out_t1errgrid_theta1_${theta1}_theta2_${theta2}.stdout.txt \
-e err/err_t1errgrid_theta1_${theta1}_theta2_${theta2}.stdout.txt \
--job-name="t1errgrid_theta1_${theta1}_theta2_${theta2}" \
batch_t1err_grid.sh
#
#
sleep 1 # pause to be kind to the scheduler
done
done
done
done
done
