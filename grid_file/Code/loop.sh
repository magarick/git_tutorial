for ITER in {1..2}; do
  echo "running $ITER"
  echo "export ITER=$ITER; echo "running ITER: $ITER TASK: \$SGE_TASK_ID"; R --no-save < loop.r" | qsub -N moms-$ITER -j y -t 1-2 -cwd
  sleep 1
done
