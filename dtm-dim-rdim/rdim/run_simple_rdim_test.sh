if [ -z $1 ]; then echo "Need to specify outdir"; exit; fi

K=2 # Number of topics
DATA=del-d6-acl # Name of test data
ODIR=${PWD}/../outputs/$1
THREADS=4

PREFIX=${PWD}/../acl_arf-test_data/ldac/${DATA}
MDATA=${PWD}/../acl_arf-test_data/paper_metadata.csv
T=17  
MEAN=1
SD=2

#export OMP_PROC_BIND=true
#export OPENBLAS_NUM_THREADS=${THREADS}
#export MKL_NUM_THREADS=${THREADS}
export KMP_AFFINITY=verbose,scatter
#export MKL_DYNAMIC=FALSE
#export OMP_NESTED=TRUE

mkdir -p $ODIR
rm -rf $ODIR/*

date +%s > $ODIR/run.timer

# Here we go...
./main \
--mode=fit \
--model=rdim \
--rng_seed=1 \
--top_obs_var=0.5 \
--top_chain_var=0.005 \
--conditioning_epsilon=12 \
--sigma_d=0.0001 \
--sigma_l=0.0001 \
--sigma_mu=0.0001 \
--alpha=0.01 \
--initial_lda_ss=${PWD}/../acl_arf-test_data/k2-initial-lda-ss.dat \
--lda_sequence_min_iter=4 \
--lda_sequence_max_iter=5 \
--save_time=-2 \
--outname=${ODIR} \
--ntopics=${K} \
--kthreads=${K} \
--corpus_prefix=${PREFIX} \
--time_resolution=1 \
--resolution=1 \
--influence_mean_years=${MEAN} \
--influence_stdev_years=${SD} \
--metafields=1919 \
--debug=1 \
--threads=${THREADS} \
--lda_max_em_iter=25 \
--initialize_lda=true \
> >(tee ${ODIR}/run.out) 2> >(tee ${ODIR}/run.err >&2)
#--checkpoint_recover=${ODIR}/../forktest/checkpoints/iter-4/ \
#--enable_checkpointing=1 \
#--checkpoint_outdir=${ODIR}/checkpoints/ \
#--read_hdf_term1=/home/gerow/files/influence_models/influence/dtm-dim/outputs/eq32_term1.h5 \

# Skip after this point for dev.
# But do the rest for checking sanity of results
# exit(0)

##### If any parameters are changed, swap out the --initial_lda_ss otion for these two. Will change results slightly with each run:

#### Some other flags for the model:

date +%s >> $ODIR/run.timer
{ tr "\n" - < $ODIR/run.timer ; echo 0; } | bc >> $ODIR/run.timer;
echo "Total seconds taken (* -1):"
tail -n1 $ODIR/run.timer

echo "Running post-hoc analysis"
cp ../../examine-acl.R ${ODIR}/
cp ../../examine_influence-corr-acl.R ${ODIR}/
cd ${ODIR}
Rscript examine-acl.R 20 ${T} ${PREFIX} > >(tee topics.txt) 2> >(tee topics.err)
Rscript examine_influence-corr-acl.R ${K} ${PREFIX} ${MDATA} > ./influence_data.csv 2> influence.err
tail -n1 ./influence_data.csv

echo "The rank correlation above should be about 0.26"

echo "All Done"
