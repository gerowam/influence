K=$1
DATA=$2
ODIR=../outputs/$3
THREADS=$4

PREFIX=${PWD}/../acl_arf-test_data/ldac/${DATA}
MDATA=${PWD}/../acl_arf-test_data/paper_metadata.csv

#PREFIX=${PWD}/../acl_arf-full_data/ldac/${DATA}
#MDATA=${PWD}/../acl_arf-full_data/paper_metadata.csv

T=17
MEAN=1
SD=2

export OMP_PROC_BIND=true
#export GOMP_CPU_AFFINITY=0-${THREADS}
#export KMP_AFFINITY=proclist=[0-${THREADS}]

mkdir -p $ODIR
rm -rf $ODIR/*

date +%s > $ODIR/run.timer

# Here we go...
# The last two flags are new and required
#gdb main --args \
./main \
--mode=fit \
--model=fixed \
--top_obs_var=0.5 \
--top_chain_var=0.005 \
--sigma_d=0.0001 \
--sigma_l=0.0001 \
--sigma_mu=0.0001 \
--alpha=0.01 \
--lda_sequence_min_iter=20 \
--initialize_lda=true \
--lda_max_em_iter=10 \
--lda_sequence_max_iter=25 \
--save_time=-1 \
--outname=${ODIR} \
--ntopics=${K} \
--corpus_prefix=${PREFIX} \
--time_resolution=1 \
--resolution=1 \
--influence_mean_years=${MEAN} \
--influence_stdev_years=${SD} \
--rng_seed=1 \
--debug=0 \
--threads=${THREADS} \
> >(tee ${ODIR}/run.out) 2> >(tee ${ODIR}/run.err >&2)
#--read_binary_term1=${PWD}/../outputs/writetest/eq32_term1.bin \
#--heldout_corpus_prefix=${HELDOUT} \
#--initial_lda_ss=${PWD}/../acl_arf-test_data/initial-lda-ss.dat \

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

echo "All Done"
