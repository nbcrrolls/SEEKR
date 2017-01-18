WORKDIR=$WORK
HOMEDIR=$HOME

uname -a

module load intel

module load fftw3

echo $PATH

cat STDERR_REDIRECT_WORKING > /dev/null

echo $WORKDIR/charm-6.7.0-verbs-linux-x86_64-iccstatic

mkdir $WORKDIR/charm-6.7.0-verbs-linux-x86_64-iccstatic

cd $WORKDIR/charm-6.7.0-verbs-linux-x86_64-iccstatic

wget --no-verbose http://www.ks.uiuc.edu/~jim/build/charm-6.7.0.tar

tar xmf $WORKDIR/charm-6.7.0-verbs-linux-x86_64-iccstatic/charm-6.7.0.tar

/bin/rm -f $WORKDIR/charm-6.7.0-verbs-linux-x86_64-iccstatic/charm-6.7.0.tar

cd $WORKDIR/charm-6.7.0-verbs-linux-x86_64-iccstatic/charm-6.7.0

touch src/util/pup_f.f90



./build charm++ verbs-linux-x86_64 iccstatic --no-build-shared --with-production 

cd verbs-linux-x86_64-iccstatic/tests/charm++/megatest

make pgm

ldd pgm

mv pgm ../../../bin/megatest

cd ../../converse/megacon

make pgm

ldd pgm

mv pgm ../../../bin/megacon

cd ../../..

/bin/rm -rf $HOMEDIR/charm-6.7.0/verbs-linux-x86_64-iccstatic

mkdir $HOMEDIR/charm-6.7.0

chmod a+rX $HOMEDIR/charm-6.7.0

mkdir $HOMEDIR/charm-6.7.0/verbs-linux-x86_64-iccstatic

cp -rL bin include $HOMEDIR/charm-6.7.0/verbs-linux-x86_64-iccstatic

cd lib

mkdir $HOMEDIR/charm-6.7.0/verbs-linux-x86_64-iccstatic/lib

mv lib*modulemsa.a lib*moduleNeighborLB.a lib*moduleHybridLB.a lib*moduleRefineLB.a lib*moduleGreedyLB.a lib*moduleCkMulticast.a lib*moduleCkIO.a lib*moduleCkLoop.a lib*ck.a lib*ckmain.a lib*conv-cplus-y.a lib*conv-core.a lib*conv-util.a lib*conv-partition.a lib*conv-ldb.a cray_tlbh*.* libmemory-default.* libmemory-os.* libthreads-default.* lib*ckqt.a lib*trace-projections.a lib*trace-summary.a libldb-rand.* ~/charm-6.7.0/verbs-linux-x86_64-iccstatic/lib

mkdir $HOMEDIR/charm-6.7.0/verbs-linux-x86_64-iccstatic/lib_so

cd; /bin/rm -rf $WORKDIR/charm-6.7.0-verbs-linux-x86_64-iccstatic

chmod -R a+rX $HOMEDIR/charm-6.7.0/verbs-linux-x86_64-iccstatic

ls -l $HOMEDIR/charm-6.7.0/verbs-linux-x86_64-iccstatic/lib

ls -l $HOMEDIR/charm-6.7.0/verbs-linux-x86_64-iccstatic/bin/charmrun $HOMEDIR/charm-6.7.0/verbs-linux-x86_64-iccstatic/bin/megatest

cat STDERR_REDIRECT_WORKING > /dev/null

echo Done building via scp.
