# load modules
module load julia cuda

# Set environments
export PROJECT_ROOT=/oscar/scratch/pthakur8/sem/seas-seme
export AMGX_DIR=$PROJECT_ROOT/deps/amgx
export JULIA_AMGX_PATH=$AMGX_DIR/build/libamgxsh.so

# Set LD LIBRARY PATH
export LD_LIBRARY_PATH=$AMGX_DIR/build:$LD_LIBRARY_PATH

echo "AMGX Environment Configured!"
