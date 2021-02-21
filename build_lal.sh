GWSYNTH_DIR=$PWD
DEPS_DIR=$GWSYNTH_DIR/deps
NPROC=`nproc`
export PKG_CONFIG_PATH="$HOME/.conda/envs/local/lib/pkgconfig:$PKG_CONFIG_PATH"

# Build lal
cd $DEPS_DIR
tar -xJf lal-7.1.2.tar.xz
cd lal-7.1.2
./configure --prefix=$GWSYNTH_DIR/lal --disable-swig
make -j $NPROC
make install
export PKG_CONFIG_PATH="$GWSYNTH_DIR/lal/lib/pkgconfig:$PKG_CONFIG_PATH"

# Build lalsimulation
cd $DEPS_DIR
tar -xJf lalsimulation-2.5.1.tar.xz
cd lalsimulation-2.5.1
./configure --prefix=$GWSYNTH_DIR/lal --disable-swig
make -j $NPROC
make install

export LIBRARY_PATH="$GWSYNTH_DIR/lal/lib:LIBRARY_PATH"
export LD_LIBRARY_PATH="$GWSYNTH_DIR/lal/lib:$LD_LIBRARY_PATH"

echo 'export LIBRARY_PATH="$GWSYNTH_DIR/lal/lib:LIBRARY_PATH'
echo 'export LD_LIBRARY_PATH="$GWSYNTH_DIR/lal/lib:$LD_LIBRARY_PATH'