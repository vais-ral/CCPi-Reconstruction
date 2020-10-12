echo $CIL_VERSION

mkdir ${SRC_DIR}/build
cp -r ${RECIPE_DIR}/../ ${SRC_DIR}/build
mkdir ${SRC_DIR}/build/build
cd ${SRC_DIR}/build/build
cmake -G "Unix Makefiles" -DLIBRARY_LIB="${CONDA_PREFIX}/lib" -DLIBRARY_INC="${CONDA_PREFIX}" -DCMAKE_INSTALL_PREFIX="${PREFIX}" -DCONDA_BUILD=ON ..
make install
