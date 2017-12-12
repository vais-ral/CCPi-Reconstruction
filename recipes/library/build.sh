if [ -z "$CIL_VERSION" ]; then
    echo "Need to set CIL_VERSION"
    exit 1
fi  





mkdir ${SRC_DIR}/build
cp -r "${RECIPE_DIR}/../../Core" ${SRC_DIR}/build
cd ${SRC_DIR}/build

cmake -G "Unix Makefiles" ${SRC_DIR}/build