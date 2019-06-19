if [ -z "$CIL_VERSION" ]; then
    echo "Need to set CIL_VERSION"
    exit 1
fi  
mkdir ${SRC_DIR}/ccpi
cp -r "${RECIPE_DIR}/../../../" ${SRC_DIR}/ccpi

cd ${SRC_DIR}/ccpi/Wrappers/Python
cmake .

echo "This is Path to the Python executable in the host prefix ${PYTHON}"
${PYTHON} setup.py install
