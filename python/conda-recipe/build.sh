cp -r ${RECIPE_DIR}/.. ${SRC_DIR}/ccpi
cp -r ${RECIPE_DIR}/../../src/* ${SRC_DIR}/ccpi
cp ${RECIPE_DIR}/../../src/Algorithms/* ${SRC_DIR}/ccpi
cp ${RECIPE_DIR}/../../src/Readers/* ${SRC_DIR}/ccpi

cd ${SRC_DIR}/ccpi
$PYTHON setup.py install
