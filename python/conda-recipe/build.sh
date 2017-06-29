mkdir ${SRC_DIR}/ccpi
cp -r ${RECIPE_DIR}/../../ ${SRC_DIR}/ccpi

cd ${SRC_DIR}/ccpi/python
$PYTHON setup.py install
