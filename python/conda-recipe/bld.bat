xcopy /e "%RECIPE_DIR%\.." "%SRC_DIR%\ccpi"
xcopy /e "%RECIPE_DIR%\..\..\src" "%SRC_DIR%\ccpi"
xcopy /e "%RECIPE_DIR%\..\..\src\Algorithms" "%SRC_DIR%\ccpi"
xcopy /e "%RECIPE_DIR%\..\..\src\Readers" "%SRC_DIR%\ccpi"
cd ccpi
set BOOST_INCLUDE="C:\Apps\Anaconda3\envs\cgls\Library\include"
%PYTHON% setup.py build_ext
if errorlevel 1 exit 1
%PYTHON% setup.py install
if errorlevel 1 exit 1