IF NOT DEFINED CIL_VERSION (
ECHO CIL_VERSION Not Defined.
exit 1
)

mkdir "%SRC_DIR%\build"
ROBOCOPY /E "%RECIPE_DIR%\..\..\Core" "%SRC_DIR%\build"
::ROBOCOPY /E "%RECIPE_DIR%\..\..\Wrappers\python\src" "%SRC_DIR%\build\module"
cd "%SRC_DIR%\build"

echo "we should be in %SRC_DIR%\build"

cmake -G "NMake Makefiles" -DBOOST_ROOT="%LIBRARY_PREFIX%" "%SRC_DIR%\build"

:: Build C library
nmake
if errorlevel 1 exit 1

:: Install step
