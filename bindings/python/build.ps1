# Build the cpp-math Python extension using pybind11 + MinGW g++.
$ErrorActionPreference = "Stop"

$src = "bindings\python\module.cpp"
$cppSources = (Get-ChildItem "src\*.cpp" | ForEach-Object { $_.FullName }) -join " "

$pyIncludes = python -c "import pybind11; print(pybind11.get_include())"
$pythonIncludes = python -c "import sysconfig; print('-I' + sysconfig.get_paths()['include'])"
$pythonLib = python -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))"
$ext = python -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
$outFile = "mathx$ext"

# Python import library (pythonXY.lib). On Windows this must be linked so
# the extension resolves symbols from the CPython DLL.
$pyLib = Get-ChildItem $pythonLib -Filter "python*.lib" | Select-Object -First 1
if (-not $pyLib) { throw "No python*.lib found in $pythonLib" }

$cmd = "g++ -O3 -march=native -fopenmp -std=c++17 -DMATHX_BUILDING_LIBRARY " +
       "-Iinclude -I`"$pyIncludes`" $pythonIncludes -shared -fPIC " +
       "$src $cppSources -o `"$outFile`" -L`"$pythonLib`" -l$([IO.Path]::GetFileNameWithoutExtension($pyLib.Name)) -O3"

Write-Host "Building $outFile ..."
Invoke-Expression $cmd
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

# Test import. On Windows the extension depends on the MinGW runtime DLLs
# (libgomp, libstdc++, ...); add their directory so CPython can resolve them.
$mingwBin = "C:\msys64\ucrt64\bin"
if (Test-Path $mingwBin) {
    python -c "import os; os.add_dll_directory(r'$mingwBin'); import mathx; print('OK:', mathx.lambert_w(1.0), mathx.__version__)"
} else {
    python -c "import mathx; print('OK:', mathx.lambert_w(1.0), mathx.__version__)"
}
