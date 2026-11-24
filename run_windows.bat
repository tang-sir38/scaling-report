@echo off
setlocal
cd /d "%~dp0"
set "HERE=%~dp0"
set "PYEXE="
set "APPVER=1.1.1"

rem ---- find a Python that already has numpy -----------------------------
call :try "%LOCALAPPDATA%\Programs\Python\Python313\python.exe"
call :try "%LOCALAPPDATA%\Programs\Python\Python312\python.exe"
call :try "%LOCALAPPDATA%\Programs\Python\Python311\python.exe"
call :try "%LOCALAPPDATA%\Programs\Python\Python310\python.exe"
call :try "%LOCALAPPDATA%\Programs\Python\Python39\python.exe"
call :try "C:/Python313/python.exe"
call :try "C:/Python312/python.exe"
call :try "C:/Python311/python.exe"
call :try "C:/ProgramData/anaconda3/python.exe"
call :try "C:/ProgramData/miniconda3/python.exe"
if defined PYEXE goto :routes
call :trycmd python
if defined PYEXE goto :routes
call :trycmd python3
if defined PYEXE goto :routes
goto :nopython

:try
if defined PYEXE exit /b 0
if not exist "%~1" exit /b 0
"%~1" -c "import numpy" >nul 2>&1
if not errorlevel 1 set "PYEXE=%~1"
exit /b 0

:trycmd
if defined PYEXE exit /b 0
where %~1 >nul 2>&1
if errorlevel 1 exit /b 0
%~1 -c "import numpy" >nul 2>&1
if not errorlevel 1 set "PYEXE=%~1"
exit /b 0

:nopython
echo.
echo   No Python with numpy was found on this machine.
echo.
echo   Install Python 3.9 or newer from
echo       https://www.python.org/downloads/
echo   and tick "Add Python to PATH" during setup, then run this again.
echo.
echo   If Python is installed already, open a command prompt and run
echo       python -m pip install numpy
echo.
pause
exit /b 1

rem ---- routes: command-line argument first, then the menu ----------------
:routes
set "ARG=%~1"
if /i "%ARG%"=="selftest" goto :selftest
if /i "%ARG%"=="demo" goto :demo
if "%ARG%"=="" goto :menu
set "TARGET=%~1"
goto :own
goto :menu

:menu
cls
echo.
echo   ==========================================================
echo    scaling_report  %APPVER%
echo    What a Poisson test really costs at your bin size
echo   ==========================================================
echo.
echo    Python: %PYEXE%
echo.
echo    [1]  Self-test        can this machine run it?      ~10 s
echo    [2]  Demo             fake data, full pipeline       ~1 min
echo    [3]  My own data      point at a Space Ranger folder
echo    [0]  Exit
echo.
choice /c 1230 /n /m "  Press 1, 2, 3 or 0: "
set "PICK=%errorlevel%"
echo.
if "%PICK%"=="1" goto :selftest
if "%PICK%"=="2" goto :demo
if "%PICK%"=="3" goto :own
if "%PICK%"=="4" goto :bye
goto :menu

rem ---- 1 self-test -----------------------------------------------------
:selftest
echo   Running the self-test - seven checks, no data needed.
echo.
"%PYEXE%" -m scaling_report --selftest
if errorlevel 1 goto :selftest_bad
echo.
echo   All checks passed. This machine is ready.
goto :hold
:selftest_bad
echo.
echo   The self-test did NOT pass. Copy the lines above when asking for help.
goto :hold

rem ---- 2 demo ----------------------------------------------------------
:demo
if not exist "%HERE%demo" mkdir "%HERE%demo"
echo   Step 1 of 2 - making synthetic Visium HD output in demo
echo.
"%PYEXE%" -m scaling_report --demo "%HERE%demo" --ds 2,4,8,16,32,64
if errorlevel 1 goto :failed
echo.
echo   Step 2 of 2 - reading it back and writing the report
echo.
"%PYEXE%" -m scaling_report --from-10x "%HERE%demo" --out "%HERE%report" --ds 2,4,8,16,32,64 --min-mu 0.05
if errorlevel 1 goto :failed
goto :done

rem ---- 3 my own data ---------------------------------------------------
:own
if defined TARGET goto :own_go
echo   Drag your Space Ranger folder into this window, then press Enter.
echo   It can be the run folder, its outs folder, or binned_outputs.
echo.
set /p TARGET=  Folder: 
if not defined TARGET goto :own_empty
:own_go
set "TARGET=%TARGET:"=%"
if not exist "%TARGET%" goto :own_missing
echo.
echo   Reading: %TARGET%
echo   A 2 um matrix can take a while; the larger bins are fast.
echo.
"%PYEXE%" -m scaling_report --from-10x "%TARGET%" --out "%HERE%report"
if errorlevel 1 goto :failed
goto :done

:own_empty
echo.
echo   Nothing given, going back to the menu.
goto :hold
:own_missing
echo.
echo   That path does not exist:
echo       %TARGET%
goto :hold

rem ---- shared endings --------------------------------------------------
:done
echo.
if not exist "%HERE%report\scaling_report.html" goto :done_nohtml
echo   Done. Opening report\scaling_report.html
start "" "%HERE%report\scaling_report.html"
goto :done_tail
:done_nohtml
echo   Done, but no HTML was written. Look in the report folder.
:done_tail
echo.
echo   Everything is in  %HERE%report
goto :hold

:failed
echo.
echo   The run stopped with an error. The message above says why.
echo   If it is about numpy, run:  "%PYEXE%" -m pip install numpy
goto :hold

:hold
echo.
pause
goto :menu

:bye
endlocal
exit /b 0
