@echo off
rem Starts the Polya Pipeline window. The console stays open only if something goes wrong.
cd /d "%~dp0"
where python >nul 2>nul
if errorlevel 1 (
    echo Python was not found. Install it from https://www.python.org/downloads/
    echo and tick "Add python.exe to PATH" during the installation.
    pause
    exit /b 1
)
python "%~dp0derivatives_gui.py"
if errorlevel 1 (
    echo.
    echo Polya Pipeline stopped with an error. If a package is missing, run:
    echo     python -m pip install -r requirements.txt
    pause
)
