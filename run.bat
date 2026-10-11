@echo off
cd /d "%~dp0"
title Flip Vacuum (POSCAR / CIF)

where python >nul 2>nul
if errorlevel 1 (
    echo [ERROR] python not found in PATH. Please install Python first.
    echo.
    pause
    exit /b 1
)

echo Starting GUI...
python flip_vacuum.py
if errorlevel 1 (
    echo.
    echo [ERROR] Program exited with an error. See message above.
    echo.
    pause
)
