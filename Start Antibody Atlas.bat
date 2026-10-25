@echo off
REM Double-click this file to start Antibody Atlas on Windows.
REM It starts the local server and opens your browser. Leave this window open
REM while you use the app; close it or press Ctrl-C to stop.
title Antibody Atlas
cd /d "%~dp0"
echo Antibody Atlas
echo Folder: %CD%
echo.
where powershell >nul 2>nul
if errorlevel 1 (
  echo ERROR: PowerShell was not found, so run_local.ps1 cannot be started.
  pause
  exit /b 1
)
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_local.ps1"
echo.
echo Antibody Atlas has stopped.
pause
