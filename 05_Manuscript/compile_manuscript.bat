@echo off
REM Batch script to compile LaTeX manuscript on Windows
REM Date: 2025-12-09
REM Requires: MiKTeX or TeX Live installed

echo ====================================
echo Compiling AML Manuscript (LaTeX)
echo ====================================
echo.

cd /d "%~dp0"

REM First pass - generate aux files
echo [1/4] Running pdflatex (first pass)...
pdflatex -interaction=nonstopmode AML_Molecular_Subtypes_Manuscript.tex
if errorlevel 1 (
    echo ERROR: First pdflatex pass failed
    pause
    exit /b 1
)

REM Second pass - resolve references
echo.
echo [2/4] Running pdflatex (second pass)...
pdflatex -interaction=nonstopmode AML_Molecular_Subtypes_Manuscript.tex
if errorlevel 1 (
    echo ERROR: Second pdflatex pass failed
    pause
    exit /b 1
)

REM Third pass - final formatting
echo.
echo [3/4] Running pdflatex (third pass)...
pdflatex -interaction=nonstopmode AML_Molecular_Subtypes_Manuscript.tex
if errorlevel 1 (
    echo ERROR: Third pdflatex pass failed
    pause
    exit /b 1
)

REM Clean up auxiliary files
echo.
echo [4/4] Cleaning up auxiliary files...
del AML_Molecular_Subtypes_Manuscript.aux 2>nul
del AML_Molecular_Subtypes_Manuscript.log 2>nul
del AML_Molecular_Subtypes_Manuscript.out 2>nul
del AML_Molecular_Subtypes_Manuscript.toc 2>nul
del AML_Molecular_Subtypes_Manuscript.bbl 2>nul
del AML_Molecular_Subtypes_Manuscript.blg 2>nul

echo.
echo ====================================
echo SUCCESS: PDF compiled successfully!
echo ====================================
echo.
echo Output: AML_Molecular_Subtypes_Manuscript.pdf
echo.

REM Open the PDF
if exist AML_Molecular_Subtypes_Manuscript.pdf (
    echo Opening PDF...
    start AML_Molecular_Subtypes_Manuscript.pdf
)

pause
