all: build/main.pdf

# hier Python-Skripte:
plots/Brom.pdf plots/Cu_Peaks.pdf plots/Gallium.pdf plots/Rubidium.pdf plots/Strontium.pdf plots/Zink.pdf plots/Zirkonium.pdf plots/Cu_Emission.pdf plots/Bragg.pdf: plot.py matplotlibrc header-matplotlib.tex | build
	TEXINPUTS=$$(pwd): python plot.py

# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/main.pdf: plots/Bragg.pdf plots/Brom.pdf plots/Cu_Peaks.pdf plots/Gallium.pdf plots/Rubidium.pdf plots/Strontium.pdf plots/Zink.pdf plots/Zirkonium.pdf plots/Cu_Emission.pdf

build/main.pdf: FORCE | build
	  TEXINPUTS=build: \
	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk -pvc\
	  --lualatex \
	  --output-directory=build \
	  --interaction=errorstopmode \
	main.tex

build:
	mkdir -p build
	mkdir -p plots

clean:
	rm -rf build
	rm -rf plots

FORCE:

.PHONY: all clean
