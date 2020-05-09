all: build/main.pdf

# hier Python-Skripte:
build/Bragg.pdf: plot.py matplotlibrc header-matplotlib.tex | build
	TEXINPUTS=$$(pwd): python plot.py
build/Cu_Emission.pdf: plot.py matplotlibrc header-matplotlib.tex | build
	TEXINPUTS=$$(pwd): python plot.py
build/Brom.pdf: plot.py matplotlibrc header-matplotlib.tex | build
	TEXINPUTS=$$(pwd): python plot.py

# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/main.pdf: build/Bragg.pdf build/Brom.pdf build/Cu_Emission.pdf

build/main.pdf: FORCE | build
	  TEXINPUTS=build: \
	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk -pvc\
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	main.tex

build:
	mkdir -p build
	mkdir -p plots

clean:
	rm -rf build
	rm -rf plots

FORCE:

.PHONY: all clean
