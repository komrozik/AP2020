all: build/main.pdf

# hier Python-Skripte:
build/CU_Spektrum.pdf: plot.py matplotlibrc header-matplotlib.tex | build
	TEXINPUTS=$$(pwd): python plot.py
build/Transmission.pdf: plot.py matplotlibrc header-matplotlib.tex | build
	TEXINPUTS=$$(pwd): python plot.py


# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/main.pdf: build/CU_Spektrum.pdf build/Transmission.pdf

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
