FC = gfortran
FCFLAGS = -g -ffpe-trap=zero,invalid,overflow,underflow
LDFLAGS = -llapack -lblas

.PHONY: all clean init

all: init rq gera

init:
	@mkdir -p build
	@cp ajuste_polinomial.py build/
	@cp run.sh build/
	@echo "Diretório build criado e arquivo copiado"

rq:
	@echo "Compilando residuos_quadraticos..."
	@$(FC) -o build/rq src/rq_linear.f90 $(FCFLAGS) $(LDFLAGS)
	@$(FC) -o build/rq_nao_linear src/rq_nao_linear.f90 $(FCFLAGS) $(LDFLAGS)

gera:
	@echo "Compilando gera_input..."
	@$(FC) $(FCFLAGS) -o build/gera src/gera_input_linear.f90
	@$(FC) $(FCFLAGS) -o build/gera src/gera_input_nao_linear.f90

clean:
	@echo "Removendo arquivos gerados..."
	@$(RM) -r build


