# Modelo Matemático da Bolsa de Valores (B3)

Este projeto consiste em um motor matemático de alto desempenho desenvolvido em **Fortran 90/95/2003** e uma interface de visualização em **Python**. O objetivo é simular e analisar o comportamento de ativos da Bolsa de Valores brasileira (B3) utilizando métodos quantitativos avançados.

## 🏗️ Estrutura do Projeto

O sistema é dividido em 38 módulos Fortran que cobrem:
- Equações Diferenciais (Estocásticas e Ordinárias)
- Modelos de Volatilidade (GARCH, EWMA)
- Análise Técnica (RSI, MACD, Bandas de Bollinger)
- Filtro de Kalman e Processamento de Sinais (Fourier, Ondaletas)
- Análise Fractal e Caos (Hurst, DFA)
- Precificação de Derivativos (Black-Scholes, Heston)
- Otimização de Portfólio (Markowitz, Algoritmos Genéticos)

## 🛠️ Como Executar

### 1. Compilação do Motor Fortran
Requer `gfortran` instalado.

```bash
cd fortran
gfortran -O3 -c *.f90
gfortran *.o -o ../bolsa_valores
cd ..
```

### 2. Execução da Simulação
```bash
./bolsa_valores
```
Isso gerará o arquivo `resultados/dados_processados.csv`.

### 3. Geração de Gráficos (Python)
Requer `pandas`, `numpy` e `matplotlib`.

```bash
cd python
python3 graf_01_precos.py
python3 graf_02_volatilidade.py
python3 graf_03_correlacao.py
```

## 📊 Visualizações
Os gráficos gerados estarão disponíveis na pasta `graficos/`.

---
**Autor:** Luiz Tiago Wilcke
**Tecnologias:** Fortran 2003, Python 3.x, Git
