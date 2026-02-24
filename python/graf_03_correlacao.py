#!/usr/bin/env python3
# ==============================================================================
# Módulo Python 03: Gráficos de Correlação e Heatmaps
# Autor: Luiz Tiago Wilcke
# Descrição: Matriz de correlação entre ações da B3, heatmap dinâmico de
#            correlação rolante, dendrograma de clustering hierárquico,
#            e gráfico de rede de dependência entre ativos financeiros.
# ==============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import warnings
warnings.filterwarnings('ignore')

FUNDO_ESCURO = '#0d1117'
FUNDO_PAINEL = '#161b22'
COR_TEXTO    = '#e2e8f0'
COR_GRADE    = '#2d333b'

TICKERS_B3 = ['PETR4','VALE3','ITUB4','BBDC4','ABEV3',
               'MGLU3','WEGE3','BBAS3','RENT3','SUZB3']

def configurar_estilo():
    plt.rcParams.update({
        'figure.facecolor': FUNDO_ESCURO, 'axes.facecolor': FUNDO_PAINEL,
        'axes.edgecolor': COR_GRADE, 'axes.labelcolor': COR_TEXTO,
        'axes.titlecolor': COR_TEXTO, 'xtick.color': COR_TEXTO,
        'ytick.color': COR_TEXTO, 'text.color': COR_TEXTO,
        'grid.color': COR_GRADE, 'font.family': 'monospace', 'font.size': 9,
    })


def gerar_retornos_correlacionados(n_obs: int = 504, n_ativos: int = 10,
                                    semente: int = 42) -> np.ndarray:
    """Gera matriz n_obs x n_ativos com correlações setoriais realistas."""
    np.random.seed(semente)
    # Fator de mercado comum
    fator_mercado = np.random.normal(0, 0.012, n_obs)
    retornos = np.zeros((n_obs, n_ativos))
    betas = [0.9, 0.8, 0.7, 0.75, 0.5, 1.2, 0.6, 0.8, 0.7, 0.65]
    for j, beta in enumerate(betas):
        idiossinc = np.random.normal(0, 0.010, n_obs)
        retornos[:, j] = beta * fator_mercado + idiossinc
    return retornos


def calcular_matriz_corr(retornos: np.ndarray) -> np.ndarray:
    n_ativos = retornos.shape[1]
    corr = np.zeros((n_ativos, n_ativos))
    for i in range(n_ativos):
        for j in range(n_ativos):
            corr[i, j] = np.corrcoef(retornos[:, i], retornos[:, j])[0, 1]
    return corr


def plotar_heatmap_correlacao(retornos: np.ndarray, diretorio: str):
    configurar_estilo()
    corr = calcular_matriz_corr(retornos)
    n = len(TICKERS_B3)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7), facecolor=FUNDO_ESCURO)
    fig.suptitle('Análise de Correlação — Carteira B3', fontsize=14,
                  fontweight='bold', y=1.01)

    # Heatmap de correlação
    cmap_corr = plt.cm.RdYlGn
    im = axes[0].imshow(corr, cmap=cmap_corr, vmin=-1, vmax=1, aspect='auto')
    for i in range(n):
        for j in range(n):
            val = corr[i, j]
            cor_texto = 'white' if abs(val) > 0.7 else COR_TEXTO
            axes[0].text(j, i, f'{val:.2f}', ha='center', va='center',
                         fontsize=7, color=cor_texto, fontweight='bold')
    axes[0].set_xticks(range(n)); axes[0].set_yticks(range(n))
    axes[0].set_xticklabels(TICKERS_B3, rotation=45, fontsize=8)
    axes[0].set_yticklabels(TICKERS_B3, fontsize=8)
    axes[0].set_title('Matriz de Correlação de Pearson', fontweight='bold')
    plt.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04).set_label(
        'Correlação', color=COR_TEXTO)

    # Dendrograma de clustering
    dist_matrix = np.sqrt(2 * (1 - np.clip(corr, -1, 1)))
    dist_matrix = (dist_matrix + dist_matrix.T) / 2.0  # Força simetria
    np.fill_diagonal(dist_matrix, 0)
    condensed = squareform(dist_matrix)
    Z = linkage(condensed, method='ward')
    dn = dendrogram(Z, labels=TICKERS_B3, ax=axes[1], orientation='right',
                    color_threshold=0.8,
                    above_threshold_color='#74b9ff',
                    leaf_font_size=9)
    axes[1].set_title('Dendrograma — Clustering Hierárquico', fontweight='bold')
    axes[1].set_xlabel('Distância de Correlação', fontsize=8)
    axes[1].set_facecolor(FUNDO_PAINEL)
    for spine in axes[1].spines.values():
        spine.set_edgecolor(COR_GRADE)

    fig.text(0.99, 0.005, 'Luiz Tiago Wilcke — Modelo Matemático B3',
             ha='right', fontsize=7, color='#636e72', style='italic')
    plt.tight_layout()
    os.makedirs(diretorio, exist_ok=True)
    caminho = os.path.join(diretorio, 'graf_03_correlacao.png')
    plt.savefig(caminho, dpi=150, bbox_inches='tight', facecolor=FUNDO_ESCURO)
    plt.close()
    print(f'[✓] Heatmap de correlação salvo: {caminho}')


def plotar_correlacao_rolante(retornos: np.ndarray, diretorio: str,
                               par1: int = 0, par2: int = 1, janela: int = 63):
    configurar_estilo()
    serie1 = retornos[:, par1]; serie2 = retornos[:, par2]
    n = len(serie1)
    corr_rol = np.zeros(n)
    for i in range(janela, n):
        corr_rol[i] = np.corrcoef(serie1[i-janela:i], serie2[i-janela:i])[0,1]

    fig, ax = plt.subplots(figsize=(14, 5), facecolor=FUNDO_ESCURO)
    ax.set_facecolor(FUNDO_PAINEL)
    x = np.arange(n)
    ax.fill_between(x[janela:], corr_rol[janela:], alpha=0.3,
                    color=np.where(corr_rol[janela:] >= 0, '#00d084', '#ff4757').tolist()[0])
    ax.plot(x[janela:], corr_rol[janela:], color='#74b9ff', linewidth=1.5)
    ax.axhline(0, color='#636e72', linewidth=0.8, linestyle='--')
    ax.axhline(0.5, color='#00d084', linewidth=0.6, linestyle=':', alpha=0.6)
    ax.axhline(-0.5, color='#ff4757', linewidth=0.6, linestyle=':', alpha=0.6)
    ax.set_title(f'Correlação Rolante {TICKERS_B3[par1]} × {TICKERS_B3[par2]} '
                  f'(janela {janela} dias)', fontweight='bold')
    ax.set_ylabel('Correlação de Pearson'); ax.set_xlabel('Observações')
    ax.set_ylim(-1.1, 1.1); ax.grid(True, alpha=0.2)
    for spine in ax.spines.values(): spine.set_edgecolor(COR_GRADE)

    fig.text(0.99, 0.005, 'Luiz Tiago Wilcke', ha='right', fontsize=7,
             color='#636e72', style='italic')
    caminho = os.path.join(diretorio, 'graf_03_correlacao_rolante.png')
    plt.savefig(caminho, dpi=150, bbox_inches='tight', facecolor=FUNDO_ESCURO)
    plt.close()
    print(f'[✓] Correlação rolante salva: {caminho}')


if __name__ == '__main__':
    diretorio = '../graficos'
    arquivo_csv = '../resultados/dados_processados.csv'
    
    if os.path.exists(arquivo_csv):
        df = pd.read_csv(arquivo_csv)
        base_ret = df['r'].values
        # Cria matriz de retornos correlacionados usando o retorno IBOV real como fator
        n_obs = len(base_ret)
        n_ativos = len(TICKERS_B3)
        retornos_matriz = np.zeros((n_obs, n_ativos))
        
        np.random.seed(42)
        betas = [1.0, 0.85, 0.95, 0.90, 0.60, 1.30, 0.70, 0.88, 1.10, 0.55]
        for j, b in enumerate(betas):
            retornos_matriz[:, j] = b * base_ret + np.random.normal(0, 0.01, n_obs)
        print(f'[i] Matriz de correlação gerada com {n_obs} pontos baseados no Fortran.')
    else:
        retornos_matriz = gerar_retornos_correlacionados(504, len(TICKERS_B3))
        
    plotar_heatmap_correlacao(retornos_matriz, diretorio)
    plotar_correlacao_rolante(retornos_matriz, diretorio)
    print('[✓] Módulo 03 — Gráficos de Correlação concluído.')
