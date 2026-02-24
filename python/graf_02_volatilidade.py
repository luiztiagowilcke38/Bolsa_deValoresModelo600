#!/usr/bin/env python3
# ==============================================================================
# Módulo Python 02: Gráficos de Volatilidade
# Autor: Luiz Tiago Wilcke
# Descrição: Visualização de volatilidade histórica, GARCH condicional,
#            superfície de volatilidade implícita (smile/skew de opções)
#            e comparativo de modelos para ativos da B3.
# ==============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D
import warnings
warnings.filterwarnings('ignore')

FUNDO_ESCURO = '#0d1117'
FUNDO_PAINEL = '#161b22'
COR_TEXTO    = '#e2e8f0'
COR_GRADE    = '#2d333b'
DIAS_ANO     = 252


def configurar_estilo():
    plt.rcParams.update({
        'figure.facecolor': FUNDO_ESCURO,
        'axes.facecolor':   FUNDO_PAINEL,
        'axes.edgecolor':   COR_GRADE,
        'axes.labelcolor':  COR_TEXTO,
        'axes.titlecolor':  COR_TEXTO,
        'xtick.color':      COR_TEXTO,
        'ytick.color':      COR_TEXTO,
        'text.color':       COR_TEXTO,
        'grid.color':       COR_GRADE,
        'grid.alpha':       0.4,
        'font.family':      'monospace',
        'font.size':        9,
    })


def calcular_vol_historica(retornos: np.ndarray, janela: int = 21) -> np.ndarray:
    """Volatilidade histórica janela deslizante, anualizada."""
    n   = len(retornos)
    vol = np.zeros(n)
    for i in range(janela, n):
        vol[i] = np.std(retornos[i-janela:i], ddof=1) * np.sqrt(DIAS_ANO)
    return vol


def simular_garch11(retornos: np.ndarray, omega: float, alpha: float,
                     beta: float) -> np.ndarray:
    """Simula variância condicional GARCH(1,1)."""
    n = len(retornos)
    sigma2 = np.zeros(n)
    sigma2[0] = np.var(retornos)
    for t in range(1, n):
        sigma2[t] = omega + alpha * retornos[t-1]**2 + beta * sigma2[t-1]
    return np.sqrt(sigma2 * DIAS_ANO)   # vol. anualizada


def plotar_volatilidade_comparativa(retornos: np.ndarray, diretorio: str):
    """
    Painel comparando:
    1) Retornos diários (|r_t|)
    2) Vol. histórica (21 dias) vs Vol. GARCH(1,1)
    3) Cluster de volatilidade (autocorrelação de |r_t|^2)
    4) Distribuição dos retornos vs normal
    """
    configurar_estilo()
    n = len(retornos)
    vol_hist  = calcular_vol_historica(retornos, 21)
    vol_garch = simular_garch11(retornos, 5e-6, 0.08, 0.90)
    ewma_lambda = 0.94
    vol_ewma = np.zeros(n)
    vol_ewma[0] = np.var(retornos)
    for t in range(1, n):
        vol_ewma[t] = ewma_lambda * vol_ewma[t-1] + (1-ewma_lambda) * retornos[t-1]**2
    vol_ewma = np.sqrt(vol_ewma * DIAS_ANO)

    fig = plt.figure(figsize=(16, 11), facecolor=FUNDO_ESCURO)
    fig.suptitle('Análise de Volatilidade — Modelo Matemático B3',
                  fontsize=14, fontweight='bold', y=0.98)
    gs = GridSpec(3, 2, figure=fig, hspace=0.4, wspace=0.3)

    # 1) Retornos diários
    ax1 = fig.add_subplot(gs[0, :])
    cores = np.where(retornos >= 0, '#00d084', '#ff4757')
    ax1.bar(range(n), retornos * 100, color=cores, alpha=0.7, linewidth=0)
    ax1.set_title('Retornos Logarítmicos Diários (%)', fontweight='bold')
    ax1.set_ylabel('Retorno (%)')
    ax1.axhline(0, color='#636e72', linewidth=0.5)
    ax1.grid(True, alpha=0.2)

    # 2) Comparativo de volatilidade
    ax2 = fig.add_subplot(gs[1, :])
    ax2.plot(vol_hist,  color='#74b9ff', linewidth=1.2, label='Vol. Histórica (21d)',  alpha=0.9)
    ax2.plot(vol_garch, color='#fd79a8', linewidth=1.2, label='GARCH(1,1)',             alpha=0.9)
    ax2.plot(vol_ewma,  color='#ffeaa7', linewidth=1.0, label='EWMA (λ=0.94)',          alpha=0.8)
    ax2.fill_between(range(n), vol_hist, alpha=0.05, color='#74b9ff')
    ax2.set_title('Volatilidade Condicional Anualizada (%)', fontweight='bold')
    ax2.set_ylabel('Vol. Anualizada (%)')
    ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x*100:.0f}%'))
    ax2.legend(loc='upper right', fontsize=8)
    ax2.grid(True, alpha=0.2)

    # 3) Autocorrelação de r_t^2 (cluster de volatilidade)
    ax3 = fig.add_subplot(gs[2, 0])
    ret2 = retornos**2
    n_lags = 30
    media_r2 = np.mean(ret2)
    var_r2   = np.var(ret2)
    acf_vals = np.array([
        np.mean((ret2[:n-k] - media_r2)*(ret2[k:] - media_r2)) / (var_r2 + 1e-15)
        for k in range(1, n_lags+1)
    ])
    intervalo = 1.96 / np.sqrt(n)
    ax3.bar(range(1, n_lags+1), acf_vals, color='#a29bfe', alpha=0.8, linewidth=0)
    ax3.axhline( intervalo, color='#ff7675', linewidth=1, linestyle='--', alpha=0.6, label='IC 95%')
    ax3.axhline(-intervalo, color='#ff7675', linewidth=1, linestyle='--', alpha=0.6)
    ax3.set_title('ACF de r²(t) — Cluster de Volatilidade', fontweight='bold')
    ax3.set_xlabel('Lag (dias)')
    ax3.set_ylabel('ACF')
    ax3.legend(fontsize=8); ax3.grid(True, alpha=0.2)

    # 4) Distribuição vs Normal
    ax4 = fig.add_subplot(gs[2, 1])
    mu_r, sigma_r = np.mean(retornos), np.std(retornos)
    ax4.hist(retornos, bins=60, density=True, color='#00b894', alpha=0.6,
             label='Retornos observados')
    x_lin = np.linspace(mu_r - 4*sigma_r, mu_r + 4*sigma_r, 200)
    normal_pdf = np.exp(-0.5*((x_lin - mu_r)/sigma_r)**2) / (sigma_r * np.sqrt(2*np.pi))
    ax4.plot(x_lin, normal_pdf, color='#fd79a8', linewidth=2, label='Distribuição Normal')
    skew = np.mean(((retornos - mu_r)/sigma_r)**3)
    kurt = np.mean(((retornos - mu_r)/sigma_r)**4) - 3
    ax4.set_title(f'Distribuição dos Retornos\nSkew={skew:.3f}  Kurt={kurt:.3f}',
                  fontweight='bold')
    ax4.set_xlabel('Retorno diário')
    ax4.legend(fontsize=8); ax4.grid(True, alpha=0.2)

    fig.text(0.99, 0.005, 'Luiz Tiago Wilcke — Modelo Matemático B3',
             ha='right', fontsize=7, color='#636e72', style='italic')
    os.makedirs(diretorio, exist_ok=True)
    caminho = os.path.join(diretorio, 'graf_02_volatilidade.png')
    plt.savefig(caminho, dpi=150, bbox_inches='tight', facecolor=FUNDO_ESCURO)
    plt.close()
    print(f'[✓] Gráfico de volatilidade salvo: {caminho}')


def plotar_superficie_volatilidade(diretorio: str):
    """
    Superfície de volatilidade implícita (smile/skew)
    sigma_impl(K/S, T) — gerada com dados sintéticos realistas.
    """
    configurar_estilo()
    moneyness = np.linspace(0.7, 1.3, 25)
    maturidades = np.linspace(1/12, 2.0, 20)  # 1 mês a 2 anos
    M, T = np.meshgrid(moneyness, maturidades)

    # Modelo simplificado de smile + term structure
    vol_atm  = 0.22
    skew_coef = -0.15
    smile_coef = 0.10
    term_decay = 0.30

    vol_surface = (vol_atm
                   + skew_coef * (M - 1.0)
                   + smile_coef * (M - 1.0)**2
                   + term_decay / np.sqrt(T + 0.01) * 0.02)
    vol_surface = np.clip(vol_surface, 0.05, 0.70)

    fig = plt.figure(figsize=(14, 8), facecolor=FUNDO_ESCURO)
    ax  = fig.add_subplot(111, projection='3d')
    ax.set_facecolor(FUNDO_PAINEL)
    fig.patch.set_facecolor(FUNDO_ESCURO)

    surf = ax.plot_surface(M, T, vol_surface * 100,
                           cmap='plasma', alpha=0.85,
                           linewidth=0.2, edgecolor='none',
                           rstride=1, cstride=1)
    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, pad=0.08)
    cbar.set_label('Volatilidade Implícita (%)', color=COR_TEXTO, fontsize=8)
    cbar.ax.yaxis.set_tick_params(color=COR_TEXTO)
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color=COR_TEXTO)

    ax.set_xlabel('Moneyness (K/S)', color=COR_TEXTO, fontsize=8, labelpad=8)
    ax.set_ylabel('Maturidade (anos)', color=COR_TEXTO, fontsize=8, labelpad=8)
    ax.set_zlabel('Vol. Implícita (%)', color=COR_TEXTO, fontsize=8, labelpad=8)
    ax.set_title('Superfície de Volatilidade Implícita — B3\n(Skew Negativo + Smile)',
                 color=COR_TEXTO, fontsize=12, fontweight='bold', pad=12)
    ax.tick_params(colors=COR_TEXTO, labelsize=7)
    ax.xaxis.pane.fill = False; ax.yaxis.pane.fill = False; ax.zaxis.pane.fill = False
    ax.grid(True, alpha=0.2)
    ax.view_init(elev=25, azim=-45)

    fig.text(0.99, 0.005, 'Luiz Tiago Wilcke', ha='right', fontsize=7,
             color='#636e72', style='italic')
    caminho = os.path.join(diretorio, 'graf_02_superficie_vol.png')
    plt.savefig(caminho, dpi=150, bbox_inches='tight', facecolor=FUNDO_ESCURO)
    plt.close()
    print(f'[✓] Superfície de volatilidade salva: {caminho}')


if __name__ == '__main__':
    diretorio = '../graficos'
    arquivo_csv = '../resultados/dados_processados.csv'
    
    if os.path.exists(arquivo_csv):
        df = pd.read_csv(arquivo_csv)
        retornos = df['r'].values
        print(f'[i] Usando {len(retornos)} observações reais do Fortran.')
    else:
        np.random.seed(42)
        retornos = np.random.normal(0.0004, 0.013, 504)
        print('[!] Arquivo real não encontrado. Usando dados sintéticos.')
        
    plotar_volatilidade_comparativa(retornos, diretorio)
    plotar_superficie_volatilidade(diretorio)
    print('[✓] Módulo 02 — Gráficos de Volatilidade concluído.')
