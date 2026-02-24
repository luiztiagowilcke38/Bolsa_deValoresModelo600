#!/usr/bin/env python3
# ==============================================================================
# Módulo Python 01: Gráficos de Preços e Indicadores Técnicos
# Autor: Luiz Tiago Wilcke
# Descrição: Geração de gráficos de candlestick, séries de preços com EMAs,
#            Bandas de Bollinger e volume para ações da B3. Integra-se com
#            os resultados numéricos dos módulos Fortran.
# ==============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from datetime import datetime, timedelta
import warnings
warnings.filterwarnings('ignore')

# Paleta de cores do tema financeiro (dark mode)
FUNDO_ESCURO   = '#0d1117'
FUNDO_PAINEL   = '#161b22'
COR_ALTA       = '#00d084'   # verde para alta
COR_BAIXA      = '#ff4757'   # vermelho para baixa
COR_EMA9       = '#ffd32a'
COR_EMA21      = '#ff9f43'
COR_EMA200     = '#a29bfe'
COR_BB_SUP     = '#74b9ff'
COR_BB_INF     = '#74b9ff'
COR_BB_MED     = '#dfe6e9'
COR_VOLUME     = '#636e72'
COR_TEXTO      = '#e2e8f0'
COR_GRADE      = '#2d333b'


def configurar_estilo_dark():
    """Configura o estilo escuro global para todos os gráficos."""
    plt.rcParams.update({
        'figure.facecolor':  FUNDO_ESCURO,
        'axes.facecolor':    FUNDO_PAINEL,
        'axes.edgecolor':    COR_GRADE,
        'axes.labelcolor':   COR_TEXTO,
        'axes.titlecolor':   COR_TEXTO,
        'xtick.color':       COR_TEXTO,
        'ytick.color':       COR_TEXTO,
        'text.color':        COR_TEXTO,
        'grid.color':        COR_GRADE,
        'grid.linewidth':    0.5,
        'grid.alpha':        0.5,
        'legend.facecolor':  FUNDO_PAINEL,
        'legend.edgecolor':  COR_GRADE,
        'font.family':       'monospace',
        'font.size':         9,
    })


def gerar_dados_sinteticos_ohlcv(n: int = 252, S0: float = 130000.0,
                                   mu: float = 0.0004, sigma: float = 0.013,
                                   semente: int = 42) -> pd.DataFrame:
    """
    Gera dados OHLCV sintéticos via MBG para simular o IBOVESPA.

    Processo:
        S_t = S_{t-1} * exp((mu - sigma^2/2)*dt + sigma*sqrt(dt)*Z)
        High = S_t * exp(|epsilon1|),  Low = S_t * exp(-|epsilon2|)
    """
    np.random.seed(semente)
    retornos_diarios = np.random.normal(mu, sigma, n)
    log_precos = np.log(S0) + np.cumsum(retornos_diarios)
    fechamento = np.exp(log_precos)

    eps1 = np.abs(np.random.normal(0, sigma * 0.5, n))
    eps2 = np.abs(np.random.normal(0, sigma * 0.5, n))

    maximo    = fechamento * np.exp(eps1)
    minimo    = fechamento * np.exp(-eps2)
    abertura  = np.roll(fechamento, 1)
    abertura[0] = S0
    volume_base = np.random.lognormal(18.5, 0.5, n)

    # Datas dos dias úteis da B3 (2020-01-02 em diante)
    data_inicio = datetime(2020, 1, 2)
    datas = []
    d = data_inicio
    while len(datas) < n:
        if d.weekday() < 5:   # segunda a sexta
            datas.append(d)
        d += timedelta(days=1)

    return pd.DataFrame({
        'data':        datas,
        'abertura':    abertura,
        'maximo':      maximo,
        'minimo':      minimo,
        'fechamento':  fechamento,
        'volume':      volume_base,
        'retorno':     np.concatenate([[0], retornos_diarios[:-1]])
    })


def calcular_ema_python(serie: np.ndarray, janela: int) -> np.ndarray:
    """EMA (média móvel exponencial) em Python para uso nos gráficos."""
    alpha = 2.0 / (janela + 1)
    ema_vals = np.zeros_like(serie)
    ema_vals[janela-1] = np.mean(serie[:janela])
    for i in range(janela, len(serie)):
        ema_vals[i] = alpha * serie[i] + (1 - alpha) * ema_vals[i-1]
    return ema_vals


def calcular_bollinger(preco: np.ndarray, janela: int = 20,
                        k: float = 2.0):
    """Calcula Bandas de Bollinger."""
    media = np.array([np.mean(preco[max(0,i-janela):i+1])
                      for i in range(len(preco))])
    desvio = np.array([np.std(preco[max(0,i-janela):i+1], ddof=1)
                       for i in range(len(preco))])
    return media, media + k*desvio, media - k*desvio


def plotar_preco_com_emas(df: pd.DataFrame, caminho_saida: str,
                           ticker: str = 'IBOV'):
    """
    Gráfico principal: preço de fechamento com EMA9, EMA21, EMA200 e
    Bandas de Bollinger no painel superior; volume no painel inferior.

    Layout: Gridspec 3:1 (preço:volume)
    """
    configurar_estilo_dark()
    precos   = df['fechamento'].values
    n        = len(precos)
    indices  = np.arange(n)

    ema9   = calcular_ema_python(precos, 9)
    ema21  = calcular_ema_python(precos, 21)
    ema200 = calcular_ema_python(precos, 200)
    bb_med, bb_sup, bb_inf = calcular_bollinger(precos, 20, 2.0)

    fig = plt.figure(figsize=(16, 9), facecolor=FUNDO_ESCURO)
    gs  = GridSpec(4, 1, figure=fig, hspace=0.05,
                   height_ratios=[3, 0.6, 0.6, 0.8])

    ax_preco  = fig.add_subplot(gs[0])
    ax_vol    = fig.add_subplot(gs[3], sharex=ax_preco)

    # --- Painel de Preços ---
    cores_vela = np.where(precos >= df['abertura'].values, COR_ALTA, COR_BAIXA)
    for i in range(n):
        c = COR_ALTA if precos[i] >= df['abertura'].values[i] else COR_BAIXA
        # Corpo da vela
        corpo_min = min(precos[i], df['abertura'].values[i])
        corpo_max = max(precos[i], df['abertura'].values[i])
        if i % 3 == 0:   # plota apenas 1 em 3 velas para legibilidade
            ax_preco.bar(i, corpo_max - corpo_min, bottom=corpo_min,
                         color=c, alpha=0.8, width=0.6, linewidth=0)
            ax_preco.plot([i, i], [df['minimo'].values[i], df['maximo'].values[i]],
                          color=c, linewidth=0.5, alpha=0.6)

    # EMAs
    ax_preco.plot(indices[200:], ema200[200:], color=COR_EMA200, linewidth=1.2,
                  label='EMA 200', alpha=0.8)
    ax_preco.plot(indices[20:],  ema21[20:],   color=COR_EMA21,  linewidth=1.2,
                  label='EMA 21',  alpha=0.9)
    ax_preco.plot(indices[9:],   ema9[9:],     color=COR_EMA9,   linewidth=1.0,
                  label='EMA 9',   alpha=0.9)

    # Bollinger
    ax_preco.fill_between(indices[20:], bb_sup[20:], bb_inf[20:],
                           alpha=0.06, color=COR_BB_SUP, label='BB ±2σ')
    ax_preco.plot(indices[20:], bb_sup[20:], color=COR_BB_SUP, linewidth=0.7, alpha=0.5)
    ax_preco.plot(indices[20:], bb_inf[20:], color=COR_BB_INF, linewidth=0.7, alpha=0.5)
    ax_preco.plot(indices[20:], bb_med[20:], color=COR_BB_MED, linewidth=0.5,
                  linestyle='--', alpha=0.4)

    ax_preco.set_title(f'  {ticker} — Preço, EMAs e Bandas de Bollinger',
                        loc='left', fontsize=13, fontweight='bold', pad=10)
    ax_preco.set_ylabel('Pontos (IBOV)', fontsize=9)
    ax_preco.legend(loc='upper left', fontsize=8, framealpha=0.6, ncol=4)
    ax_preco.grid(True, alpha=0.25)
    ax_preco.yaxis.set_major_formatter(plt.FuncFormatter(
        lambda x, _: f'{x:,.0f}'))
    plt.setp(ax_preco.get_xticklabels(), visible=False)

    # --- Painel de Volume ---
    cores_vol = np.where(df['retorno'].values >= 0, COR_ALTA, COR_BAIXA)
    ax_vol.bar(indices, df['volume'].values / 1e9, color=cores_vol,
               alpha=0.7, linewidth=0)
    ax_vol.set_ylabel('Vol. (bi R$)', fontsize=8)
    ax_vol.grid(True, alpha=0.25)
    ax_vol.set_ylim(bottom=0)

    # Anotação de autoria
    fig.text(0.99, 0.005, 'Luiz Tiago Wilcke — Modelo Matemático B3',
             ha='right', fontsize=7, color='#636e72', style='italic')

    os.makedirs(os.path.dirname(caminho_saida) or '.', exist_ok=True)
    plt.savefig(caminho_saida, dpi=150, bbox_inches='tight',
                facecolor=FUNDO_ESCURO)
    plt.close()
    print(f'[✓] Gráfico de preços salvo: {caminho_saida}')


def plotar_dashboard_precos(df: pd.DataFrame, diretorio: str):
    """Gera o painel completo de preços para múltiplos ativos simulados."""
    configurar_estilo_dark()
    tickers_b3 = ['PETR4', 'VALE3', 'ITUB4', 'BBDC4', 'ABEV3',
                   'MGLU3', 'WEGE3', 'BBAS3', 'RENT3', 'SUZB3']
    np.random.seed(99)
    precos_iniciais = [35, 85, 32, 28, 14, 8, 45, 55, 65, 50]

    fig, axes = plt.subplots(2, 5, figsize=(20, 8), facecolor=FUNDO_ESCURO)
    fig.suptitle('Painel de Ações — B3 (Simulação)',
                  fontsize=14, fontweight='bold', color=COR_TEXTO, y=1.01)

    for ax, ticker, p0 in zip(axes.flat, tickers_b3, precos_iniciais):
        n = len(df)
        ret  = np.random.normal(0.0004, 0.015, n)
        preco = p0 * np.exp(np.cumsum(ret))
        ema20 = calcular_ema_python(preco, 20)
        cor_linha = COR_ALTA if preco[-1] > p0 else COR_BAIXA
        variacao = (preco[-1]/p0 - 1) * 100

        ax.plot(preco, color=cor_linha, linewidth=1.2, alpha=0.9)
        ax.plot(ema20, color='#fdcb6e', linewidth=0.8, alpha=0.7, linestyle='--')
        ax.fill_between(range(n), preco, p0, alpha=0.08, color=cor_linha)
        ax.set_title(f'{ticker}\n{preco[-1]:.2f}  ({variacao:+.1f}%)',
                     fontsize=9, fontweight='bold',
                     color=COR_ALTA if variacao > 0 else COR_BAIXA)
        ax.set_facecolor(FUNDO_PAINEL)
        ax.grid(True, alpha=0.2)
        ax.tick_params(labelsize=7)
        for spine in ax.spines.values():
            spine.set_edgecolor(COR_GRADE)

    plt.tight_layout()
    caminho = os.path.join(diretorio, 'graf_01_dashboard_precos.png')
    os.makedirs(diretorio, exist_ok=True)
    plt.savefig(caminho, dpi=150, bbox_inches='tight', facecolor=FUNDO_ESCURO)
    plt.close()
    print(f'[✓] Dashboard de preços salvo: {caminho}')


def carregar_dados_reais(caminho):
    """Carrega dados processados pelo motor Fortran."""
    if not os.path.exists(caminho):
        print(f"[!] Erro: Arquivo {caminho} não encontrado. Usando sintéticos.")
        return gerar_dados_sinteticos_ohlcv(n=504)
    
    df = pd.read_csv(caminho)
    # Converte colunas para numerico, forçando erros para NaN
    for col in df.columns:
        if col != 'idx':
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Preenche NaNs resultantes de asteriscos do Fortran
    df = df.ffill().bfill().fillna(0)
    
    # Mapeia colunas do Fortran para o que o script espera
    df = df.rename(columns={
        'p': 'fechamento', 'r': 'retorno', 'bu': 'bb_sup',
        'bl': 'bb_inf', 'm': 'bb_med', 'k': 'kalman'
    })
    # Gera OHLC fictício
    df['abertura'] = df['fechamento'].shift(1).fillna(df['fechamento'])
    df['maximo']   = df['fechamento'] * 1.005
    df['minimo']   = df['fechamento'] * 0.995
    df['volume']   = np.random.lognormal(18.5, 0.5, len(df))
    
    # Garante que temos datas validas
    if 'data' not in df.columns:
        d_start = datetime(2020, 1, 2)
        df['data'] = [d_start + timedelta(days=i) for i in range(len(df))]
        
    return df

if __name__ == '__main__':
    diretorio_graficos = '../graficos'
    arquivo_csv = '../resultados/dados_processados.csv'
    
    # Carrega dados reais do Fortran
    df = carregar_dados_reais(arquivo_csv)
    
    # Gera Gráficos
    plotar_preco_com_emas(df, os.path.join(diretorio_graficos, 'graf_01_precos_emas.png'))
    plotar_dashboard_precos(df, diretorio_graficos)
    print('[✓] Módulo 01 — Gráficos de Preços (com dados reais) concluído.')
