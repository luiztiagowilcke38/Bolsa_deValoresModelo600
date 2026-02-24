!==============================================================================
! Módulo 38: Programa Principal — Orquestrador do Modelo Financeiro B3
! Autor: Luiz Tiago Wilcke
!==============================================================================
MODULE mod_principal
   USE mod_series_temporais
   USE mod_equacoes_diferenciais
   USE mod_black_scholes
   USE mod_volatilidade
   USE mod_media_movel
   USE mod_rsi
   USE mod_macd
   USE mod_bollinger
   USE mod_regressao
   USE mod_correlacao
   USE mod_covariancia
   USE mod_portfolio
   USE mod_monte_carlo
   USE mod_var
   USE mod_garch
   USE mod_arima
   USE mod_cadeias_markov
   USE mod_redes_neurais
   USE mod_kalman
   USE mod_ondas
   USE mod_fourier
   USE mod_entropia
   USE mod_fractal
   USE mod_caos
   USE mod_estocastico
   USE mod_difusao
   USE mod_heston
   USE mod_taxas_juros
   USE mod_derivativos
   USE mod_liquidez_sentimento_fatorial
   USE mod_otimizacao
   IMPLICIT NONE

   INTEGER, PARAMETER :: DP_L = KIND(1.0D0)
   INTEGER, PARAMETER :: DIAS_ANU_L = 252
   REAL(DP_L), PARAMETER :: SELIC_L = 0.1375_8
   REAL(DP_L), PARAMETER :: IPCA_L  = 0.050_8

CONTAINS

   SUBROUTINE gerar_dados_b3(serie, n_obs_in)
      TYPE(SerieFinanceira), INTENT(INOUT) :: serie
      INTEGER, INTENT(IN) :: n_obs_in
      REAL(DP_L), ALLOCATABLE :: traj(:,:)
      INTEGER :: i_c
      REAL(DP_L), PARAMETER :: S0_V = 130000.0_8
      REAL(DP_L), PARAMETER :: MU_V = 0.12_8
      REAL(DP_L), PARAMETER :: SG_V = 0.20_8

      IF (ALLOCATED(serie%preco)) DEALLOCATE(serie%preco)
      IF (ALLOCATED(serie%datas)) DEALLOCATE(serie%datas)
      ALLOCATE(serie%preco(n_obs_in), serie%datas(n_obs_in))
      ALLOCATE(traj(1, n_obs_in))

      CALL simular_mbg(S0_V, MU_V / REAL(DIAS_ANU_L, DP_L), &
         SG_V / SQRT(REAL(DIAS_ANU_L, DP_L)), &
         1.0_8, n_obs_in-1, 1, traj)

      serie%preco = traj(1, :)
      serie%num_observacoes = n_obs_in
      serie%ticker = 'IBOV'
      DO i_c = 1, n_obs_in
         serie%datas(i_c) = '2024-01-01'
      END DO
      IF (ALLOCATED(traj)) DEALLOCATE(traj)
   END SUBROUTINE gerar_dados_b3

   SUBROUTINE executar_pipeline(n_tot, f_out)
      INTEGER, INTENT(IN) :: n_tot
      CHARACTER(LEN=*), INTENT(IN) :: f_out
      TYPE(SerieFinanceira) :: ibov_s

      REAL(DP_L), ALLOCATABLE :: r(:), g(:), m50(:), rs(:), ml(:), ms(:), mh(:)
      REAL(DP_L), ALLOCATABLE :: bsu(:), bin(:), bme(:), kf(:), kp(:)
      REAL(DP_L), ALLOCATABLE :: fr(:), fi(:), fa(:), e9(:), e21(:), e200(:)

      REAL(DP_L) :: h_v, dfa_v, ent_h, efic_v, mdd_v, va_v, cva_v, sha_v
      INTEGER :: n_r, u_o, k_i

      WRITE(*,*) 'Iniciando Pipeline B3...'
      CALL gerar_dados_b3(ibov_s, n_tot)
      n_r = n_tot - 1

      ALLOCATE(r(n_r), g(n_r), m50(n_tot), rs(n_tot), ml(n_tot), ms(n_tot), mh(n_tot))
      ALLOCATE(bsu(n_tot), bin(n_tot), bme(n_tot), kf(n_tot), kp(n_tot))
      ALLOCATE(fr(n_r), fi(n_r), fa(n_r), e9(n_tot), e21(n_tot), e200(n_tot))

      CALL calcular_retornos_log(ibov_s)
      r = ibov_s%retorno

      CALL garch11(r, 0.00001_8, 0.08_8, 0.90_8, g)
      CALL sma(ibov_s%preco, 50, m50)
      CALL ema(ibov_s%preco, 9, e9)
      CALL ema(ibov_s%preco, 21, e21)
      CALL ema(ibov_s%preco, 200, e200)
      CALL calcular_rsi(ibov_s%preco, 14, rs)
      CALL calcular_macd(ibov_s%preco, 12, 26, 9, ml, ms, mh)
      CALL bandas_bollinger(ibov_s%preco, 20, 2.0_8, bme, bsu, bin)
      CALL filtro_kalman_escalar(ibov_s%preco, 0.99_8, 1000.0_8, 5000.0_8, &
         ibov_s%preco(1), 1.0D6, kf, kp)
      CALL dft(r, n_r, fr, fi)
      CALL amplitude_espectral(fr, fi, n_r, fa)
      CALL expoente_hurst_rs(r, n_r, 10, h_v)
      CALL dfa(r, n_r, 10, n_r/4, dfa_v)

      ent_h = entropia_retornos(r, n_r, 20)
      efic_v = eficiencia_mercado(ent_h, 20)
      mdd_v = maximo_drawdown(ibov_s%preco)
      va_v = var_historico(r, n_r, 0.95_8, 1.0_8)
      cva_v = cvar(r, n_r, 0.95_8, 1.0_8)
      sha_v = razao_sharpe(SUM(r)/REAL(n_r,DP_L), SELIC_L/REAL(DIAS_ANU_L,DP_L), &
         SQRT(SUM(r**2)/REAL(n_r,DP_L)*REAL(DIAS_ANU_L,DP_L)))

      OPEN(NEWUNIT=u_o, FILE=TRIM(f_out), STATUS='REPLACE', ACTION='WRITE')
      WRITE(u_o,'(A)') 'idx,p,r,v,m,e9,e21,rsi,macd,bu,bl,k'
      DO k_i = 1, n_r
         WRITE(u_o,'(I6,11(",",ES15.6))') k_i, ibov_s%preco(k_i), r(k_i), &
            SQRT(ABS(g(k_i)))*SQRT(REAL(DIAS_ANU_L,DP_L)), &
            m50(k_i), e9(k_i), e21(k_i), rs(k_i), ml(k_i), bsu(k_i), bin(k_i), kf(k_i)
      END DO
      CLOSE(u_o)

      WRITE(*,*) 'Sucesso: ', TRIM(f_out)
      CALL liberar_serie(ibov_s)
      DEALLOCATE(r, g, m50, rs, ml, ms, mh, bsu, bin, bme, kf, kp, fr, fi, fa, e9, e21, e200)
   END SUBROUTINE executar_pipeline

END MODULE mod_principal

PROGRAM main_prog
   USE mod_principal
   IMPLICIT NONE
   CALL executar_pipeline(1260, 'resultados/dados_processados.csv')
END PROGRAM main_prog
