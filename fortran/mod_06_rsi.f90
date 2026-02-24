!==============================================================================
! Módulo 06: RSI - Índice de Força Relativa
! Autor: Luiz Tiago Wilcke
! Descrição: Cálculo completo do RSI (Relative Strength Index) de Wilder,
!            RSI estocástico, RSI suavizado e interpretação de divergências
!            para análise técnica de ações brasileiras (B3).
!==============================================================================
MODULE mod_rsi
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! RSI de Wilder (14 períodos padrão)
   ! RS = EMA_ganhos / EMA_perdas
   ! RSI = 100 - 100/(1 + RS)
   !----------------------------------------------------------------------------
   SUBROUTINE calcular_rsi(precos, periodo, rsi_out)
      REAL(dp), INTENT(IN)  :: precos(:)
      INTEGER,  INTENT(IN)  :: periodo
      REAL(dp), INTENT(OUT) :: rsi_out(:)
      REAL(dp) :: ganho_med, perda_med, RS, ganho, perda
      REAL(dp) :: alpha
      INTEGER :: n, i
      n = SIZE(precos)
      rsi_out = 50.0_dp   ! valor neutro como padrão
      IF (n < periodo + 1) RETURN
      alpha = 1.0_dp / REAL(periodo, dp)
      ! Média inicial (SMA dos primeiros 'periodo' ganhos/perdas)
      ganho_med = 0.0_dp;  perda_med = 0.0_dp
      DO i = 2, periodo + 1
         ganho = MAX(precos(i) - precos(i-1), 0.0_dp)
         perda = MAX(precos(i-1) - precos(i), 0.0_dp)
         ganho_med = ganho_med + ganho
         perda_med = perda_med + perda
      END DO
      ganho_med = ganho_med / REAL(periodo, dp)
      perda_med = perda_med / REAL(periodo, dp)
      ! Primeiro RSI calculável
      IF (perda_med > 1.0D-10) THEN
         RS = ganho_med / perda_med
         rsi_out(periodo+1) = 100.0_dp - 100.0_dp / (1.0_dp + RS)
      ELSE
         rsi_out(periodo+1) = 100.0_dp
      END IF
      ! EMA de Wilder para o restante
      DO i = periodo + 2, n
         ganho = MAX(precos(i) - precos(i-1), 0.0_dp)
         perda = MAX(precos(i-1) - precos(i), 0.0_dp)
         ganho_med = (1.0_dp - alpha) * ganho_med + alpha * ganho
         perda_med = (1.0_dp - alpha) * perda_med + alpha * perda
         IF (perda_med > 1.0D-10) THEN
            RS = ganho_med / perda_med
            rsi_out(i) = 100.0_dp - 100.0_dp / (1.0_dp + RS)
         ELSE
            rsi_out(i) = 100.0_dp
         END IF
      END DO
   END SUBROUTINE calcular_rsi

   !----------------------------------------------------------------------------
   ! RSI Estocástico: aplica o estocástico sobre o RSI
   ! StochRSI = (RSI - RSI_min_periodo) / (RSI_max_periodo - RSI_min_periodo)
   !----------------------------------------------------------------------------
   SUBROUTINE rsi_estocastico(rsi_vals, periodo, stoch_rsi)
      REAL(dp), INTENT(IN)  :: rsi_vals(:)
      INTEGER,  INTENT(IN)  :: periodo
      REAL(dp), INTENT(OUT) :: stoch_rsi(:)
      REAL(dp) :: rsi_min, rsi_max, amplitude
      INTEGER :: n, i
      n = SIZE(rsi_vals)
      stoch_rsi = 0.0_dp
      DO i = periodo, n
         rsi_min   = MINVAL(rsi_vals(i-periodo+1:i))
         rsi_max   = MAXVAL(rsi_vals(i-periodo+1:i))
         amplitude = rsi_max - rsi_min
         IF (amplitude > 1.0D-10) THEN
            stoch_rsi(i) = (rsi_vals(i) - rsi_min) / amplitude
         ELSE
            stoch_rsi(i) = 0.5_dp
         END IF
      END DO
   END SUBROUTINE rsi_estocastico

   !----------------------------------------------------------------------------
   ! Zona de sobrecompra / sobrevenda
   ! retorna: +1 (sobrecompra RSI > 70), -1 (sobrevenda RSI < 30), 0 (neutro)
   !----------------------------------------------------------------------------
   INTEGER FUNCTION zona_rsi(valor_rsi, limiar_alto, limiar_baixo) RESULT(zona)
      REAL(dp), INTENT(IN) :: valor_rsi, limiar_alto, limiar_baixo
      IF (valor_rsi >= limiar_alto) THEN
         zona = 1
      ELSE IF (valor_rsi <= limiar_baixo) THEN
         zona = -1
      ELSE
         zona = 0
      END IF
   END FUNCTION zona_rsi

   !----------------------------------------------------------------------------
   ! Detecção de divergência bullish: preço faz mínimos mais baixos, RSI não
   !----------------------------------------------------------------------------
   SUBROUTINE detectar_divergencia_bullish(precos, rsi_vals, n_janela, divergencias)
      REAL(dp), INTENT(IN)  :: precos(:), rsi_vals(:)
      INTEGER,  INTENT(IN)  :: n_janela
      LOGICAL,  INTENT(OUT) :: divergencias(:)
      INTEGER :: n, i
      n = SIZE(precos)
      divergencias = .FALSE.
      DO i = n_janela + 1, n
         ! Preço faz novo mínimo mas RSI não confirma
         IF (precos(i) < MINVAL(precos(i-n_janela:i-1)) .AND. &
            rsi_vals(i) > MINVAL(rsi_vals(i-n_janela:i-1))) THEN
            divergencias(i) = .TRUE.
         END IF
      END DO
   END SUBROUTINE detectar_divergencia_bullish

   !----------------------------------------------------------------------------
   ! Detecção de divergência bearish: preço faz máximos mais altos, RSI não
   !----------------------------------------------------------------------------
   SUBROUTINE detectar_divergencia_bearish(precos, rsi_vals, n_janela, divergencias)
      REAL(dp), INTENT(IN)  :: precos(:), rsi_vals(:)
      INTEGER,  INTENT(IN)  :: n_janela
      LOGICAL,  INTENT(OUT) :: divergencias(:)
      INTEGER :: n, i
      n = SIZE(precos)
      divergencias = .FALSE.
      DO i = n_janela + 1, n
         IF (precos(i) > MAXVAL(precos(i-n_janela:i-1)) .AND. &
            rsi_vals(i) < MAXVAL(rsi_vals(i-n_janela:i-1))) THEN
            divergencias(i) = .TRUE.
         END IF
      END DO
   END SUBROUTINE detectar_divergencia_bearish

   !----------------------------------------------------------------------------
   ! RSI suavizado com EMA de 3 períodos (RSI Signal)
   !----------------------------------------------------------------------------
   SUBROUTINE rsi_suavizado(rsi_vals, periodo_suav, rsi_sinal)
      REAL(dp), INTENT(IN)  :: rsi_vals(:)
      INTEGER,  INTENT(IN)  :: periodo_suav
      REAL(dp), INTENT(OUT) :: rsi_sinal(:)
      REAL(dp) :: alpha
      INTEGER :: n, i
      n = SIZE(rsi_vals)
      alpha = 2.0_dp / REAL(periodo_suav + 1, dp)
      rsi_sinal = 0.0_dp
      rsi_sinal(periodo_suav) = SUM(rsi_vals(1:periodo_suav)) / REAL(periodo_suav, dp)
      DO i = periodo_suav + 1, n
         rsi_sinal(i) = alpha * rsi_vals(i) + (1.0_dp - alpha) * rsi_sinal(i-1)
      END DO
   END SUBROUTINE rsi_suavizado

   !----------------------------------------------------------------------------
   ! Histograma RSI (diferença entre RSI e sua linha de sinal)
   !----------------------------------------------------------------------------
   SUBROUTINE histograma_rsi(rsi_vals, rsi_sinal, histograma)
      REAL(dp), INTENT(IN)  :: rsi_vals(:), rsi_sinal(:)
      REAL(dp), INTENT(OUT) :: histograma(:)
      INTEGER :: n
      n = MIN(SIZE(rsi_vals), SIZE(rsi_sinal), SIZE(histograma))
      histograma(1:n) = rsi_vals(1:n) - rsi_sinal(1:n)
   END SUBROUTINE histograma_rsi

   !----------------------------------------------------------------------------
   ! RSI Ponderado por Volume (RSI-V)
   ! Pondera os ganhos e perdas pelo volume negociado
   !----------------------------------------------------------------------------
   SUBROUTINE rsi_ponderado_volume(precos, volumes, periodo, rsiv_out)
      REAL(dp), INTENT(IN)  :: precos(:), volumes(:)
      INTEGER,  INTENT(IN)  :: periodo
      REAL(dp), INTENT(OUT) :: rsiv_out(:)
      REAL(dp) :: ganho_vol, perda_vol, RS, ganho, perda, alpha
      INTEGER :: n, i
      n = MIN(SIZE(precos), SIZE(volumes))
      rsiv_out = 50.0_dp
      IF (n < periodo + 1) RETURN
      alpha = 1.0_dp / REAL(periodo, dp)
      ganho_vol = 0.0_dp;  perda_vol = 0.0_dp
      DO i = 2, periodo + 1
         ganho = MAX(precos(i) - precos(i-1), 0.0_dp) * volumes(i)
         perda = MAX(precos(i-1) - precos(i), 0.0_dp) * volumes(i)
         ganho_vol = ganho_vol + ganho
         perda_vol = perda_vol + perda
      END DO
      ganho_vol = ganho_vol / REAL(periodo, dp)
      perda_vol = perda_vol / REAL(periodo, dp)
      IF (perda_vol > 1.0D-10) THEN
         RS = ganho_vol / perda_vol
         rsiv_out(periodo+1) = 100.0_dp - 100.0_dp / (1.0_dp + RS)
      ELSE
         rsiv_out(periodo+1) = 100.0_dp
      END IF
      DO i = periodo + 2, n
         ganho = MAX(precos(i) - precos(i-1), 0.0_dp) * volumes(i)
         perda = MAX(precos(i-1) - precos(i), 0.0_dp) * volumes(i)
         ganho_vol = (1.0_dp - alpha) * ganho_vol + alpha * ganho
         perda_vol = (1.0_dp - alpha) * perda_vol + alpha * perda
         IF (perda_vol > 1.0D-10) THEN
            RS = ganho_vol / perda_vol
            rsiv_out(i) = 100.0_dp - 100.0_dp / (1.0_dp + RS)
         ELSE
            rsiv_out(i) = 100.0_dp
         END IF
      END DO
   END SUBROUTINE rsi_ponderado_volume

END MODULE mod_rsi
