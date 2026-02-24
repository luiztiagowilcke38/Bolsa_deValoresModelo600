!==============================================================================
! Módulo 07: MACD - Convergência e Divergência de Médias Móveis
! Autor: Luiz Tiago Wilcke
! Descrição: Implementação do MACD (Moving Average Convergence Divergence)
!            e indicadores derivados: linha de sinal, histograma, MACD-V,
!            MACD normalizado e detecção de divergências para B3.
!==============================================================================
MODULE mod_macd
   USE mod_media_movel, ONLY: ema
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! Calcula MACD clássico (Appel)
   ! MACD_linha = EMA_rapida - EMA_lenta  (padrão: 12 e 26)
   ! Sinal      = EMA(MACD_linha, 9)
   ! Histograma = MACD_linha - Sinal
   !----------------------------------------------------------------------------
   SUBROUTINE calcular_macd(precos, per_rapido, per_lento, per_sinal, &
      macd_linha, linha_sinal, histograma)
      REAL(dp), INTENT(IN)  :: precos(:)
      INTEGER,  INTENT(IN)  :: per_rapido, per_lento, per_sinal
      REAL(dp), INTENT(OUT) :: macd_linha(:), linha_sinal(:), histograma(:)
      REAL(dp), ALLOCATABLE :: ema_rap(:), ema_len(:)
      INTEGER :: n
      n = SIZE(precos)
      ALLOCATE(ema_rap(n), ema_len(n))
      CALL ema(precos,     per_rapido, ema_rap)
      CALL ema(precos,     per_lento,  ema_len)
      macd_linha   = ema_rap - ema_len
      CALL ema(macd_linha, per_sinal,  linha_sinal)
      histograma   = macd_linha - linha_sinal
      DEALLOCATE(ema_rap, ema_len)
   END SUBROUTINE calcular_macd

   !----------------------------------------------------------------------------
   ! Normalização do histograma pelo ATR (Average True Range)
   ! MACD-V = histograma / ATR * 100
   !----------------------------------------------------------------------------
   SUBROUTINE macd_normalizado(histograma, atr, macd_v)
      REAL(dp), INTENT(IN)  :: histograma(:), atr(:)
      REAL(dp), INTENT(OUT) :: macd_v(:)
      INTEGER :: n, i
      n = MIN(SIZE(histograma), SIZE(atr), SIZE(macd_v))
      DO i = 1, n
         IF (ABS(atr(i)) > 1.0D-10) THEN
            macd_v(i) = histograma(i) / atr(i) * 100.0_dp
         ELSE
            macd_v(i) = 0.0_dp
         END IF
      END DO
   END SUBROUTINE macd_normalizado

   !----------------------------------------------------------------------------
   ! Detecta cruzamento da linha MACD com a linha de sinal
   ! +1: MACD cruza sinal para cima (compra)
   ! -1: MACD cruza sinal para baixo (venda)
   !  0: sem cruzamento
   !----------------------------------------------------------------------------
   SUBROUTINE sinal_macd(macd_linha, linha_sinal, sinais)
      REAL(dp), INTENT(IN)  :: macd_linha(:), linha_sinal(:)
      INTEGER,  INTENT(OUT) :: sinais(:)
      INTEGER :: n, i
      n = MIN(SIZE(macd_linha), SIZE(linha_sinal), SIZE(sinais))
      sinais = 0
      DO i = 2, n
         IF (macd_linha(i-1) <= linha_sinal(i-1) .AND. &
            macd_linha(i)   >  linha_sinal(i)) THEN
            sinais(i) = 1    ! cruzamento bullish
         ELSE IF (macd_linha(i-1) >= linha_sinal(i-1) .AND. &
            macd_linha(i)   <  linha_sinal(i)) THEN
            sinais(i) = -1   ! cruzamento bearish
         END IF
      END DO
   END SUBROUTINE sinal_macd

   !----------------------------------------------------------------------------
   ! Cruzamento do MACD com a linha zero
   ! +1: cruzamento para cima (impulso bullish)
   ! -1: cruzamento para baixo (impulso bearish)
   !----------------------------------------------------------------------------
   SUBROUTINE cruzamento_zero_macd(macd_linha, sinais_zero)
      REAL(dp), INTENT(IN)  :: macd_linha(:)
      INTEGER,  INTENT(OUT) :: sinais_zero(:)
      INTEGER :: n, i
      n = SIZE(macd_linha)
      sinais_zero = 0
      DO i = 2, n
         IF (macd_linha(i-1) < 0.0_dp .AND. macd_linha(i) >= 0.0_dp) THEN
            sinais_zero(i) = 1
         ELSE IF (macd_linha(i-1) > 0.0_dp .AND. macd_linha(i) <= 0.0_dp) THEN
            sinais_zero(i) = -1
         END IF
      END DO
   END SUBROUTINE cruzamento_zero_macd

   !----------------------------------------------------------------------------
   ! Divergência bullish no histograma MACD:
   ! Histograma faz mínimos menos negativos enquanto preço faz mínimos mais baixos
   !----------------------------------------------------------------------------
   SUBROUTINE divergencia_bullish_macd(precos, histograma, janela, divergencias)
      REAL(dp), INTENT(IN)  :: precos(:), histograma(:)
      INTEGER,  INTENT(IN)  :: janela
      LOGICAL,  INTENT(OUT) :: divergencias(:)
      INTEGER :: n, i
      n = MIN(SIZE(precos), SIZE(histograma))
      divergencias = .FALSE.
      DO i = janela + 1, n
         IF (precos(i) < MINVAL(precos(i-janela:i-1)) .AND. &
            histograma(i) > MINVAL(histograma(i-janela:i-1)) .AND. &
            histograma(i) < 0.0_dp) THEN
            divergencias(i) = .TRUE.
         END IF
      END DO
   END SUBROUTINE divergencia_bullish_macd

   !----------------------------------------------------------------------------
   ! MACD de longo prazo (configuração semanal: 5, 10, 3 dias úteis)
   !----------------------------------------------------------------------------
   SUBROUTINE macd_semanal(precos, macd_linha, linha_sinal, histograma)
      REAL(dp), INTENT(IN)  :: precos(:)
      REAL(dp), INTENT(OUT) :: macd_linha(:), linha_sinal(:), histograma(:)
      CALL calcular_macd(precos, 5, 10, 3, macd_linha, linha_sinal, histograma)
   END SUBROUTINE macd_semanal

   !----------------------------------------------------------------------------
   ! Razão de força do momentum a partir do histograma
   ! MR = histograma / abs(macd_linha)  (em %)
   !----------------------------------------------------------------------------
   SUBROUTINE razao_momentum(macd_linha, histograma, razao)
      REAL(dp), INTENT(IN)  :: macd_linha(:), histograma(:)
      REAL(dp), INTENT(OUT) :: razao(:)
      INTEGER :: n, i
      n = MIN(SIZE(macd_linha), SIZE(histograma), SIZE(razao))
      DO i = 1, n
         IF (ABS(macd_linha(i)) > 1.0D-10) THEN
            razao(i) = (histograma(i) / ABS(macd_linha(i))) * 100.0_dp
         ELSE
            razao(i) = 0.0_dp
         END IF
      END DO
   END SUBROUTINE razao_momentum

   !----------------------------------------------------------------------------
   ! Número de períodos consecutivos com histograma positivo (momentum sustentado)
   !----------------------------------------------------------------------------
   INTEGER FUNCTION periodos_hist_positivo(histograma, indice_final) RESULT(np)
      REAL(dp), INTENT(IN) :: histograma(:)
      INTEGER,  INTENT(IN) :: indice_final
      INTEGER :: i
      np = 0
      i  = indice_final
      DO WHILE (i >= 1)
         IF (histograma(i) > 0.0_dp) THEN
            np = np + 1
            i  = i - 1
         ELSE
            EXIT
         END IF
      END DO
   END FUNCTION periodos_hist_positivo

END MODULE mod_macd
