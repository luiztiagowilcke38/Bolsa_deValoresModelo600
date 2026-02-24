!==============================================================================
! Módulo 05: Médias Móveis e Suavização
! Autor: Luiz Tiago Wilcke
! Descrição: Implementa SMA, EMA, WMA, HMA, TEMA e médias Kaufman adaptativas
!            para identificação de tendências em ativos da B3.
!==============================================================================
MODULE mod_media_movel
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! Média Móvel Simples (SMA) — janela deslizante
   ! SMA_t = (1/n) * sum_{i=0}^{n-1} P_{t-i}
   !----------------------------------------------------------------------------
   SUBROUTINE sma(precos, janela, sma_out)
      REAL(dp), INTENT(IN)  :: precos(:)
      INTEGER,  INTENT(IN)  :: janela
      REAL(dp), INTENT(OUT) :: sma_out(:)
      INTEGER :: n, i
      n = SIZE(precos)
      sma_out = 0.0_dp
      DO i = janela, n
         sma_out(i) = SUM(precos(i-janela+1:i)) / REAL(janela, dp)
      END DO
   END SUBROUTINE sma

   !----------------------------------------------------------------------------
   ! Média Móvel Exponencial (EMA)
   ! EMA_t = alpha*P_t + (1-alpha)*EMA_{t-1},  alpha = 2/(n+1)
   !----------------------------------------------------------------------------
   SUBROUTINE ema(precos, janela, ema_out)
      REAL(dp), INTENT(IN)  :: precos(:)
      INTEGER,  INTENT(IN)  :: janela
      REAL(dp), INTENT(OUT) :: ema_out(:)
      REAL(dp) :: alpha
      INTEGER :: n, i
      n = SIZE(precos)
      alpha = 2.0_dp / REAL(janela + 1, dp)
      ema_out = 0.0_dp
      ! Inicializa com SMA dos primeiros 'janela' pontos
      ema_out(janela) = SUM(precos(1:janela)) / REAL(janela, dp)
      DO i = janela+1, n
         ema_out(i) = alpha * precos(i) + (1.0_dp - alpha) * ema_out(i-1)
      END DO
   END SUBROUTINE ema

   !----------------------------------------------------------------------------
   ! Média Móvel Ponderada (WMA) — pesos lineares (mais recente tem peso maior)
   ! WMA_t = sum_{i=1}^{n} i*P_{t-n+i} / sum_{i=1}^{n} i
   !----------------------------------------------------------------------------
   SUBROUTINE wma(precos, janela, wma_out)
      REAL(dp), INTENT(IN)  :: precos(:)
      INTEGER,  INTENT(IN)  :: janela
      REAL(dp), INTENT(OUT) :: wma_out(:)
      REAL(dp) :: soma_pesos, soma_ponderada
      INTEGER :: n, i, j
      n = SIZE(precos)
      soma_pesos = REAL(janela * (janela + 1), dp) / 2.0_dp
      wma_out = 0.0_dp
      DO i = janela, n
         soma_ponderada = 0.0_dp
         DO j = 1, janela
            soma_ponderada = soma_ponderada + REAL(j, dp) * precos(i - janela + j)
         END DO
         wma_out(i) = soma_ponderada / soma_pesos
      END DO
   END SUBROUTINE wma

   !----------------------------------------------------------------------------
   ! Média Móvel de Kaufman Adaptativa (KAMA)
   ! Ajusta a velocidade de resposta com base no índice de eficiência (ER)
   ! ER = |P_t - P_{t-n}| / sum|P_i - P_{i-1}|
   ! SC = [ER*(fast - slow) + slow]^2
   ! KAMA_t = KAMA_{t-1} + SC*(P_t - KAMA_{t-1})
   !----------------------------------------------------------------------------
   SUBROUTINE kama(precos, n_ef, n_rapido, n_lento, kama_out)
      REAL(dp), INTENT(IN)  :: precos(:)
      INTEGER,  INTENT(IN)  :: n_ef, n_rapido, n_lento
      REAL(dp), INTENT(OUT) :: kama_out(:)
      REAL(dp) :: fast_sc, slow_sc, ER, SC, soma_vol
      INTEGER :: n, i, j
      n = SIZE(precos)
      fast_sc = 2.0_dp / REAL(n_rapido + 1, dp)
      slow_sc = 2.0_dp / REAL(n_lento  + 1, dp)
      kama_out = 0.0_dp
      kama_out(n_ef) = precos(n_ef)
      DO i = n_ef+1, n
         soma_vol = 0.0_dp
         DO j = i-n_ef+1, i
            soma_vol = soma_vol + ABS(precos(j) - precos(j-1))
         END DO
         IF (soma_vol > 1.0D-10) THEN
            ER = ABS(precos(i) - precos(i-n_ef)) / soma_vol
         ELSE
            ER = 0.0_dp
         END IF
         SC = (ER * (fast_sc - slow_sc) + slow_sc)**2
         kama_out(i) = kama_out(i-1) + SC * (precos(i) - kama_out(i-1))
      END DO
   END SUBROUTINE kama

   !----------------------------------------------------------------------------
   ! Triple EMA (TEMA) — reduz lag da EMA
   ! TEMA = 3*EMA1 - 3*EMA2 + EMA3
   !----------------------------------------------------------------------------
   SUBROUTINE tema(precos, janela, tema_out)
      REAL(dp), INTENT(IN)  :: precos(:)
      INTEGER,  INTENT(IN)  :: janela
      REAL(dp), INTENT(OUT) :: tema_out(:)
      REAL(dp), ALLOCATABLE :: ema1(:), ema2(:), ema3(:)
      INTEGER :: n
      n = SIZE(precos)
      ALLOCATE(ema1(n), ema2(n), ema3(n))
      CALL ema(precos, janela, ema1)
      CALL ema(ema1,   janela, ema2)
      CALL ema(ema2,   janela, ema3)
      tema_out = 3.0_dp * ema1 - 3.0_dp * ema2 + ema3
      DEALLOCATE(ema1, ema2, ema3)
   END SUBROUTINE tema

   !----------------------------------------------------------------------------
   ! Hull Moving Average (HMA) — ainda mais rápida que WMA
   ! HMA(n) = WMA(2*WMA(n/2) - WMA(n), sqrt(n))
   !----------------------------------------------------------------------------
   SUBROUTINE hma(precos, janela, hma_out)
      REAL(dp), INTENT(IN)  :: precos(:)
      INTEGER,  INTENT(IN)  :: janela
      REAL(dp), INTENT(OUT) :: hma_out(:)
      REAL(dp), ALLOCATABLE :: wma_n(:), wma_meia(:), serie_diff(:)
      INTEGER :: n, janela_hma
      n = SIZE(precos)
      janela_hma = INT(SQRT(REAL(janela, dp)))
      ALLOCATE(wma_n(n), wma_meia(n), serie_diff(n))
      CALL wma(precos, janela,        wma_n)
      CALL wma(precos, janela/2,      wma_meia)
      serie_diff = 2.0_dp * wma_meia - wma_n
      CALL wma(serie_diff, janela_hma, hma_out)
      DEALLOCATE(wma_n, wma_meia, serie_diff)
   END SUBROUTINE hma

   !----------------------------------------------------------------------------
   ! Cruzamento de médias: retorna sinal de compra/venda
   ! +1 (compra) quando SMA_rapida cruza SMA_lenta para cima
   ! -1 (venda) quando cruza para baixo
   !  0 sem cruzamento
   !----------------------------------------------------------------------------
   SUBROUTINE sinal_cruzamento(media_rapida, media_lenta, sinais)
      REAL(dp), INTENT(IN)  :: media_rapida(:), media_lenta(:)
      INTEGER,  INTENT(OUT) :: sinais(:)
      INTEGER :: n, i
      n = MIN(SIZE(media_rapida), SIZE(media_lenta), SIZE(sinais))
      sinais = 0
      DO i = 2, n
         IF (media_rapida(i-1) < media_lenta(i-1) .AND. &
            media_rapida(i)   > media_lenta(i)) THEN
            sinais(i) = 1      ! cruzamento para cima — compra
         ELSE IF (media_rapida(i-1) > media_lenta(i-1) .AND. &
            media_rapida(i)   < media_lenta(i)) THEN
            sinais(i) = -1     ! cruzamento para baixo — venda
         END IF
      END DO
   END SUBROUTINE sinal_cruzamento

   !----------------------------------------------------------------------------
   ! Bandas de Médias: calcula canal superior e inferior
   ! superior = media + k*sigma
   ! inferior = media - k*sigma
   !----------------------------------------------------------------------------
   SUBROUTINE bandas_media(precos, janela, k, media, superior, inferior)
      REAL(dp), INTENT(IN)  :: precos(:)
      INTEGER,  INTENT(IN)  :: janela
      REAL(dp), INTENT(IN)  :: k
      REAL(dp), INTENT(OUT) :: media(:), superior(:), inferior(:)
      REAL(dp) :: sigma
      INTEGER :: n, i
      n = SIZE(precos)
      CALL sma(precos, janela, media)
      DO i = janela, n
         sigma = SQRT(SUM((precos(i-janela+1:i) - media(i))**2) / REAL(janela-1, dp))
         superior(i) = media(i) + k * sigma
         inferior(i) = media(i) - k * sigma
      END DO
   END SUBROUTINE bandas_media

   !----------------------------------------------------------------------------
   ! Detecta Golden Cross (EMA9 > EMA21 > EMA200)
   !----------------------------------------------------------------------------
   LOGICAL FUNCTION golden_cross(ema_rapida, ema_media, ema_lenta, indice)
      REAL(dp), INTENT(IN) :: ema_rapida(:), ema_media(:), ema_lenta(:)
      INTEGER,  INTENT(IN) :: indice
      golden_cross = (ema_rapida(indice) > ema_media(indice)) .AND. &
         (ema_media(indice)  > ema_lenta(indice))
   END FUNCTION golden_cross

END MODULE mod_media_movel
