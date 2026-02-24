!==============================================================================
! Módulo 08: Bandas de Bollinger e Canais de Volatilidade
! Autor: Luiz Tiago Wilcke
! Descrição: Implementa Bandas de Bollinger, Keltner Channels, Donchian,
!            %B, Largura de Banda e ATR para análise de volatilidade
!            relativa e reversão à média em ações da B3.
!==============================================================================
MODULE mod_bollinger
   USE mod_media_movel, ONLY: sma, ema
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! Bandas de Bollinger clássicas
   ! Banda_sup = SMA + k*dp_sigma
   ! Banda_inf = SMA - k*dp_sigma
   ! k = 2.0 (padrão Bollinger)
   !----------------------------------------------------------------------------
   SUBROUTINE bandas_bollinger(precos, janela, k_sigma, &
      banda_media, banda_sup, banda_inf)
      REAL(dp), INTENT(IN)  :: precos(:)
      INTEGER,  INTENT(IN)  :: janela
      REAL(dp), INTENT(IN)  :: k_sigma
      REAL(dp), INTENT(OUT) :: banda_media(:), banda_sup(:), banda_inf(:)
      REAL(dp) :: desvio
      INTEGER :: n, i
      n = SIZE(precos)
      CALL sma(precos, janela, banda_media)
      banda_sup = 0.0_dp;  banda_inf = 0.0_dp
      DO i = janela, n
         desvio = SQRT(SUM((precos(i-janela+1:i) - banda_media(i))**2) / &
            REAL(janela - 1, dp))
         banda_sup(i) = banda_media(i) + k_sigma * desvio
         banda_inf(i) = banda_media(i) - k_sigma * desvio
      END DO
   END SUBROUTINE bandas_bollinger

   !----------------------------------------------------------------------------
   ! %B: posição do preço dentro das bandas
   ! %B = (P - Banda_inf) / (Banda_sup - Banda_inf)
   ! %B < 0: abaixo da banda inferior; %B > 1: acima da superior
   !----------------------------------------------------------------------------
   SUBROUTINE percentual_b(precos, banda_sup, banda_inf, pct_b)
      REAL(dp), INTENT(IN)  :: precos(:), banda_sup(:), banda_inf(:)
      REAL(dp), INTENT(OUT) :: pct_b(:)
      REAL(dp) :: amplitude
      INTEGER :: n, i
      n = MIN(SIZE(precos), SIZE(banda_sup), SIZE(banda_inf), SIZE(pct_b))
      pct_b = 0.0_dp
      DO i = 1, n
         amplitude = banda_sup(i) - banda_inf(i)
         IF (amplitude > 1.0D-10) THEN
            pct_b(i) = (precos(i) - banda_inf(i)) / amplitude
         ELSE
            pct_b(i) = 0.5_dp
         END IF
      END DO
   END SUBROUTINE percentual_b

   !----------------------------------------------------------------------------
   ! Largura de Banda (Bandwidth) — normalizada
   ! BW = (Banda_sup - Banda_inf) / Banda_media
   !----------------------------------------------------------------------------
   SUBROUTINE largura_banda(banda_sup, banda_inf, banda_media, bw)
      REAL(dp), INTENT(IN)  :: banda_sup(:), banda_inf(:), banda_media(:)
      REAL(dp), INTENT(OUT) :: bw(:)
      INTEGER :: n, i
      n = MIN(SIZE(banda_sup), SIZE(banda_inf), SIZE(banda_media), SIZE(bw))
      DO i = 1, n
         IF (ABS(banda_media(i)) > 1.0D-10) THEN
            bw(i) = (banda_sup(i) - banda_inf(i)) / banda_media(i)
         ELSE
            bw(i) = 0.0_dp
         END IF
      END DO
   END SUBROUTINE largura_banda

   !----------------------------------------------------------------------------
   ! ATR — Average True Range de Wilder
   ! TR = max(Hi-Lo, |Hi-Ci-1|, |Lo-Ci-1|)
   ! ATR = EMA(TR, periodo)
   !----------------------------------------------------------------------------
   SUBROUTINE atr_wilder(maximos, minimos, fechamentos, periodo, atr)
      REAL(dp), INTENT(IN)  :: maximos(:), minimos(:), fechamentos(:)
      INTEGER,  INTENT(IN)  :: periodo
      REAL(dp), INTENT(OUT) :: atr(:)
      REAL(dp) :: TR, alpha
      REAL(dp), ALLOCATABLE :: tr_serie(:)
      INTEGER :: n, i
      n = SIZE(fechamentos)
      ALLOCATE(tr_serie(n))
      tr_serie(1) = maximos(1) - minimos(1)
      DO i = 2, n
         TR = MAX(maximos(i) - minimos(i), &
            ABS(maximos(i) - fechamentos(i-1)), &
            ABS(minimos(i) - fechamentos(i-1)))
         tr_serie(i) = TR
      END DO
      alpha = 1.0_dp / REAL(periodo, dp)
      atr(periodo) = SUM(tr_serie(1:periodo)) / REAL(periodo, dp)
      DO i = periodo + 1, n
         atr(i) = (1.0_dp - alpha) * atr(i-1) + alpha * tr_serie(i)
      END DO
      DEALLOCATE(tr_serie)
   END SUBROUTINE atr_wilder

   !----------------------------------------------------------------------------
   ! Canal de Keltner (EMA ± multiplicador * ATR)
   ! Sup = EMA + fator * ATR
   ! Inf = EMA - fator * ATR
   !----------------------------------------------------------------------------
   SUBROUTINE canal_keltner(precos, maximos, minimos, fechamentos, &
      per_ema, per_atr, fator, kc_sup, kc_med, kc_inf)
      REAL(dp), INTENT(IN)  :: precos(:), maximos(:), minimos(:), fechamentos(:)
      INTEGER,  INTENT(IN)  :: per_ema, per_atr
      REAL(dp), INTENT(IN)  :: fator
      REAL(dp), INTENT(OUT) :: kc_sup(:), kc_med(:), kc_inf(:)
      REAL(dp), ALLOCATABLE :: atr_vals(:)
      INTEGER :: n
      n = SIZE(precos)
      ALLOCATE(atr_vals(n))
      atr_vals = 0.0_dp
      CALL ema(precos, per_ema, kc_med)
      CALL atr_wilder(maximos, minimos, fechamentos, per_atr, atr_vals)
      kc_sup = kc_med + fator * atr_vals
      kc_inf = kc_med - fator * atr_vals
      DEALLOCATE(atr_vals)
   END SUBROUTINE canal_keltner

   !----------------------------------------------------------------------------
   ! Canal de Donchian (máximo e mínimo dos últimos N períodos)
   !----------------------------------------------------------------------------
   SUBROUTINE canal_donchian(maximos, minimos, janela, don_sup, don_inf, don_med)
      REAL(dp), INTENT(IN)  :: maximos(:), minimos(:)
      INTEGER,  INTENT(IN)  :: janela
      REAL(dp), INTENT(OUT) :: don_sup(:), don_inf(:), don_med(:)
      INTEGER :: n, i
      n = MIN(SIZE(maximos), SIZE(minimos))
      don_sup = 0.0_dp;  don_inf = 0.0_dp;  don_med = 0.0_dp
      DO i = janela, n
         don_sup(i) = MAXVAL(maximos(i-janela+1:i))
         don_inf(i) = MINVAL(minimos(i-janela+1:i))
         don_med(i) = (don_sup(i) + don_inf(i)) / 2.0_dp
      END DO
   END SUBROUTINE canal_donchian

   !----------------------------------------------------------------------------
   ! Squeeze do Bollinger dentro de Keltner (período de baixa volatilidade)
   ! True se BB estiver dentro do KC (compressão de volatilidade)
   !----------------------------------------------------------------------------
   LOGICAL FUNCTION squeeze_momentun(bb_sup, bb_inf, kc_sup, kc_inf) RESULT(sq)
      REAL(dp), INTENT(IN) :: bb_sup, bb_inf, kc_sup, kc_inf
      sq = (bb_sup < kc_sup) .AND. (bb_inf > kc_inf)
   END FUNCTION squeeze_momentun

   !----------------------------------------------------------------------------
   ! Índice de Volatilidade Relativa (VIX simplificado para ativo individual)
   ! VRI = BW / media_movel(BW, n)
   !----------------------------------------------------------------------------
   SUBROUTINE vri(largura_banda_vals, janela, vri_out)
      REAL(dp), INTENT(IN)  :: largura_banda_vals(:)
      INTEGER,  INTENT(IN)  :: janela
      REAL(dp), INTENT(OUT) :: vri_out(:)
      REAL(dp), ALLOCATABLE :: media_bw(:)
      INTEGER :: n
      n = SIZE(largura_banda_vals)
      ALLOCATE(media_bw(n))
      CALL sma(largura_banda_vals, janela, media_bw)
      WHERE (ABS(media_bw) > 1.0D-10)
         vri_out = largura_banda_vals / media_bw
      ELSEWHERE
         vri_out = 1.0_dp
      END WHERE
      DEALLOCATE(media_bw)
   END SUBROUTINE vri

END MODULE mod_bollinger
