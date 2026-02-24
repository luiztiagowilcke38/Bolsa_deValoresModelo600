!==============================================================================
! Módulo 04: Volatilidade Financeira
! Autor: Luiz Tiago Wilcke
! Descrição: Modelos de volatilidade histórica, realizada, condicional e
!            parkinsoniana para ativos da B3 (Bolsa brasileira).
!==============================================================================
MODULE mod_volatilidade
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)
   REAL(dp), PARAMETER :: DIAS_ANO = 252.0_dp   ! dias úteis brasileiros

CONTAINS

   !----------------------------------------------------------------------------
   ! Volatilidade histórica por janela deslizante (desvio padrão dos retornos)
   ! sigma = std(r_t) * sqrt(252)
   !----------------------------------------------------------------------------
   SUBROUTINE volatilidade_historica(retornos, janela, vol_hist)
      REAL(dp), INTENT(IN)  :: retornos(:)
      INTEGER,  INTENT(IN)  :: janela
      REAL(dp), INTENT(OUT) :: vol_hist(:)
      INTEGER :: n, i
      REAL(dp) :: media, variancia
      n = SIZE(retornos)
      vol_hist = 0.0_dp
      DO i = janela, n
         media = SUM(retornos(i-janela+1:i)) / REAL(janela, dp)
         variancia = SUM((retornos(i-janela+1:i) - media)**2) / REAL(janela-1, dp)
         vol_hist(i) = SQRT(variancia * DIAS_ANO)
      END DO
   END SUBROUTINE volatilidade_historica

   !----------------------------------------------------------------------------
   ! Volatilidade de Parkinson (estimador de range diário)
   ! sigma_P = sqrt[ (1/(4*N*ln2)) * sum(ln(Hi/Li))^2 ]
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION volatilidade_parkinson(maximos, minimos, n) RESULT(vol_P)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: maximos(n), minimos(n)
      REAL(dp) :: soma
      INTEGER :: i
      soma = 0.0_dp
      DO i = 1, n
         IF (minimos(i) > 0.0_dp) &
            soma = soma + (LOG(maximos(i) / minimos(i)))**2
      END DO
      vol_P = SQRT(soma / (4.0_dp * REAL(n, dp) * LOG(2.0_dp))) * SQRT(DIAS_ANO)
   END FUNCTION volatilidade_parkinson

   !----------------------------------------------------------------------------
   ! Estimador Garman-Klass (usa open, high, low, close)
   ! sigma^2 = (1/N) * sum[ 0.5*(ln(H/L))^2 - (2ln2-1)*(ln(C/O))^2 ]
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION volatilidade_garman_klass(abert, maxim, minim, fech, n) &
      RESULT(vol_GK)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: abert(n), maxim(n), minim(n), fech(n)
      REAL(dp) :: soma, termo1, termo2
      INTEGER :: i
      soma = 0.0_dp
      DO i = 1, n
         IF (abert(i) > 0.0_dp .AND. minim(i) > 0.0_dp) THEN
            termo1 = 0.5_dp * (LOG(maxim(i) / minim(i)))**2
            termo2 = (2.0_dp * LOG(2.0_dp) - 1.0_dp) * (LOG(fech(i) / abert(i)))**2
            soma = soma + termo1 - termo2
         END IF
      END DO
      vol_GK = SQRT(soma / REAL(n, dp) * DIAS_ANO)
   END FUNCTION volatilidade_garman_klass

   !----------------------------------------------------------------------------
   ! Volatilidade realizada (soma de retornos ao quadrado intradiários)
   ! RV_t = sum_{i=1}^{M} r_{t,i}^2
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION volatilidade_realizada(retornos_intra, M) RESULT(RV)
      INTEGER,  INTENT(IN) :: M
      REAL(dp), INTENT(IN) :: retornos_intra(M)
      RV = SQRT(SUM(retornos_intra**2) * DIAS_ANO)
   END FUNCTION volatilidade_realizada

   !----------------------------------------------------------------------------
   ! Índice VIX simplificado (variância esperada via aproximação)
   ! VIX = 100 * sqrt(E[sigma^2] * 252)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION indice_vix_estimado(variancia_condicional) RESULT(vix)
      REAL(dp), INTENT(IN) :: variancia_condicional
      vix = 100.0_dp * SQRT(variancia_condicional * DIAS_ANO)
   END FUNCTION indice_vix_estimado

   !----------------------------------------------------------------------------
   ! Decaimento exponencial de volatilidade (EWMA - RiskMetrics)
   ! sigma_t^2 = lambda*sigma_{t-1}^2 + (1-lambda)*r_{t-1}^2
   !----------------------------------------------------------------------------
   SUBROUTINE ewma_volatilidade(retornos, lambda, sigma0_sq, sigma_sq)
      REAL(dp), INTENT(IN)  :: retornos(:)
      REAL(dp), INTENT(IN)  :: lambda, sigma0_sq
      REAL(dp), INTENT(OUT) :: sigma_sq(:)
      INTEGER :: i, n
      n = SIZE(retornos)
      sigma_sq(1) = sigma0_sq
      DO i = 2, n
         sigma_sq(i) = lambda * sigma_sq(i-1) + (1.0_dp - lambda) * retornos(i-1)**2
      END DO
   END SUBROUTINE ewma_volatilidade

   !----------------------------------------------------------------------------
   ! Razão de Sharpe (retorno ajustado ao risco)
   ! S = (mu_port - r_livre) / sigma_port
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION razao_sharpe(retorno_medio, taxa_livre, vol_anual) RESULT(S)
      REAL(dp), INTENT(IN) :: retorno_medio, taxa_livre, vol_anual
      S = (retorno_medio * DIAS_ANO - taxa_livre) / (vol_anual + 1.0D-15)
   END FUNCTION razao_sharpe

   !----------------------------------------------------------------------------
   ! Razão de Sortino (penaliza apenas volatilidade negativa)
   ! S_sortino = (mu - r_livre) / sigma_downside
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION razao_sortino(retornos, taxa_livre) RESULT(sortino)
      REAL(dp), INTENT(IN) :: retornos(:), taxa_livre
      REAL(dp) :: media, downside_var, retorno_diario_livre
      INTEGER :: n, i
      n = SIZE(retornos)
      media = SUM(retornos) / REAL(n, dp)
      retorno_diario_livre = taxa_livre / DIAS_ANO
      downside_var = 0.0_dp
      DO i = 1, n
         IF (retornos(i) < retorno_diario_livre) &
            downside_var = downside_var + (retornos(i) - retorno_diario_livre)**2
      END DO
      downside_var = downside_var / REAL(n, dp)
      sortino = (media * DIAS_ANO - taxa_livre) / (SQRT(downside_var * DIAS_ANO) + 1.0D-15)
   END FUNCTION razao_sortino

   !----------------------------------------------------------------------------
   ! Máximo Drawdown (maior queda do pico ao vale)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION maximo_drawdown(precos) RESULT(mdd)
      REAL(dp), INTENT(IN) :: precos(:)
      REAL(dp) :: pico_atual, drawdown_atual
      INTEGER :: i
      pico_atual = precos(1)
      mdd = 0.0_dp
      DO i = 2, SIZE(precos)
         IF (precos(i) > pico_atual) pico_atual = precos(i)
         IF (pico_atual > 0.0_dp) THEN
            drawdown_atual = (pico_atual - precos(i)) / pico_atual
            IF (drawdown_atual > mdd) mdd = drawdown_atual
         END IF
      END DO
   END FUNCTION maximo_drawdown

   !----------------------------------------------------------------------------
   ! Beta de mercado (sensibilidade do ativo ao índice Ibovespa)
   ! beta = Cov(r_ativo, r_ibov) / Var(r_ibov)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION beta_mercado(ret_ativo, ret_ibov) RESULT(beta)
      REAL(dp), INTENT(IN) :: ret_ativo(:), ret_ibov(:)
      REAL(dp) :: media_a, media_b, covariancia, variancia_b
      INTEGER :: n, i
      n = MIN(SIZE(ret_ativo), SIZE(ret_ibov))
      media_a = SUM(ret_ativo(1:n)) / REAL(n, dp)
      media_b = SUM(ret_ibov(1:n))  / REAL(n, dp)
      covariancia = 0.0_dp;  variancia_b = 0.0_dp
      DO i = 1, n
         covariancia  = covariancia  + (ret_ativo(i) - media_a) * (ret_ibov(i) - media_b)
         variancia_b  = variancia_b  + (ret_ibov(i) - media_b)**2
      END DO
      beta = covariancia / (variancia_b + 1.0D-15)
   END FUNCTION beta_mercado

   !----------------------------------------------------------------------------
   ! Alfa de Jensen (retorno anormal em relação ao CAPM)
   ! alpha = retorno_ativo - [r_livre + beta*(r_ibov - r_livre)]
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION alfa_jensen(ret_medio_ativo, beta, ret_ibov_medio, &
      taxa_livre) RESULT(alfa)
      REAL(dp), INTENT(IN) :: ret_medio_ativo, beta, ret_ibov_medio, taxa_livre
      alfa = ret_medio_ativo - (taxa_livre + beta * (ret_ibov_medio - taxa_livre))
   END FUNCTION alfa_jensen

END MODULE mod_volatilidade
