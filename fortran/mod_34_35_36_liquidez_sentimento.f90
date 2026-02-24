!==============================================================================
! Módulo 34-35-36: Liquidez, Sentimento e Análise Fatorial
! Autor: Luiz Tiago Wilcke
! Descrição: Métricas de liquidez de mercado (Amihud, Roll, Kyle),
!            indicadores de sentimento do investidor (IHH, put/call ratio,
!            índice de medo), Análise de Componentes Principais (PCA) e
!            modelo fatorial (Fama-French 3 e 5 fatores) para a B3.
!==============================================================================
MODULE mod_liquidez_sentimento_fatorial
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !============ LIQUIDEZ ======================================================

   !----------------------------------------------------------------------------
   ! Illiquidez de Amihud (2002)
   ! ILLIQ_t = |r_t| / (Volume_t * Preco_t)  * 10^6
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION illiquidez_amihud(retorno, volume, preco) RESULT(ILLIQ)
      REAL(dp), INTENT(IN) :: retorno, volume, preco
      REAL(dp) :: volume_R
      volume_R = volume * preco   ! volume em R$
      IF (volume_R > 1.0D-10) THEN
         ILLIQ = ABS(retorno) / volume_R * 1.0D6
      ELSE
         ILLIQ = 0.0_dp
      END IF
   END FUNCTION illiquidez_amihud

   !----------------------------------------------------------------------------
   ! Spread efetivo de Roll (1984) — estimado pela autocovariância dos retornos
   ! Spread = 2 * sqrt(max(-cov(r_t, r_{t-1}), 0))
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION spread_roll(retornos, n) RESULT(spread)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: retornos(n)
      REAL(dp) :: media, cov_lag1
      INTEGER :: i
      media = SUM(retornos) / REAL(n, dp)
      cov_lag1 = 0.0_dp
      DO i = 2, n
         cov_lag1 = cov_lag1 + (retornos(i) - media) * (retornos(i-1) - media)
      END DO
      cov_lag1 = cov_lag1 / REAL(n - 1, dp)
      spread   = 2.0_dp * SQRT(MAX(-cov_lag1, 0.0_dp))
   END FUNCTION spread_roll

   !----------------------------------------------------------------------------
   ! Impacto de mercado de Kyle (lambda de Kyle)
   ! Estima lambda por regressão: delta_P_t = lambda * volume_t + epsilon_t
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION lambda_kyle(delta_preco, volume, n) RESULT(lambda)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: delta_preco(n), volume(n)
      REAL(dp) :: sum_vol2, sum_pv, m_v, m_p
      m_v = SUM(volume) / REAL(n, dp)
      m_p = SUM(delta_preco) / REAL(n, dp)
      sum_vol2 = SUM((volume - m_v)**2)
      sum_pv   = DOT_PRODUCT(volume - m_v, delta_preco - m_p)
      IF (sum_vol2 > 1.0D-15) THEN
         lambda = sum_pv / sum_vol2
      ELSE
         lambda = 0.0_dp
      END IF
   END FUNCTION lambda_kyle

   !============ SENTIMENTO ====================================================

   !----------------------------------------------------------------------------
   ! Razão Put/Call (indicador de sentimento)
   ! > 1.2: pessimismo extremo (sinal de compra contrário)
   ! < 0.7: otimismo extremo (sinal de venda contrário)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION razao_put_call(vol_puts, vol_calls) RESULT(rpc)
      REAL(dp), INTENT(IN) :: vol_puts, vol_calls
      rpc = vol_puts / (vol_calls + 1.0D-10)
   END FUNCTION razao_put_call

   !----------------------------------------------------------------------------
   ! Índice de Força do Mercado (IFM) como proxy de sentimento
   ! IFM = (ações que sobem - ações que descem) / total de ações
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION indice_forca_mercado(n_alta, n_baixa, n_neutro) RESULT(IFM)
      INTEGER, INTENT(IN) :: n_alta, n_baixa, n_neutro
      REAL(dp) :: n_total
      n_total = REAL(n_alta + n_baixa + n_neutro, dp)
      IF (n_total > 0.0_dp) THEN
         IFM = REAL(n_alta - n_baixa, dp) / n_total
      ELSE
         IFM = 0.0_dp
      END IF
   END FUNCTION indice_forca_mercado

   !============ ANÁLISE FATORIAL =============================================

   !----------------------------------------------------------------------------
   ! Modelo Fama-French 3 Fatores
   ! ri - rf = alpha_i + beta_i*(rm - rf) + si*SMB + hi*HML + epsilon_i
   ! SMB = retorno Small-Minus-Big (tamanho)
   ! HML = retorno High-Minus-Low (valor/book-to-market)
   !----------------------------------------------------------------------------
   SUBROUTINE fama_french_3f(r_ativo, r_mkt, SMB, HML, n, &
      alpha, beta_mkt, beta_smb, beta_hml, r2)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: r_ativo(n), r_mkt(n), SMB(n), HML(n)
      REAL(dp), INTENT(OUT) :: alpha, beta_mkt, beta_smb, beta_hml, r2
      REAL(dp) :: X(n, 4), Y(n), XtX(4,4), XtY(4), beta(4), aug(4,5), fator
      REAL(dp) :: residuos(n), media_y, ss_tot, ss_res
      INTEGER :: i, j, k
      ! Monta matriz design [1, rm-rf, SMB, HML]
      DO i = 1, n
         X(i,1) = 1.0_dp
         X(i,2) = r_mkt(i)
         X(i,3) = SMB(i)
         X(i,4) = HML(i)
         Y(i)   = r_ativo(i)
      END DO
      XtX = MATMUL(TRANSPOSE(X), X)
      XtY = MATMUL(TRANSPOSE(X), Y)
      ! Aumentada
      DO i = 1, 4
         DO j = 1, 4; aug(i,j) = XtX(i,j); END DO
         aug(i,5) = XtY(i)
      END DO
      ! Gauss
      DO k = 1, 4
         DO i = k+1, 4
            fator = aug(i,k) / (aug(k,k) + 1.0D-15)
            aug(i,k:5) = aug(i,k:5) - fator * aug(k,k:5)
         END DO
      END DO
      DO i = 4, 1, -1
         beta(i) = aug(i,5)
         DO j = i+1, 4
            beta(i) = beta(i) - aug(i,j) * beta(j)
         END DO
         beta(i) = beta(i) / (aug(i,i) + 1.0D-15)
      END DO
      alpha    = beta(1);  beta_mkt = beta(2)
      beta_smb = beta(3);  beta_hml = beta(4)
      ! R²
      residuos = Y - MATMUL(X, beta)
      media_y  = SUM(Y) / REAL(n, dp)
      ss_tot   = SUM((Y - media_y)**2)
      ss_res   = SUM(residuos**2)
      r2 = MERGE(1.0_dp - ss_res/ss_tot, 0.0_dp, ss_tot > 1.0D-15)
   END SUBROUTINE fama_french_3f

   !----------------------------------------------------------------------------
   ! PCA simplificada — encontra primeiro componente principal via iteração
   ! de potência (approx): v1 = eigenvector correspondente ao maior autovalor
   !----------------------------------------------------------------------------
   SUBROUTINE pca_primeiro_componente(S, n, v1, lambda1)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: S(n, n)
      REAL(dp), INTENT(OUT) :: v1(n), lambda1
      REAL(dp) :: v(n), v_novo(n), norma
      INTEGER :: iter, i
      ! Inicialização aleatória
      v = 1.0_dp / SQRT(REAL(n, dp))
      DO iter = 1, 200
         v_novo = MATMUL(S, v)
         norma  = SQRT(DOT_PRODUCT(v_novo, v_novo))
         IF (norma > 1.0D-15) THEN
            v = v_novo / norma
         ELSE
            EXIT
         END IF
      END DO
      v1     = v
      lambda1 = DOT_PRODUCT(v, MATMUL(S, v))
   END SUBROUTINE pca_primeiro_componente

END MODULE mod_liquidez_sentimento_fatorial
