!==============================================================================
! Módulo 11: Covariância e Matrizes de Risco
! Autor: Luiz Tiago Wilcke
! Descrição: Estimação de matrizes de covariância (amostral, EWMA, Ledoit-Wolf),
!            decomposição espectral, decomposição de Cholesky e aplicações
!            em gestão de risco de portfólios para o mercado da B3.
!==============================================================================
MODULE mod_covariancia
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! Matriz de covariância amostral S = (1/(n-1)) * (X - X̄)'(X - X̄)
   !----------------------------------------------------------------------------
   SUBROUTINE covariancia_amostral(retornos, n_obs, n_ativos, S)
      INTEGER,  INTENT(IN)  :: n_obs, n_ativos
      REAL(dp), INTENT(IN)  :: retornos(n_obs, n_ativos)
      REAL(dp), INTENT(OUT) :: S(n_ativos, n_ativos)
      REAL(dp) :: medias(n_ativos), X_cent(n_obs, n_ativos)
      INTEGER :: i, j, k
      DO j = 1, n_ativos
         medias(j) = SUM(retornos(:, j)) / REAL(n_obs, dp)
         DO i = 1, n_obs
            X_cent(i, j) = retornos(i, j) - medias(j)
         END DO
      END DO
      DO j = 1, n_ativos
         DO k = j, n_ativos
            S(j, k) = DOT_PRODUCT(X_cent(:, j), X_cent(:, k)) / REAL(n_obs - 1, dp)
            S(k, j) = S(j, k)
         END DO
      END DO
   END SUBROUTINE covariancia_amostral

   !----------------------------------------------------------------------------
   ! Matriz de covariância EWMA (RiskMetrics)
   ! Sigma_t = lambda*Sigma_{t-1} + (1-lambda)*r_{t-1}*r_{t-1}'
   !----------------------------------------------------------------------------
   SUBROUTINE covariancia_ewma(retornos, n_obs, n_ativos, lambda, Sigma_final)
      INTEGER,  INTENT(IN)  :: n_obs, n_ativos
      REAL(dp), INTENT(IN)  :: retornos(n_obs, n_ativos), lambda
      REAL(dp), INTENT(OUT) :: Sigma_final(n_ativos, n_ativos)
      REAL(dp) :: Sigma(n_ativos, n_ativos), r_outer(n_ativos, n_ativos)
      INTEGER :: i, j, k
      ! Inicializa com covariância dos primeiros 20 períodos
      CALL covariancia_amostral(retornos(1:MIN(20,n_obs),:), MIN(20,n_obs), &
         n_ativos, Sigma)
      DO i = 21, n_obs
         DO j = 1, n_ativos
            DO k = 1, n_ativos
               r_outer(j, k) = retornos(i-1, j) * retornos(i-1, k)
            END DO
         END DO
         Sigma = lambda * Sigma + (1.0_dp - lambda) * r_outer
      END DO
      Sigma_final = Sigma
   END SUBROUTINE covariancia_ewma

   !----------------------------------------------------------------------------
   ! Encolhimento de Ledoit-Wolf (simplificado — target é identidade escalada)
   ! S_LW = (1-delta)*S + delta*mu_bar*I
   ! onde mu_bar = trace(S)/p e delta é coeficiente de encolhimento
   !----------------------------------------------------------------------------
   SUBROUTINE covariancia_ledoit_wolf(S, n_ativos, delta, S_lw)
      INTEGER,  INTENT(IN)  :: n_ativos
      REAL(dp), INTENT(IN)  :: S(n_ativos, n_ativos), delta
      REAL(dp), INTENT(OUT) :: S_lw(n_ativos, n_ativos)
      REAL(dp) :: mu_bar, traco
      INTEGER :: i, j
      traco = 0.0_dp
      DO i = 1, n_ativos
         traco = traco + S(i, i)
      END DO
      mu_bar = traco / REAL(n_ativos, dp)
      DO i = 1, n_ativos
         DO j = 1, n_ativos
            IF (i == j) THEN
               S_lw(i, j) = (1.0_dp - delta) * S(i, j) + delta * mu_bar
            ELSE
               S_lw(i, j) = (1.0_dp - delta) * S(i, j)
            END IF
         END DO
      END DO
   END SUBROUTINE covariancia_ledoit_wolf

   !----------------------------------------------------------------------------
   ! Decomposição de Cholesky: S = L * L'  (para S definida positiva)
   ! Usada para geração de variáveis correlacionadas em Monte Carlo
   !----------------------------------------------------------------------------
   SUBROUTINE cholesky(S, L, n)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: S(n, n)
      REAL(dp), INTENT(OUT) :: L(n, n)
      REAL(dp) :: soma
      INTEGER :: i, j, k
      L = 0.0_dp
      DO i = 1, n
         DO j = 1, i
            soma = S(i, j)
            DO k = 1, j - 1
               soma = soma - L(i, k) * L(j, k)
            END DO
            IF (i == j) THEN
               IF (soma > 0.0_dp) THEN
                  L(i, i) = SQRT(soma)
               ELSE
                  L(i, i) = 1.0D-8
               END IF
            ELSE
               IF (ABS(L(j, j)) > 1.0D-15) THEN
                  L(i, j) = soma / L(j, j)
               END IF
            END IF
         END DO
      END DO
   END SUBROUTINE cholesky

   !----------------------------------------------------------------------------
   ! Variância do portfólio: sigma_p^2 = w' * Sigma * w
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION variancia_portfolio(pesos, Sigma, n) RESULT(var_p)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: pesos(n), Sigma(n, n)
      REAL(dp) :: Sw(n)
      Sw = MATMUL(Sigma, pesos)
      var_p = DOT_PRODUCT(pesos, Sw)
   END FUNCTION variancia_portfolio

   !----------------------------------------------------------------------------
   ! Contribuição marginal ao risco de cada ativo no portfólio
   ! MCR_i = (Sigma * w)_i / sigma_p
   !----------------------------------------------------------------------------
   SUBROUTINE contribuicao_marginal_risco(pesos, Sigma, n, mcr)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: pesos(n), Sigma(n, n)
      REAL(dp), INTENT(OUT) :: mcr(n)
      REAL(dp) :: Sw(n), var_p
      Sw    = MATMUL(Sigma, pesos)
      var_p = DOT_PRODUCT(pesos, Sw)
      IF (var_p > 1.0D-15) THEN
         mcr = Sw / SQRT(var_p)
      ELSE
         mcr = 0.0_dp
      END IF
   END SUBROUTINE contribuicao_marginal_risco

   !----------------------------------------------------------------------------
   ! Diversificação efetiva: número efetivo de ativos (ENP)
   ! ENP = 1 / sum(w_i^2)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION numero_efetivo_ativos(pesos, n) RESULT(ENP)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: pesos(n)
      ENP = 1.0_dp / (DOT_PRODUCT(pesos, pesos) + 1.0D-15)
   END FUNCTION numero_efetivo_ativos

   !----------------------------------------------------------------------------
   ! Determinante de matriz (expansão de Laplace, apenas até 4x4 para eficiência)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION determinante(A, n) RESULT(det)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: A(n, n)
      REAL(dp) :: A_copia(n, n), fator
      INTEGER :: i, j, k
      A_copia = A
      det = 1.0_dp
      DO k = 1, n
         DO i = k+1, n
            IF (ABS(A_copia(k,k)) > 1.0D-15) THEN
               fator = A_copia(i, k) / A_copia(k, k)
               A_copia(i, k:n) = A_copia(i, k:n) - fator * A_copia(k, k:n)
            END IF
         END DO
         det = det * A_copia(k, k)
      END DO
   END FUNCTION determinante

END MODULE mod_covariancia
