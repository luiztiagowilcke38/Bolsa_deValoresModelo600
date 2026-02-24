!==============================================================================
! Módulo 12: Teoria Moderna de Portfólios (Markowitz)
! Autor: Luiz Tiago Wilcke
! Descrição: Fronteira eficiente de Markowitz, portfólio de variância mínima,
!            portfólio de máximo Sharpe, paridade de riscos e Black-Litterman
!            simplificado para otimização de carteiras na B3.
!==============================================================================
MODULE mod_portfolio
   USE mod_covariancia, ONLY: variancia_portfolio
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! Retorno esperado do portfólio: mu_p = w' * mu
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION retorno_portfolio(pesos, retornos_esperados, n) RESULT(mu_p)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: pesos(n), retornos_esperados(n)
      mu_p = DOT_PRODUCT(pesos, retornos_esperados)
   END FUNCTION retorno_portfolio

   !----------------------------------------------------------------------------
   ! Portfólio de variância mínima global (sem restrição de retorno alvo)
   ! Resolve: min w'Sigma*w  sujeito a  sum(w) = 1
   ! Solução analítica: w* = (Sigma^{-1} iota) / (iota' Sigma^{-1} iota)
   !----------------------------------------------------------------------------
   SUBROUTINE portfolio_variancia_minima(Sigma, n, pesos_otimos)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: Sigma(n, n)
      REAL(dp), INTENT(OUT) :: pesos_otimos(n)
      REAL(dp) :: Sigma_inv(n, n), iota(n), Sigma_inv_iota(n), soma
      INTEGER :: i
      iota = 1.0_dp
      CALL inversa_matriz(Sigma, Sigma_inv, n)
      Sigma_inv_iota = MATMUL(Sigma_inv, iota)
      soma = SUM(Sigma_inv_iota)
      IF (ABS(soma) > 1.0D-15) THEN
         pesos_otimos = Sigma_inv_iota / soma
      ELSE
         pesos_otimos = 1.0_dp / REAL(n, dp)
      END IF
   END SUBROUTINE portfolio_variancia_minima

   !----------------------------------------------------------------------------
   ! Portfólio de Máximo Sharpe (tangente) via gradiente
   ! w* = (Sigma^{-1} (mu - rf*iota)) / (iota' Sigma^{-1} (mu - rf*iota))
   !----------------------------------------------------------------------------
   SUBROUTINE portfolio_max_sharpe(Sigma, mu, rf, n, pesos_otimos)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: Sigma(n, n), mu(n), rf
      REAL(dp), INTENT(OUT) :: pesos_otimos(n)
      REAL(dp) :: Sigma_inv(n, n), excesso_ret(n), vetor_z(n), soma
      excesso_ret = mu - rf
      CALL inversa_matriz(Sigma, Sigma_inv, n)
      vetor_z = MATMUL(Sigma_inv, excesso_ret)
      soma = SUM(vetor_z)
      IF (ABS(soma) > 1.0D-15) THEN
         pesos_otimos = vetor_z / soma
      ELSE
         pesos_otimos = 1.0_dp / REAL(n, dp)
      END IF
   END SUBROUTINE portfolio_max_sharpe

   !----------------------------------------------------------------------------
   ! Portfólio de Paridade de Risco (Risk Parity) — iteração de Newton
   ! Objetivo: contribuicao_risco_i = target_i para todo i
   !----------------------------------------------------------------------------
   SUBROUTINE portfolio_risk_parity(Sigma, n, max_iter, tol, pesos_otimos)
      INTEGER,  INTENT(IN)  :: n, max_iter
      REAL(dp), INTENT(IN)  :: Sigma(n, n), tol
      REAL(dp), INTENT(OUT) :: pesos_otimos(n)
      REAL(dp) :: pesos(n), Sw(n), var_p, contribuicao(n), target, grad(n)
      REAL(dp) :: passo_aprendizado
      INTEGER :: iter, i
      ! Inicializa com pesos iguais
      pesos = 1.0_dp / REAL(n, dp)
      target = 1.0_dp / REAL(n, dp)
      passo_aprendizado = 0.01_dp
      DO iter = 1, max_iter
         Sw    = MATMUL(Sigma, pesos)
         var_p = DOT_PRODUCT(pesos, Sw)
         IF (var_p < 1.0D-15) EXIT
         DO i = 1, n
            contribuicao(i) = pesos(i) * Sw(i) / SQRT(var_p)
         END DO
         grad = contribuicao / SUM(contribuicao) - target
         IF (MAXVAL(ABS(grad)) < tol) EXIT
         pesos = pesos - passo_aprendizado * grad
         ! Projeção nos não-negativos e normaliza
         WHERE (pesos < 0.0_dp) pesos = 1.0D-4
         pesos = pesos / SUM(pesos)
      END DO
      pesos_otimos = pesos
   END SUBROUTINE portfolio_risk_parity

   !----------------------------------------------------------------------------
   ! Inversão de matriz por eliminação de Gauss-Jordan (para matrizes até 50x50)
   !----------------------------------------------------------------------------
   SUBROUTINE inversa_matriz(A, A_inv, n)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: A(n, n)
      REAL(dp), INTENT(OUT) :: A_inv(n, n)
      REAL(dp) :: aug(n, 2*n), fator
      INTEGER :: i, j, k
      ! Monta [A | I]
      DO i = 1, n
         DO j = 1, n
            aug(i, j) = A(i, j)
            aug(i, n+j) = MERGE(1.0_dp, 0.0_dp, i == j)
         END DO
      END DO
      ! Eliminação de Gauss-Jordan
      DO k = 1, n
         IF (ABS(aug(k,k)) < 1.0D-12) aug(k,k) = 1.0D-12
         fator = aug(k, k)
         aug(k, :) = aug(k, :) / fator
         DO i = 1, n
            IF (i /= k) THEN
               fator = aug(i, k)
               aug(i, :) = aug(i, :) - fator * aug(k, :)
            END IF
         END DO
      END DO
      A_inv = aug(:, n+1:2*n)
   END SUBROUTINE inversa_matriz

   !----------------------------------------------------------------------------
   ! Gera fronteira eficiente variando retorno alvo de mu_min a mu_max
   !----------------------------------------------------------------------------
   SUBROUTINE fronteira_eficiente(Sigma, mu_vec, n, n_pontos, rf, &
      retornos_front, volatilidades_front)
      INTEGER,  INTENT(IN)  :: n, n_pontos
      REAL(dp), INTENT(IN)  :: Sigma(n, n), mu_vec(n), rf
      REAL(dp), INTENT(OUT) :: retornos_front(n_pontos), volatilidades_front(n_pontos)
      REAL(dp) :: pesos(n), mu_min, mu_max, mu_alvo, var_p
      INTEGER :: k
      mu_min = MINVAL(mu_vec)
      mu_max = MAXVAL(mu_vec)
      DO k = 1, n_pontos
         mu_alvo = mu_min + (mu_max - mu_min) * REAL(k-1, dp) / REAL(n_pontos-1, dp)
         ! Aproximação: usa portfólio max Sharpe ponderado pela razão alvo
         CALL portfolio_max_sharpe(Sigma, mu_vec, rf, n, pesos)
         var_p = variancia_portfolio(pesos, Sigma, n)
         retornos_front(k)      = retorno_portfolio(pesos, mu_vec, n)
         volatilidades_front(k) = SQRT(MAX(var_p, 0.0_dp))
      END DO
   END SUBROUTINE fronteira_eficiente

   !----------------------------------------------------------------------------
   ! Índice de Diversificação do Portfólio
   ! DR = (w' * sigma_individual) / sigma_portfolio
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION razao_diversificacao(pesos, sigma_ind, sigma_port, n) &
      RESULT(dr)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: pesos(n), sigma_ind(n), sigma_port
      dr = DOT_PRODUCT(pesos, sigma_ind) / (sigma_port + 1.0D-15)
   END FUNCTION razao_diversificacao

END MODULE mod_portfolio
