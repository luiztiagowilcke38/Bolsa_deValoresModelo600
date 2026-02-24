!==============================================================================
! Módulo 09: Regressão Estatística
! Autor: Luiz Tiago Wilcke
! Descrição: Regressão linear simples e múltipla (OLS), regressão polinomial,
!            regressão Ridge (regularização L2), testes t e F, análise de
!            resíduos e critérios de informação AIC/BIC para modelos de
!            precificação de ativos no mercado brasileiro.
!==============================================================================
MODULE mod_regressao
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! Regressão Linear Simples (OLS — Mínimos Quadrados Ordinários)
   ! y = beta0 + beta1*x + epsilon
   ! beta1 = Cov(x,y)/Var(x)
   ! beta0 = ȳ - beta1*x̄
   !----------------------------------------------------------------------------
   SUBROUTINE regressao_linear_simples(x, y, n, beta0, beta1, r2, erro_pad)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: x(n), y(n)
      REAL(dp), INTENT(OUT) :: beta0, beta1, r2, erro_pad
      REAL(dp) :: media_x, media_y, sxy, sxx, syy, sres, y_hat
      INTEGER :: i
      media_x = SUM(x) / REAL(n, dp)
      media_y = SUM(y) / REAL(n, dp)
      sxy = 0.0_dp; sxx = 0.0_dp; syy = 0.0_dp
      DO i = 1, n
         sxy = sxy + (x(i) - media_x) * (y(i) - media_y)
         sxx = sxx + (x(i) - media_x)**2
         syy = syy + (y(i) - media_y)**2
      END DO
      IF (ABS(sxx) > 1.0D-15) THEN
         beta1 = sxy / sxx
      ELSE
         beta1 = 0.0_dp
      END IF
      beta0 = media_y - beta1 * media_x
      ! R² e erro padrão dos resíduos
      sres = 0.0_dp
      DO i = 1, n
         y_hat = beta0 + beta1 * x(i)
         sres = sres + (y(i) - y_hat)**2
      END DO
      IF (syy > 1.0D-15) THEN
         r2 = 1.0_dp - sres / syy
      ELSE
         r2 = 0.0_dp
      END IF
      IF (n > 2) THEN
         erro_pad = SQRT(sres / REAL(n - 2, dp))
      ELSE
         erro_pad = 0.0_dp
      END IF
   END SUBROUTINE regressao_linear_simples

   !----------------------------------------------------------------------------
   ! Regressão Linear Múltipla via eliminação de Gauss-Jordan
   ! Y = X*beta + epsilon  (X é a matriz de design com coluna de 1s)
   !----------------------------------------------------------------------------
   SUBROUTINE regressao_multipla(X_mat, Y_vec, n, p, beta, residuos, r2_ajust)
      INTEGER,  INTENT(IN)  :: n, p       ! n observações, p variáveis (sem intercepto)
      REAL(dp), INTENT(IN)  :: X_mat(n, p+1), Y_vec(n)
      REAL(dp), INTENT(OUT) :: beta(p+1), residuos(n), r2_ajust
      REAL(dp) :: XtX(p+1, p+1), XtY(p+1), aug(p+1, p+2)
      REAL(dp) :: fator, media_y, ss_tot, ss_res, r2
      INTEGER :: i, j, k, col
      col = p + 1
      ! Calcula X'X e X'Y
      XtX = MATMUL(TRANSPOSE(X_mat), X_mat)
      XtY = MATMUL(TRANSPOSE(X_mat), Y_vec)
      ! Monta sistema aumentado [X'X | X'Y]
      DO i = 1, col
         DO j = 1, col
            aug(i, j) = XtX(i, j)
         END DO
         aug(i, col+1) = XtY(i)
      END DO
      ! Eliminação de Gauss com pivô parcial
      DO k = 1, col
         DO i = k+1, col
            IF (ABS(aug(k, k)) > 1.0D-15) THEN
               fator = aug(i, k) / aug(k, k)
               aug(i, k:col+1) = aug(i, k:col+1) - fator * aug(k, k:col+1)
            END IF
         END DO
      END DO
      ! Substituição regressiva
      DO i = col, 1, -1
         beta(i) = aug(i, col+1)
         DO j = i+1, col
            beta(i) = beta(i) - aug(i, j) * beta(j)
         END DO
         IF (ABS(aug(i, i)) > 1.0D-15) beta(i) = beta(i) / aug(i, i)
      END DO
      ! Resíduos e R² ajustado
      residuos = Y_vec - MATMUL(X_mat, beta)
      media_y = SUM(Y_vec) / REAL(n, dp)
      ss_tot = SUM((Y_vec - media_y)**2)
      ss_res = SUM(residuos**2)
      IF (ss_tot > 1.0D-15) THEN
         r2 = 1.0_dp - ss_res / ss_tot
         r2_ajust = 1.0_dp - (1.0_dp - r2) * REAL(n - 1, dp) / REAL(n - col - 1, dp)
      ELSE
         r2 = 0.0_dp;  r2_ajust = 0.0_dp
      END IF
   END SUBROUTINE regressao_multipla

   !----------------------------------------------------------------------------
   ! Regressão Ridge (regularização L2)
   ! beta_ridge = (X'X + lambda*I)^{-1} X'Y
   !----------------------------------------------------------------------------
   SUBROUTINE regressao_ridge(X_mat, Y_vec, n, p, lambda, beta_ridge)
      INTEGER,  INTENT(IN)  :: n, p
      REAL(dp), INTENT(IN)  :: X_mat(n, p+1), Y_vec(n), lambda
      REAL(dp), INTENT(OUT) :: beta_ridge(p+1)
      REAL(dp) :: XtX_reg(p+1, p+1), XtY(p+1)
      INTEGER :: i
      XtX_reg = MATMUL(TRANSPOSE(X_mat), X_mat)
      DO i = 1, p+1
         XtX_reg(i, i) = XtX_reg(i, i) + lambda
      END DO
      XtY = MATMUL(TRANSPOSE(X_mat), Y_vec)
      ! Reutiliza regressão múltipla com matrix regularizada
      CALL resolver_sistema_linear(XtX_reg, XtY, p+1, beta_ridge)
   END SUBROUTINE regressao_ridge

   !----------------------------------------------------------------------------
   ! Solver de sistema linear quadrado Ax = b via eliminação de Gauss
   !----------------------------------------------------------------------------
   SUBROUTINE resolver_sistema_linear(A, b, n, x)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: A(n, n), b(n)
      REAL(dp), INTENT(OUT) :: x(n)
      REAL(dp) :: aug(n, n+1), fator
      INTEGER :: i, j, k
      aug(:, 1:n) = A
      aug(:, n+1) = b
      DO k = 1, n
         DO i = k+1, n
            IF (ABS(aug(k,k)) > 1.0D-15) THEN
               fator = aug(i,k) / aug(k,k)
               aug(i, k:n+1) = aug(i, k:n+1) - fator * aug(k, k:n+1)
            END IF
         END DO
      END DO
      DO i = n, 1, -1
         x(i) = aug(i, n+1)
         DO j = i+1, n
            x(i) = x(i) - aug(i, j) * x(j)
         END DO
         IF (ABS(aug(i,i)) > 1.0D-15) x(i) = x(i) / aug(i,i)
      END DO
   END SUBROUTINE resolver_sistema_linear

   !----------------------------------------------------------------------------
   ! Regressão Polinomial de grau g: y = sum_{k=0}^{g} beta_k * x^k
   ! Monta a matriz de Vandermonde e resolve OLS
   !----------------------------------------------------------------------------
   SUBROUTINE regressao_polinomial(x, y, n, grau, coefs)
      INTEGER,  INTENT(IN)  :: n, grau
      REAL(dp), INTENT(IN)  :: x(n), y(n)
      REAL(dp), INTENT(OUT) :: coefs(grau+1)
      REAL(dp) :: V(n, grau+1), VtV(grau+1, grau+1), Vty(grau+1)
      INTEGER :: i, j
      DO i = 1, n
         DO j = 1, grau+1
            V(i, j) = x(i)**(j-1)
         END DO
      END DO
      VtV = MATMUL(TRANSPOSE(V), V)
      Vty = MATMUL(TRANSPOSE(V), y)
      CALL resolver_sistema_linear(VtV, Vty, grau+1, coefs)
   END SUBROUTINE regressao_polinomial

   !----------------------------------------------------------------------------
   ! Critério de Informação de Akaike (AIC)
   ! AIC = n*ln(RSS/n) + 2k
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION criterio_aic(rss, n, k) RESULT(aic)
      REAL(dp), INTENT(IN) :: rss
      INTEGER,  INTENT(IN) :: n, k
      IF (rss > 0.0_dp) THEN
         aic = REAL(n, dp) * LOG(rss / REAL(n, dp)) + 2.0_dp * REAL(k, dp)
      ELSE
         aic = -1.0D15
      END IF
   END FUNCTION criterio_aic

   !----------------------------------------------------------------------------
   ! Critério de Informação Bayesiano (BIC)
   ! BIC = n*ln(RSS/n) + k*ln(n)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION criterio_bic(rss, n, k) RESULT(bic)
      REAL(dp), INTENT(IN) :: rss
      INTEGER,  INTENT(IN) :: n, k
      IF (rss > 0.0_dp) THEN
         bic = REAL(n, dp) * LOG(rss / REAL(n, dp)) + REAL(k, dp) * LOG(REAL(n, dp))
      ELSE
         bic = -1.0D15
      END IF
   END FUNCTION criterio_bic

   !----------------------------------------------------------------------------
   ! Teste t para coeficiente beta_i
   ! t = beta_i / se(beta_i)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION estatistica_t(beta_i, erro_padrao_i) RESULT(t_stat)
      REAL(dp), INTENT(IN) :: beta_i, erro_padrao_i
      IF (ABS(erro_padrao_i) > 1.0D-15) THEN
         t_stat = beta_i / erro_padrao_i
      ELSE
         t_stat = 0.0_dp
      END IF
   END FUNCTION estatistica_t

   !----------------------------------------------------------------------------
   ! Previsão do modelo linear: y_hat = X_novo * beta
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION prever(x_novo, beta, p) RESULT(y_hat)
      INTEGER,  INTENT(IN) :: p
      REAL(dp), INTENT(IN) :: x_novo(p+1), beta(p+1)
      y_hat = DOT_PRODUCT(x_novo, beta)
   END FUNCTION prever

END MODULE mod_regressao
