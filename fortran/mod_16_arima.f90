!==============================================================================
! Módulo 16: Modelo ARIMA e Previsão de Séries Temporais
! Autor: Luiz Tiago Wilcke
! Descrição: Modelos AR, MA, ARMA e ARIMA completos com estimação por
!            mínimos quadrados condicionais. Inclui diferenciação, inversão,
!            critérios de seleção de ordem e previsão para retornos da B3.
!==============================================================================
MODULE mod_arima
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! Diferenciação de ordem d: nabla^d(X_t)
   ! nabla X_t = X_t - X_{t-1}
   !----------------------------------------------------------------------------
   SUBROUTINE diferenciar(serie, d, serie_diff, n_diff)
      REAL(dp), INTENT(IN)  :: serie(:)
      INTEGER,  INTENT(IN)  :: d
      REAL(dp), INTENT(OUT) :: serie_diff(:)
      INTEGER,  INTENT(OUT) :: n_diff
      REAL(dp), ALLOCATABLE :: temp(:)
      INTEGER :: n, i, k
      n = SIZE(serie)
      ALLOCATE(temp(n))
      temp = serie
      DO k = 1, d
         DO i = 2, n - k + 1
            temp(i-1) = temp(i) - temp(i-1)
         END DO
      END DO
      n_diff = n - d
      serie_diff(1:n_diff) = temp(1:n_diff)
      DEALLOCATE(temp)
   END SUBROUTINE diferenciar

   !----------------------------------------------------------------------------
   ! Modelo AR(p): X_t = c + phi_1*X_{t-1} + ... + phi_p*X_{t-p} + eps_t
   ! Estimação por OLS (Yule-Walker equações)
   !----------------------------------------------------------------------------
   SUBROUTINE estimar_ar(serie, p, phi, constante, residuos)
      REAL(dp), INTENT(IN)  :: serie(:)
      INTEGER,  INTENT(IN)  :: p
      REAL(dp), INTENT(OUT) :: phi(p), constante, residuos(:)
      REAL(dp) :: XtX(p+1, p+1), XtY(p+1), x_linha(p+1), aug(p+1, p+2), fator
      INTEGER :: n, t, i, j, k
      n = SIZE(serie)
      XtX = 0.0_dp;  XtY = 0.0_dp
      DO t = p+1, n
         x_linha(1) = 1.0_dp
         DO i = 1, p
            x_linha(i+1) = serie(t-i)
         END DO
         DO i = 1, p+1
            DO j = 1, p+1
               XtX(i,j) = XtX(i,j) + x_linha(i) * x_linha(j)
            END DO
            XtY(i) = XtY(i) + x_linha(i) * serie(t)
         END DO
      END DO
      ! Solve via Gauss com aumentada [XtX | XtY]
      DO i = 1, p+1
         DO j = 1, p+1
            aug(i,j) = XtX(i,j)
         END DO
         aug(i, p+2) = XtY(i)
      END DO
      DO k = 1, p+1
         IF (ABS(aug(k,k)) < 1.0D-15) aug(k,k) = 1.0D-12
         DO i = k+1, p+1
            fator = aug(i,k) / aug(k,k)
            aug(i, k:p+2) = aug(i, k:p+2) - fator * aug(k, k:p+2)
         END DO
      END DO
      DO i = p+1, 1, -1
         phi(MAX(i-1,1)) = aug(i, p+2)   ! temp reutilizado
         DO j = i+1, p+1
            IF (j-1 <= p) phi(MAX(i-1,1)) = phi(MAX(i-1,1)) - aug(i,j) * phi(MAX(j-1,1))
         END DO
         IF (ABS(aug(i,i)) > 1.0D-15) phi(MAX(i-1,1)) = phi(MAX(i-1,1)) / aug(i,i)
      END DO
      constante = phi(1)
      IF (p > 0) phi(1:p) = phi(2:MIN(p+1,p+1))
      ! Calcula resíduos
      DO t = p+1, n
         residuos(t) = serie(t) - constante
         DO i = 1, p
            residuos(t) = residuos(t) - phi(i) * serie(t-i)
         END DO
      END DO
   END SUBROUTINE estimar_ar

   !----------------------------------------------------------------------------
   ! Previsão h passos à frente do modelo AR(p)
   ! X_{t+h} = c + phi_1*X_{t+h-1} + ... + phi_p*X_{t+h-p}
   !----------------------------------------------------------------------------
   SUBROUTINE previsao_ar(serie, n, phi, constante, p, h_passos, previsoes)
      INTEGER,  INTENT(IN)  :: n, p, h_passos
      REAL(dp), INTENT(IN)  :: serie(n), phi(p), constante
      REAL(dp), INTENT(OUT) :: previsoes(h_passos)
      REAL(dp) :: historico(n + h_passos)
      INTEGER :: h, i
      historico(1:n) = serie
      DO h = 1, h_passos
         historico(n+h) = constante
         DO i = 1, p
            IF (n + h - i >= 1) &
               historico(n+h) = historico(n+h) + phi(i) * historico(n+h-i)
         END DO
         previsoes(h) = historico(n+h)
      END DO
   END SUBROUTINE previsao_ar

   !----------------------------------------------------------------------------
   ! Filtragem do componente MA(q) a partir de resíduos estimados
   ! eps_t = X_t - c - phi*X_{t-1} - theta_1*eps_{t-1} - ... - theta_q*eps_{t-q}
   !----------------------------------------------------------------------------
   SUBROUTINE filtrar_ma(serie, constante, phi, p, theta, q, residuos)
      REAL(dp), INTENT(IN)  :: serie(:), constante, phi(:), theta(:)
      INTEGER,  INTENT(IN)  :: p, q
      REAL(dp), INTENT(OUT) :: residuos(:)
      INTEGER :: n, t, i
      n = SIZE(serie)
      residuos = 0.0_dp
      DO t = MAX(p,q)+1, n
         residuos(t) = serie(t) - constante
         DO i = 1, p
            IF (t-i >= 1) residuos(t) = residuos(t) - phi(i) * serie(t-i)
         END DO
         DO i = 1, q
            IF (t-i >= 1) residuos(t) = residuos(t) - theta(i) * residuos(t-i)
         END DO
      END DO
   END SUBROUTINE filtrar_ma

   !----------------------------------------------------------------------------
   ! AIC/BIC para modelo ARIMA
   ! AIC = 2k - 2*loglike
   ! BIC = k*ln(n) - 2*loglike
   !----------------------------------------------------------------------------
   SUBROUTINE criterios_arima(residuos, n, p, q, aic, bic)
      INTEGER,  INTENT(IN)  :: n, p, q
      REAL(dp), INTENT(IN)  :: residuos(n)
      REAL(dp), INTENT(OUT) :: aic, bic
      REAL(dp) :: sigma2, loglike
      INTEGER :: k
      k = p + q + 1   ! número de parâmetros (incluindo constante)
      sigma2  = SUM(residuos**2) / REAL(n, dp)
      IF (sigma2 > 1.0D-15) THEN
         loglike = -0.5_dp * REAL(n, dp) * LOG(sigma2)
      ELSE
         loglike = 0.0_dp
      END IF
      aic = 2.0_dp * REAL(k, dp)   - 2.0_dp * loglike
      bic = REAL(k, dp) * LOG(REAL(n, dp)) - 2.0_dp * loglike
   END SUBROUTINE criterios_arima

   !----------------------------------------------------------------------------
   ! Teste de portmanteau de Box-Pierce para resíduos de ARIMA
   ! Q = n * sum_{k=1}^{m} rho_k^2   ~ chi2(m-p-q)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION teste_box_pierce(residuos, n, m) RESULT(Q)
      INTEGER,  INTENT(IN) :: n, m
      REAL(dp), INTENT(IN) :: residuos(n)
      REAL(dp) :: media, var, cov, rho_k
      INTEGER :: k, i
      media = SUM(residuos) / REAL(n, dp)
      var   = SUM((residuos - media)**2) / REAL(n, dp)
      Q     = 0.0_dp
      DO k = 1, m
         cov = 0.0_dp
         DO i = 1, n - k
            cov = cov + (residuos(i) - media) * (residuos(i+k) - media)
         END DO
         cov = cov / REAL(n, dp)
         IF (ABS(var) > 1.0D-15) THEN
            rho_k = cov / var
            Q = Q + rho_k**2
         END IF
      END DO
      Q = REAL(n, dp) * Q
   END FUNCTION teste_box_pierce

   !----------------------------------------------------------------------------
   ! Integração (inversão da diferenciação) para retornar à escala original
   ! X_t = X_{t-1} + nabla(X_t)
   !----------------------------------------------------------------------------
   SUBROUTINE integrar_serie(serie_diff, valor_inicial, d, n, serie_orig)
      INTEGER,  INTENT(IN)  :: n, d
      REAL(dp), INTENT(IN)  :: serie_diff(n), valor_inicial
      REAL(dp), INTENT(OUT) :: serie_orig(n+d)
      INTEGER :: i, k
      serie_orig(1:d) = valor_inicial
      serie_orig(d+1:n+d) = serie_diff(1:n)
      DO k = 1, d
         DO i = 2, n + d
            serie_orig(i) = serie_orig(i) + serie_orig(i-1)
         END DO
      END DO
   END SUBROUTINE integrar_serie

END MODULE mod_arima
