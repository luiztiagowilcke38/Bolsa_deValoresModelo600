!==============================================================================
! Módulo 15: Modelo GARCH (Heterocedasticidade Condicional)
! Autor: Luiz Tiago Wilcke
! Descrição: Estimação e previsão de GARCH(p,q), EGARCH (assimétrico),
!            GJR-GARCH (efeito leverage), GARCH-M e cálculo de volatilidade
!            condicional para ativos brasileiros na B3.
!==============================================================================
MODULE mod_garch
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! GARCH(1,1): sigma_t^2 = omega + alpha*epsilon_{t-1}^2 + beta*sigma_{t-1}^2
   ! Restrição: omega > 0, alpha >= 0, beta >= 0, alpha+beta < 1
   !----------------------------------------------------------------------------
   SUBROUTINE garch11(retornos, omega, alpha, beta, sigma2_out)
      REAL(dp), INTENT(IN)  :: retornos(:), omega, alpha, beta
      REAL(dp), INTENT(OUT) :: sigma2_out(:)
      REAL(dp) :: variancia_longo_prazo
      INTEGER :: n, t
      n = SIZE(retornos)
      ! Variância de longo prazo como valor inicial
      variancia_longo_prazo = SUM(retornos**2) / REAL(n, dp)
      sigma2_out(1) = variancia_longo_prazo
      DO t = 2, n
         sigma2_out(t) = omega + alpha * retornos(t-1)**2 + beta * sigma2_out(t-1)
      END DO
   END SUBROUTINE garch11

   !----------------------------------------------------------------------------
   ! GARCH(p,q) geral
   ! sigma_t^2 = omega + sum_{i=1}^{q} alpha_i*eps_{t-i}^2 +
   !                     sum_{j=1}^{p} beta_j*sigma_{t-j}^2
   !----------------------------------------------------------------------------
   SUBROUTINE garch_pq(retornos, omega, alphas, betas, p, q, sigma2_out)
      REAL(dp), INTENT(IN)  :: retornos(:), omega, alphas(q), betas(p)
      INTEGER,  INTENT(IN)  :: p, q
      REAL(dp), INTENT(OUT) :: sigma2_out(:)
      REAL(dp) :: var_inicializ, soma
      INTEGER :: n, t, i, inicio
      n = SIZE(retornos)
      var_inicializ = SUM(retornos**2) / REAL(n, dp)
      inicio = MAX(p, q) + 1
      DO t = 1, inicio - 1
         sigma2_out(t) = var_inicializ
      END DO
      DO t = inicio, n
         soma = omega
         DO i = 1, q
            IF (t - i >= 1) soma = soma + alphas(i) * retornos(t-i)**2
         END DO
         DO i = 1, p
            IF (t - i >= 1) soma = soma + betas(i) * sigma2_out(t-i)
         END DO
         sigma2_out(t) = MAX(soma, 1.0D-12)
      END DO
   END SUBROUTINE garch_pq

   !----------------------------------------------------------------------------
   ! EGARCH(1,1) — Nelson (1991): permite efeitos assimétricos (leverage)
   ! ln(sigma_t^2) = omega + alpha*[|z_{t-1}| - E|z|] + gamma*z_{t-1} + beta*ln(sigma_{t-1}^2)
   ! z_t = epsilon_t / sigma_t
   !----------------------------------------------------------------------------
   SUBROUTINE egarch11(retornos, omega, alpha, gamma_asim, beta, log_sigma2)
      REAL(dp), INTENT(IN)  :: retornos(:), omega, alpha, gamma_asim, beta
      REAL(dp), INTENT(OUT) :: log_sigma2(:)
      REAL(dp), PARAMETER :: E_ABS_Z = 0.7978845608_dp  ! sqrt(2/pi)
      REAL(dp) :: z_prev, ln_s2_prev
      INTEGER :: n, t
      n = SIZE(retornos)
      ln_s2_prev  = LOG(MAX(SUM(retornos**2) / REAL(n, dp), 1.0D-12))
      log_sigma2(1) = ln_s2_prev
      DO t = 2, n
         z_prev = retornos(t-1) / SQRT(EXP(log_sigma2(t-1)) + 1.0D-15)
         log_sigma2(t) = omega + alpha * (ABS(z_prev) - E_ABS_Z) + &
            gamma_asim * z_prev + beta * log_sigma2(t-1)
      END DO
   END SUBROUTINE egarch11

   !----------------------------------------------------------------------------
   ! GJR-GARCH(1,1) — Glosten, Jagannathan, Runkle (1993)
   ! sigma_t^2 = omega + (alpha + gamma*I_{t-1})*eps_{t-1}^2 + beta*sigma_{t-1}^2
   ! I_{t-1} = 1 se eps_{t-1} < 0 (noticias negativas têm maior impacto)
   !----------------------------------------------------------------------------
   SUBROUTINE gjr_garch11(retornos, omega, alpha, gamma_lev, beta, sigma2_out)
      REAL(dp), INTENT(IN)  :: retornos(:), omega, alpha, gamma_lev, beta
      REAL(dp), INTENT(OUT) :: sigma2_out(:)
      REAL(dp) :: indicador
      INTEGER :: n, t
      n = SIZE(retornos)
      sigma2_out(1) = SUM(retornos**2) / REAL(n, dp)
      DO t = 2, n
         indicador = MERGE(1.0_dp, 0.0_dp, retornos(t-1) < 0.0_dp)
         sigma2_out(t) = omega + (alpha + gamma_lev * indicador) * retornos(t-1)**2 + &
            beta * sigma2_out(t-1)
         sigma2_out(t) = MAX(sigma2_out(t), 1.0D-12)
      END DO
   END SUBROUTINE gjr_garch11

   !----------------------------------------------------------------------------
   ! Log-verossimilhança gaussiana do GARCH (para estimação MLE)
   ! L = -0.5 * sum_t [ln(sigma_t^2) + epsilon_t^2/sigma_t^2]
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION loglike_garch(retornos, sigma2, n) RESULT(LL)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: retornos(n), sigma2(n)
      REAL(dp) :: PI_val
      INTEGER :: t
      PI_val = 3.14159265358979323846_dp
      LL = 0.0_dp
      DO t = 1, n
         IF (sigma2(t) > 1.0D-15) &
            LL = LL - 0.5_dp * (LOG(2.0_dp * PI_val) + LOG(sigma2(t)) + &
            retornos(t)**2 / sigma2(t))
      END DO
   END FUNCTION loglike_garch

   !----------------------------------------------------------------------------
   ! Previsão multi-passo da variância condicional (h passos à frente)
   ! sigma_{t+h}^2 = omega/(1-alpha-beta) + (alpha+beta)^{h-1}*(sigma_t^2 - var_lp)
   !----------------------------------------------------------------------------
   SUBROUTINE previsao_variancia_garch(omega, alpha, beta, sigma2_t, h, previsoes)
      REAL(dp), INTENT(IN)  :: omega, alpha, beta, sigma2_t
      INTEGER,  INTENT(IN)  :: h
      REAL(dp), INTENT(OUT) :: previsoes(h)
      REAL(dp) :: var_lp, soma_ab
      INTEGER :: k
      soma_ab = alpha + beta
      IF (soma_ab < 1.0_dp) THEN
         var_lp = omega / (1.0_dp - soma_ab)
      ELSE
         var_lp = SUM(previsoes) / REAL(h, dp)  ! fallback
      END IF
      DO k = 1, h
         previsoes(k) = var_lp + soma_ab**(k-1) * (sigma2_t - var_lp)
         previsoes(k) = MAX(previsoes(k), 1.0D-12)
      END DO
   END SUBROUTINE previsao_variancia_garch

   !----------------------------------------------------------------------------
   ! Persistência do GARCH(1,1): alpha + beta (quanto > 1, mais persistente)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION persistencia_garch(alpha, beta) RESULT(pers)
      REAL(dp), INTENT(IN) :: alpha, beta
      pers = alpha + beta
   END FUNCTION persistencia_garch

   !----------------------------------------------------------------------------
   ! Half-life da volatilidade (meia-vida do choque de volatilidade)
   ! hl = ln(0.5) / ln(alpha + beta)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION meia_vida_garch(alpha, beta) RESULT(hl)
      REAL(dp), INTENT(IN) :: alpha, beta
      REAL(dp) :: soma_ab
      soma_ab = alpha + beta
      IF (soma_ab < 1.0_dp .AND. soma_ab > 0.0_dp) THEN
         hl = LOG(0.5_dp) / LOG(soma_ab)
      ELSE
         hl = 0.0_dp
      END IF
   END FUNCTION meia_vida_garch

END MODULE mod_garch
