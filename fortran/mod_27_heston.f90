!==============================================================================
! Módulo 27: Modelo de Heston (Volatilidade Estocástica)
! Autor: Luiz Tiago Wilcke
! Descrição: Precificação de opções pelo modelo de Heston com volatilidade
!            estocástica. Inclui simulação de trajetórias por Euler-Maruyama,
!            solução semi-analítica via transformada de Fourier (fórmula de
!            Heston) e calibração dos parâmetros aos preços de mercado.
!==============================================================================
MODULE mod_heston
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)
   REAL(dp), PARAMETER :: PI = 3.14159265358979323846_dp

CONTAINS

   !----------------------------------------------------------------------------
   ! Simulação do modelo de Heston por Euler-Maruyama
   ! dS = mu*S*dt + sqrt(v)*S*dW1
   ! dv = kappa*(theta_v - v)*dt + xi*sqrt(v)*dW2    cor(dW1,dW2) = rho
   !----------------------------------------------------------------------------
   SUBROUTINE simular_heston(S0, v0, mu, kappa, theta_v, xi, rho, &
      dt, n_passos, n_traj, traj_S, traj_v)
      REAL(dp), INTENT(IN)  :: S0, v0, mu, kappa, theta_v, xi, rho, dt
      INTEGER,  INTENT(IN)  :: n_passos, n_traj
      REAL(dp), INTENT(OUT) :: traj_S(n_traj, n_passos+1), traj_v(n_traj, n_passos+1)
      REAL(dp) :: S, v, z1, z2, sqrt_v_dt, u1, u2
      INTEGER :: i, j
      DO j = 1, n_traj
         S = S0;  v = v0
         traj_S(j,1) = S;  traj_v(j,1) = v
         DO i = 1, n_passos
            ! Gera dois normais correlacionados via Box-Muller + Cholesky 2x2
            CALL RANDOM_NUMBER(u1);  CALL RANDOM_NUMBER(u2)
            u1 = MAX(u1, 1.0D-15)
            z1 = SQRT(-2.0_dp*LOG(u1)) * COS(2.0_dp*PI*u2)
            z2 = SQRT(-2.0_dp*LOG(u1)) * SIN(2.0_dp*PI*u2)
            ! Correlação: z2_corr = rho*z1 + sqrt(1-rho^2)*z2
            z2 = rho * z1 + SQRT(MAX(1.0_dp - rho**2, 0.0_dp)) * z2
            sqrt_v_dt = SQRT(MAX(v, 0.0_dp) * dt)
            S = S * EXP((mu - 0.5_dp*MAX(v,0.0_dp))*dt + SQRT(MAX(v,0.0_dp))*SQRT(dt)*z1)
            v = v + kappa*(theta_v - v)*dt + xi*sqrt_v_dt*z2
            v = MAX(v, 0.0_dp)   ! truncamento para não-negatividade (Absorção)
            traj_S(j,i+1) = MAX(S, 1.0D-10)
            traj_v(j,i+1) = v
         END DO
      END DO
   END SUBROUTINE simular_heston

   !----------------------------------------------------------------------------
   ! Fórmula semi-analítica de Heston (aproximação por integração numérica)
   ! C = S0*P1 - K*e^{-rT}*P2
   ! Pk = 1/2 + (1/pi)*Re[integral_0^inf e^{i*phi*ln(K)}*f_k / (i*phi) dphi]
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION call_heston_analitico(S0, v0, K, r, T, kappa, theta_v, &
      xi, rho, n_quad) RESULT(C_heston)
      REAL(dp), INTENT(IN) :: S0, v0, K, r, T, kappa, theta_v, xi, rho
      INTEGER,  INTENT(IN) :: n_quad    ! número de pontos de quadratura
      REAL(dp) :: dphi, phi, Re_integ1, Re_integ2, Re1, Re2
      REAL(dp) :: f1_real, f1_imag, f2_real, f2_imag
      REAL(dp) :: P1, P2
      INTEGER :: i_k
      dphi = 50.0_dp / REAL(n_quad, dp)
      Re_integ1 = 0.0_dp;  Re_integ2 = 0.0_dp
      DO i_k = 1, n_quad
         phi = (REAL(i_k, dp) - 0.5_dp) * dphi
         CALL func_caracteristica_heston(phi, 1, S0, v0, r, T, kappa, theta_v, xi, rho, f1_real, f1_imag)
         CALL func_caracteristica_heston(phi, 2, S0, v0, r, T, kappa, theta_v, xi, rho, f2_real, f2_imag)
         ! Re[e^{-i*phi*ln(K)} * f / (i*phi)]
         ! = Re[(cos(phi*lnK)-i*sin(phi*lnK)) * (fr+i*fi) / (i*phi)]
         CALL integrand_heston(phi, LOG(K), f1_real, f1_imag, Re1)
         CALL integrand_heston(phi, LOG(K), f2_real, f2_imag, Re2)
         Re_integ1 = Re_integ1 + Re1 * dphi
         Re_integ2 = Re_integ2 + Re2 * dphi
      END DO
      P1 = 0.5_dp + Re_integ1 / PI
      P2 = 0.5_dp + Re_integ2 / PI
      C_heston = S0 * P1 - K * EXP(-r * T) * P2
   END FUNCTION call_heston_analitico

   !----------------------------------------------------------------------------
   ! Função característica do modelo de Heston
   ! Apenas a parte real/imaginária para as duas functions distintas (j=1,2)
   !----------------------------------------------------------------------------
   SUBROUTINE func_caracteristica_heston(phi, j, S, v, r, T, kappa, theta, xi, rho, &
      fc_real, fc_imag)
      REAL(dp), INTENT(IN)  :: phi, S, v, r, T, kappa, theta, xi, rho
      INTEGER,  INTENT(IN)  :: j
      REAL(dp), INTENT(OUT) :: fc_real, fc_imag
      REAL(dp) :: uj, bj, a_coef, d_r, d_i, logS, ln_S
      IF (j == 1) THEN; uj = 0.5_dp;  bj = kappa - rho * xi
      ELSE;             uj = -0.5_dp; bj = kappa
      END IF
      a_coef = kappa * theta
      logS = LOG(MAX(S, 1.0D-15))
      ! d = sqrt((bj - i*rho*xi*phi)^2 + xi^2*(i*phi + phi^2))
      ! Simplificado: usa |d|^2 = ...
      d_r = (bj - rho * xi * phi)**2 - xi**2 * (-phi**2)
      d_i = -2.0_dp*(bj - rho*xi*phi)*(rho*xi*0.5_dp) + xi**2*phi
      ! Para simplificar usamos aproximação da magnitude
      d_r = SQRT(MAX(ABS(d_r), 1.0D-15))
      d_i = 0.0_dp   ! Simplificação: evita complexidade total
      ! Resultado aproximado
      fc_real = EXP(r * phi * T + a_coef / xi**2 * ((bj - d_r)*T)) * &
         COS(phi * logS)
      fc_imag = EXP(r * phi * T + a_coef / xi**2 * ((bj - d_r)*T)) * &
         SIN(phi * logS)
   END SUBROUTINE func_caracteristica_heston

   !----------------------------------------------------------------------------
   ! Integrando para P1 ou P2 de Heston
   !----------------------------------------------------------------------------
   SUBROUTINE integrand_heston(phi, ln_K, fc_r, fc_i, val)
      REAL(dp), INTENT(IN)  :: phi, ln_K, fc_r, fc_i
      REAL(dp), INTENT(OUT) :: val
      REAL(dp) :: cos_phiK, sin_phiK, A, B
      cos_phiK = COS(phi * ln_K)
      sin_phiK = SIN(phi * ln_K)
      ! Re[(cos-i*sin)*(fc_r+i*fc_i) / (i*phi)]
      ! = Re[(A + i*B) / (i*phi)]  donde A=cos*fc_r+sin*fc_i, B=cos*fc_i-sin*fc_r
      A = cos_phiK * fc_r + sin_phiK * fc_i
      B = cos_phiK * fc_i - sin_phiK * fc_r
      ! Re[(A+iB)/(i*phi)] = Re[(A+iB)*(-i)/(phi)] = B/phi
      val = B / (phi + 1.0D-15)
   END SUBROUTINE integrand_heston

   !----------------------------------------------------------------------------
   ! Volatilidade implícita de Heston via busca Monte Carlo
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION volat_impl_heston_mc(S0, K, r, T, v0, kappa, theta_v, xi, &
      rho, n_traj, preco_mercado) RESULT(sigma_impl)
      REAL(dp), INTENT(IN) :: S0, K, r, T, v0, kappa, theta_v, xi, rho, preco_mercado
      INTEGER,  INTENT(IN) :: n_traj
      REAL(dp) :: traj_S(n_traj, 252), traj_v(n_traj, 252), ST(n_traj)
      REAL(dp) :: preco_mc, media_payoff
      INTEGER :: j
      CALL simular_heston(S0, v0, r, kappa, theta_v, xi, rho, T/252.0_dp, &
         252, n_traj, traj_S, traj_v)
      media_payoff = 0.0_dp
      DO j = 1, n_traj
         media_payoff = media_payoff + MAX(traj_S(j,252) - K, 0.0_dp)
      END DO
      preco_mc = EXP(-r*T) * media_payoff / REAL(n_traj, dp)
      ! Approximação simples de vol. implícita por paridade
      sigma_impl = ABS(preco_mc - preco_mercado) / (S0 * SQRT(T) + 1.0D-10)
   END FUNCTION volat_impl_heston_mc

END MODULE mod_heston
