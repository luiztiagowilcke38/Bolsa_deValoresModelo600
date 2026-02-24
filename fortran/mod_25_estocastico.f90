!==============================================================================
! Módulo 25: Processos Estocásticos Financeiros
! Autor: Luiz Tiago Wilcke
! Descrição: Processos de Ornstein-Uhlenbeck (OU), Processo de Vasicek,
!            Processo de Cox-Ingersoll-Ross (CIR), processos de salto
!            (Merton jump-diffusion), processos de variância estocástica
!            e processos de Levy para modelagem de taxas e crédito na B3.
!==============================================================================
MODULE mod_estocastico
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)
   REAL(dp), PARAMETER :: PI = 3.14159265358979323846_dp

CONTAINS

   !----------------------------------------------------------------------------
   ! Processo de Ornstein-Uhlenbeck (OU)
   ! dX = kappa*(theta - X)*dt + sigma*dW
   ! Solução exata: X_t = theta + (X_0 - theta)*e^{-kappa*t} + volatilidade
   !----------------------------------------------------------------------------
   SUBROUTINE processo_ou(X0, kappa, theta, sigma, dt, n_passos, n_traj, trajetorias)
      REAL(dp), INTENT(IN)  :: X0, kappa, theta, sigma, dt
      INTEGER,  INTENT(IN)  :: n_passos, n_traj
      REAL(dp), INTENT(OUT) :: trajetorias(n_traj, n_passos+1)
      REAL(dp) :: X, z, sigma_eff, theta_eff
      INTEGER :: i, j
      DO j = 1, n_traj
         X = X0
         trajetorias(j,1) = X
         DO i = 1, n_passos
            sigma_eff = sigma * SQRT((1.0_dp - EXP(-2.0_dp*kappa*dt)) / (2.0_dp*kappa+1.0D-15))
            theta_eff = theta + (X - theta) * EXP(-kappa * dt)
            CALL gerar_normal_escalar(z)
            X = theta_eff + sigma_eff * z
            trajetorias(j, i+1) = X
         END DO
      END DO
   END SUBROUTINE processo_ou

   !----------------------------------------------------------------------------
   ! Processo de difusão com saltos de Merton (Jump-Diffusion)
   ! dS = (mu - lambda*mu_J)*S*dt + sigma*S*dW + (e^J - 1)*S*dN
   ! Onde N é processo de Poisson com intensidade lambda
   !          J ~ N(mu_J, sigma_J^2) é o tamanho do salto log-normal
   !----------------------------------------------------------------------------
   SUBROUTINE merton_jump_diffusion(S0, mu, sigma, lambda_salto, &
      mu_j, sigma_j, dt, n_passos, n_traj, traj)
      REAL(dp), INTENT(IN)  :: S0, mu, sigma, lambda_salto, mu_j, sigma_j, dt
      INTEGER,  INTENT(IN)  :: n_passos, n_traj
      REAL(dp), INTENT(OUT) :: traj(n_traj, n_passos+1)
      REAL(dp) :: S, drift, vol, z, J_tamanho, u_poisson, prob_salto
      INTEGER :: i, j
      prob_salto = 1.0_dp - EXP(-lambda_salto * dt)
      DO j = 1, n_traj
         S = S0
         traj(j,1) = S
         DO i = 1, n_passos
            CALL gerar_normal_escalar(z)
            drift = (mu - lambda_salto * (EXP(mu_j + 0.5_dp*sigma_j**2) - 1.0_dp) &
               - 0.5_dp * sigma**2) * dt
            vol   = sigma * SQRT(dt) * z
            S     = S * EXP(drift + vol)
            ! Componente de salto
            CALL RANDOM_NUMBER(u_poisson)
            IF (u_poisson < prob_salto) THEN
               CALL gerar_normal_escalar(J_tamanho)
               J_tamanho = mu_j + sigma_j * J_tamanho
               S = S * EXP(J_tamanho)
            END IF
            traj(j, i+1) = MAX(S, 1.0D-10)
         END DO
      END DO
   END SUBROUTINE merton_jump_diffusion

   !----------------------------------------------------------------------------
   ! Processo Cox-Ingersoll-Ross (CIR) para taxas de juros
   ! dr = kappa*(theta - r)*dt + sigma*sqrt(r)*dW
   ! Discretização de Milstein para estabilidade
   !----------------------------------------------------------------------------
   SUBROUTINE processo_cir(r0, kappa, theta, sigma, dt, n_passos, n_traj, traj_r)
      REAL(dp), INTENT(IN)  :: r0, kappa, theta, sigma, dt
      INTEGER,  INTENT(IN)  :: n_passos, n_traj
      REAL(dp), INTENT(OUT) :: traj_r(n_traj, n_passos+1)
      REAL(dp) :: r, z, drift, vol, milstein_corr
      INTEGER :: i, j
      DO j = 1, n_traj
         r = MAX(r0, 0.0_dp)
         traj_r(j,1) = r
         DO i = 1, n_passos
            CALL gerar_normal_escalar(z)
            drift        = kappa * (theta - r) * dt
            vol          = sigma * SQRT(MAX(r, 0.0_dp) * dt) * z
            milstein_corr = 0.25_dp * sigma**2 * dt * (z**2 - 1.0_dp)
            r = r + drift + vol + milstein_corr
            r = MAX(r, 0.0_dp)   ! garante não-negatividade
            traj_r(j, i+1) = r
         END DO
      END DO
   END SUBROUTINE processo_cir

   !----------------------------------------------------------------------------
   ! Gerador de varíavel aleatória N(0,1) pelo método de Box-Muller
   !----------------------------------------------------------------------------
   SUBROUTINE gerar_normal_escalar(z)
      REAL(dp), INTENT(OUT) :: z
      REAL(dp) :: u1, u2
      CALL RANDOM_NUMBER(u1);  CALL RANDOM_NUMBER(u2)
      u1 = MAX(u1, 1.0D-15)
      z  = SQRT(-2.0_dp * LOG(u1)) * COS(2.0_dp * PI * u2)
   END SUBROUTINE gerar_normal_escalar

   !----------------------------------------------------------------------------
   ! Simulação de spread de crédito (Vasicek duplo — modelo de dois fatores)
   ! dr = kappa1*(theta1-r)*dt + sigma1*dW1
   ! ds = kappa2*(theta2-s)*dt + sigma2*dW2 + rho*dW1
   !----------------------------------------------------------------------------
   SUBROUTINE spread_credito_dois_fatores(r0, s0, kappa1, theta1, sigma1, &
      kappa2, theta2, sigma2, rho, &
      dt, n_passos, traj_r, traj_s)
      REAL(dp), INTENT(IN)  :: r0, s0, kappa1, theta1, sigma1
      REAL(dp), INTENT(IN)  :: kappa2, theta2, sigma2, rho, dt
      INTEGER,  INTENT(IN)  :: n_passos
      REAL(dp), INTENT(OUT) :: traj_r(n_passos+1), traj_s(n_passos+1)
      REAL(dp) :: r, s, z1, z2
      INTEGER :: i
      r = r0;  s = s0
      traj_r(1) = r;  traj_s(1) = s
      DO i = 1, n_passos
         CALL gerar_normal_escalar(z1);  CALL gerar_normal_escalar(z2)
         z2 = rho * z1 + SQRT(MAX(1.0_dp - rho**2, 0.0_dp)) * z2
         r  = r + kappa1*(theta1 - r)*dt + sigma1*SQRT(dt)*z1
         s  = s + kappa2*(theta2 - s)*dt + sigma2*SQRT(dt)*z2
         traj_r(i+1) = r;  traj_s(i+1) = s
      END DO
   END SUBROUTINE spread_credito_dois_fatores

   !----------------------------------------------------------------------------
   ! Variância condicional pelo modelo SABR (para vol estocástica implícita)
   ! dF = sigma*F^beta*dW1
   ! dsigma = alpha*sigma*dW2   com corr(dW1,dW2) = rho_sabr
   !----------------------------------------------------------------------------
   SUBROUTINE simular_sabr(F0, sigma0, alpha_sabr, beta_sabr, rho_sabr, &
      dt, n_passos, traj_F, traj_sigma)
      REAL(dp), INTENT(IN)  :: F0, sigma0, alpha_sabr, beta_sabr, rho_sabr, dt
      INTEGER,  INTENT(IN)  :: n_passos
      REAL(dp), INTENT(OUT) :: traj_F(n_passos+1), traj_sigma(n_passos+1)
      REAL(dp) :: F, sig, z1, z2
      INTEGER :: i
      F = F0;  sig = sigma0
      traj_F(1) = F;  traj_sigma(1) = sig
      DO i = 1, n_passos
         CALL gerar_normal_escalar(z1);  CALL gerar_normal_escalar(z2)
         z2 = rho_sabr * z1 + SQRT(MAX(1.0_dp - rho_sabr**2, 0.0_dp)) * z2
         F   = F   + sig * F**beta_sabr * SQRT(dt) * z1
         sig = sig + alpha_sabr * sig * SQRT(dt) * z2
         F   = MAX(F, 1.0D-8)
         sig = MAX(sig, 1.0D-8)
         traj_F(i+1) = F;  traj_sigma(i+1) = sig
      END DO
   END SUBROUTINE simular_sabr

END MODULE mod_estocastico
