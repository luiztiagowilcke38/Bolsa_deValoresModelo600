!==============================================================================
! Módulo 28-29-30: Modelos de Taxas de Juro (Hull-White, Vasicek, CIR)
! Autor: Luiz Tiago Wilcke
! Descrição: Modelos de estrutura a termo de taxas de juros relevantes para
!            o mercado brasileiro (DI, CDI, Selic). Hull-White estendido
!            (calibrado à curva de juros), Vasicek analítico e Cox-
!            Ingersoll-Ross (CIR) com precificação de títulos e opções.
!==============================================================================
MODULE mod_taxas_juros
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)
   REAL(dp), PARAMETER :: PI = 3.14159265358979323846_dp

CONTAINS

   !============ MODELO DE VASICEK ============================================

   !----------------------------------------------------------------------------
   ! Preço analítico de título zero-cupom no modelo de Vasicek
   ! P(t,T) = A(t,T) * exp(-B(t,T)*r_t)
   ! B(t,T) = (1 - e^{-kappa*(T-t)}) / kappa
   ! A(t,T) = exp[(B-tau)*(kappa^2*theta-sigma^2/2)/kappa^2 - sigma^2*B^2/(4*kappa)]
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION preco_zero_cupom_vasicek(r_t, kappa_v, theta_v, sigma_v, tau) &
      RESULT(P)
      REAL(dp), INTENT(IN) :: r_t, kappa_v, theta_v, sigma_v, tau
      REAL(dp) :: B_tau, expoente_A, A
      B_tau = (1.0_dp - EXP(-kappa_v * tau)) / (kappa_v + 1.0D-15)
      expoente_A = (B_tau - tau) * (kappa_v**2 * theta_v - sigma_v**2/2.0_dp) / &
         (kappa_v**2 + 1.0D-15) - sigma_v**2 * B_tau**2 / (4.0_dp * kappa_v + 1.0D-15)
      A = EXP(expoente_A)
      P = A * EXP(-B_tau * r_t)
   END FUNCTION preco_zero_cupom_vasicek

   !----------------------------------------------------------------------------
   ! Taxa forward instantânea no modelo Vasicek
   ! f(0,T) = -d ln P(0,T) / dT
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION taxa_forward_vasicek(r0, kappa_v, theta_v, sigma_v, T) &
      RESULT(f)
      REAL(dp), INTENT(IN) :: r0, kappa_v, theta_v, sigma_v, T
      REAL(dp) :: exp_kT
      exp_kT = EXP(-kappa_v * T)
      f = theta_v + (r0 - theta_v) * exp_kT + &
         sigma_v**2 / (2.0_dp * kappa_v**2 + 1.0D-15) * (1.0_dp - exp_kT)**2
   END FUNCTION taxa_forward_vasicek

   !----------------------------------------------------------------------------
   ! Simulação de trajetória Vasicek (solução exata)
   !----------------------------------------------------------------------------
   SUBROUTINE simular_vasicek(r0, kappa_v, theta_v, sigma_v, dt, n_passos, traj)
      REAL(dp), INTENT(IN)  :: r0, kappa_v, theta_v, sigma_v, dt
      INTEGER,  INTENT(IN)  :: n_passos
      REAL(dp), INTENT(OUT) :: traj(n_passos+1)
      REAL(dp) :: r, e_kdt, media_cond, sigma_cond, z, u1, u2
      INTEGER :: i
      r = r0;  traj(1) = r
      e_kdt      = EXP(-kappa_v * dt)
      sigma_cond = sigma_v * SQRT((1.0_dp - e_kdt**2) / (2.0_dp*kappa_v + 1.0D-15))
      DO i = 1, n_passos
         media_cond = theta_v + (r - theta_v) * e_kdt
         CALL RANDOM_NUMBER(u1);  CALL RANDOM_NUMBER(u2)
         u1 = MAX(u1, 1.0D-15)
         z  = SQRT(-2.0_dp * LOG(u1)) * COS(2.0_dp * PI * u2)
         r  = media_cond + sigma_cond * z
         traj(i+1) = r
      END DO
   END SUBROUTINE simular_vasicek

   !============ MODELO CIR ====================================================

   !----------------------------------------------------------------------------
   ! Preço de título zero-cupom CIR
   ! P(t,T) = A_cir(tau) * exp(-B_cir(tau)*r_t)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION preco_zero_cupom_cir(r_t, kappa_c, theta_c, sigma_c, tau) &
      RESULT(P)
      REAL(dp), INTENT(IN) :: r_t, kappa_c, theta_c, sigma_c, tau
      REAL(dp) :: gamma_c, B_cir, num_A, den_A, A_cir, expoente
      gamma_c = SQRT(kappa_c**2 + 2.0_dp * sigma_c**2)
      num_A   = 2.0_dp * gamma_c * EXP(0.5_dp*(kappa_c + gamma_c)*tau)
      den_A   = (gamma_c + kappa_c)*(EXP(gamma_c*tau) - 1.0_dp) + 2.0_dp*gamma_c
      IF (ABS(den_A) < 1.0D-15) den_A = 1.0D-15
      expoente = 2.0_dp * kappa_c * theta_c / sigma_c**2
      A_cir = (num_A / den_A)**expoente
      B_cir = 2.0_dp*(EXP(gamma_c*tau) - 1.0_dp) / (den_A + 1.0D-15)
      P = A_cir * EXP(-B_cir * r_t)
   END FUNCTION preco_zero_cupom_cir

   !============ MODELO HULL-WHITE =============================================

   !----------------------------------------------------------------------------
   ! Hull-White estendido (theta(t) calibrado à curva de juros inicial)
   ! dr = [theta(t) - kappa*r]*dt + sigma*dW
   ! theta(t) = df(0,t)/dt + kappa*f(0,t) + sigma^2/(2*kappa)*(1-e^{-2*kappa*t})
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION theta_hull_white(t, kappa_hw, sigma_hw, df_dt, f0t) &
      RESULT(theta_t)
      REAL(dp), INTENT(IN) :: t, kappa_hw, sigma_hw, df_dt, f0t
      theta_t = df_dt + kappa_hw * f0t + &
         sigma_hw**2 / (2.0_dp * kappa_hw + 1.0D-15) * &
         (1.0_dp - EXP(-2.0_dp * kappa_hw * t))
   END FUNCTION theta_hull_white

   !----------------------------------------------------------------------------
   ! Preço de TITULO zero-cupom em Hull-White
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION preco_zcb_hw(r_t, tau, kappa_hw, sigma_hw, &
      f0t_base, f0T_curr) RESULT(P)
      REAL(dp), INTENT(IN) :: r_t, tau, kappa_hw, sigma_hw, f0t_base, f0T_curr
      REAL(dp) :: B_hw, A_hw
      B_hw = (1.0_dp - EXP(-kappa_hw * tau)) / (kappa_hw + 1.0D-15)
      A_hw = EXP(-f0T_curr * tau + 0.5_dp * sigma_hw**2 / (kappa_hw**2 + 1.0D-15) * &
         (EXP(-2.0_dp*kappa_hw*(0.0_dp)) - 1.0_dp + 2.0_dp*B_hw - B_hw**2*kappa_hw))
      P = A_hw * EXP(-B_hw * r_t)
   END FUNCTION preco_zcb_hw

   !----------------------------------------------------------------------------
   ! Curva de juros zero IPCA+ (simulada para calibração com valores típicos da B3)
   ! Retorna taxa spot para maturidade tau (em anos) — dados aproximados 2025
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION curva_juros_b3(tau) RESULT(r_spot)
      REAL(dp), INTENT(IN) :: tau
      ! Curva simplificada baseada em taxas típicas do DI futuro Brasil
      r_spot = 0.105_dp + 0.008_dp * (1.0_dp - EXP(-tau / 5.0_dp)) - &
         0.003_dp * (1.0_dp - EXP(-tau / 20.0_dp))
   END FUNCTION curva_juros_b3

   !----------------------------------------------------------------------------
   ! Duração modificada de portfólio de títulos
   ! D_mod = -(1/P) * dP/dr
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION duracao_modificada(P_menos, P_mais, delta_r) RESULT(D_mod)
      REAL(dp), INTENT(IN) :: P_menos, P_mais, delta_r
      REAL(dp) :: P_meio
      P_meio  = (P_menos + P_mais) / 2.0_dp
      IF (ABS(P_meio * delta_r) > 1.0D-15) THEN
         D_mod = -(P_mais - P_menos) / (2.0_dp * P_meio * delta_r)
      ELSE
         D_mod = 0.0_dp
      END IF
   END FUNCTION duracao_modificada

   !----------------------------------------------------------------------------
   ! Convexidade do título em relação à taxa de juros
   ! C = (1/P) * d²P/dr²
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION convexidade(P_menos, P_centro, P_mais, delta_r) RESULT(C)
      REAL(dp), INTENT(IN) :: P_menos, P_centro, P_mais, delta_r
      IF (ABS(P_centro * delta_r**2) > 1.0D-15) THEN
         C = (P_mais - 2.0_dp*P_centro + P_menos) / (P_centro * delta_r**2)
      ELSE
         C = 0.0_dp
      END IF
   END FUNCTION convexidade

END MODULE mod_taxas_juros
