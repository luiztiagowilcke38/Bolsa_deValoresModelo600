!==============================================================================
! Módulo 14: VaR e Modelos de Risco
! Autor: Luiz Tiago Wilcke
! Descrição: Value at Risk (VaR), Conditional VaR (CVaR/ES), VaR delta-gama,
!            stress testing, backtesting de Kupiec e Christoffersen,
!            e medidas coerentes de risco para carteiras da B3.
!==============================================================================
MODULE mod_var
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)
   REAL(dp), PARAMETER :: PI = 3.14159265358979323846_dp

CONTAINS

   !----------------------------------------------------------------------------
   ! VaR Delta-Normal (método analítico paramétrico)
   ! VaR = -W * [mu*h - z_alpha * sigma*sqrt(h)]
   ! h = horizonte temporal (dias)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION var_delta_normal(W, mu, sigma, z_alpha, h) RESULT(var_dn)
      REAL(dp), INTENT(IN) :: W, mu, sigma, z_alpha, h
      var_dn = -W * (mu * h - z_alpha * sigma * SQRT(h))
   END FUNCTION var_delta_normal

   !----------------------------------------------------------------------------
   ! VaR Delta-Gama para opções (correção por convexidade da opção)
   ! VaR_DG ≈ -delta*dS - 0.5*gamma*dS^2
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION var_delta_gamma(delta, gamma_opcao, dS) RESULT(var_dg)
      REAL(dp), INTENT(IN) :: delta, gamma_opcao, dS
      var_dg = -(delta * dS + 0.5_dp * gamma_opcao * dS**2)
   END FUNCTION var_delta_gamma

   !----------------------------------------------------------------------------
   ! Teste de Proporção de Falhas de Kupiec
   ! LR_uc = -2*ln[(1-p)^{T-n}*p^n] + 2*ln[(1-n/T)^{T-n}*(n/T)^n]
   ! H0: proporção de exceções = 1 - nível de confiança
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION teste_kupiec(n_excecoes, T_total, nivel_conf) RESULT(LR_uc)
      INTEGER,  INTENT(IN) :: n_excecoes, T_total
      REAL(dp), INTENT(IN) :: nivel_conf
      REAL(dp) :: p, p_hat, n_real, T_real
      n_real = REAL(n_excecoes, dp)
      T_real = REAL(T_total, dp)
      p      = 1.0_dp - nivel_conf    ! probabilidade esperada de exceder
      p_hat  = n_real / T_real        ! proporção observada
      IF (p_hat < 1.0D-10 .OR. ABS(1.0_dp - p_hat) < 1.0D-10) THEN
         LR_uc = 0.0_dp
         RETURN
      END IF
      LR_uc = -2.0_dp * (n_real * LOG(p / p_hat) + &
         (T_real - n_real) * LOG((1.0_dp - p) / (1.0_dp - p_hat)))
   END FUNCTION teste_kupiec

   !----------------------------------------------------------------------------
   ! Stress Testing: calcula perda sob cenário adverso
   ! Aplica choque de sigma_choque desvios padrão ao retorno
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION perda_stress(W, sigma_carteira, n_sigmas) RESULT(perda)
      REAL(dp), INTENT(IN) :: W, sigma_carteira, n_sigmas
      perda = W * sigma_carteira * n_sigmas
   END FUNCTION perda_stress

   !----------------------------------------------------------------------------
   ! Risco Condicional Implícito na opção (model-free)
   ! Estimativa de risco neutro via preços de opções (simplificado)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION variancia_risco_neutro(pesos_puts, retornos_puts, &
      pesos_calls, retornos_calls, &
      n_puts, n_calls, S, r, T) &
      RESULT(var_rn)
      INTEGER,  INTENT(IN) :: n_puts, n_calls
      REAL(dp), INTENT(IN) :: pesos_puts(n_puts), retornos_puts(n_puts)
      REAL(dp), INTENT(IN) :: pesos_calls(n_calls), retornos_calls(n_calls)
      REAL(dp), INTENT(IN) :: S, r, T
      REAL(dp) :: soma
      soma    = 2.0_dp * EXP(r * T) * (DOT_PRODUCT(pesos_puts, retornos_puts) + &
         DOT_PRODUCT(pesos_calls, retornos_calls))
      var_rn  = soma / (S * S)
   END FUNCTION variancia_risco_neutro

   !----------------------------------------------------------------------------
   ! Janela rolante de VaR Histórico
   !----------------------------------------------------------------------------
   SUBROUTINE var_historico_rolante(retornos, n, janela, nivel_conf, W, var_serie)
      INTEGER,  INTENT(IN)  :: n, janela
      REAL(dp), INTENT(IN)  :: retornos(n), nivel_conf, W
      REAL(dp), INTENT(OUT) :: var_serie(n)
      REAL(dp) :: ret_janela(janela), temp
      INTEGER :: i, j, k, idx_var
      var_serie = 0.0_dp
      DO i = janela, n
         ret_janela = retornos(i-janela+1:i)
         ! Ordenação
         DO j = 2, janela
            temp = ret_janela(j);  k = j - 1
            DO WHILE (k >= 1 .AND. ret_janela(k) > temp)
               ret_janela(k+1) = ret_janela(k);  k = k - 1
            END DO
            ret_janela(k+1) = temp
         END DO
         idx_var = MAX(1, INT((1.0_dp - nivel_conf) * REAL(janela, dp)))
         var_serie(i) = -W * ret_janela(idx_var)
      END DO
   END SUBROUTINE var_historico_rolante

   !----------------------------------------------------------------------------
   ! Risk Adjusted Return on Capital (RAROC)
   ! RAROC = Retorno_Liquido / VaR
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION raroc(retorno_liquido, var_calculado) RESULT(rc)
      REAL(dp), INTENT(IN) :: retorno_liquido, var_calculado
      IF (ABS(var_calculado) > 1.0D-15) THEN
         rc = retorno_liquido / var_calculado
      ELSE
         rc = 0.0_dp
      END IF
   END FUNCTION raroc

   !----------------------------------------------------------------------------
   ! Cálculo de Capital Econômico sob Basileia II (modelo IRM interno)
   ! CE = LGD * EAD * N[N^{-1}(PD) + rho^{0.5}*N^{-1}(0.999)] / sqrt(1-rho)
   ! Simplificado: usa aproximação normal inversa
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION capital_economico_basileia(PD, LGD, EAD, rho) RESULT(CE)
      REAL(dp), INTENT(IN) :: PD, LGD, EAD, rho
      REAL(dp) :: q_pd, q_999, phi_arg
      ! Aproximação da inversa da normal cumulativa (probit)
      q_pd  = SQRT(2.0_dp) * 1.645_dp * (PD - 0.5_dp)   ! simplificado
      q_999 = 3.09_dp      ! quantil 99.9% normal padrão
      phi_arg = (q_pd + SQRT(rho) * q_999) / SQRT(MAX(1.0_dp - rho, 1.0D-6))
      ! N(phi_arg) aproximado
      CE = LGD * EAD * (0.5_dp * (1.0_dp + ERF_APPROX(phi_arg / SQRT(2.0_dp))) - PD)
   END FUNCTION capital_economico_basileia

   !----------------------------------------------------------------------------
   ! Aproximação da função ERF (error function) para uso interno
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION ERF_APPROX(x) RESULT(erf_val)
      REAL(dp), INTENT(IN) :: x
      REAL(dp) :: t
      t = 1.0_dp / (1.0_dp + 0.47047_dp * ABS(x))
      erf_val = 1.0_dp - t * (0.3480242_dp + t * (-0.0958798_dp + t * 0.7478556_dp)) * &
         EXP(-x * x)
      IF (x < 0.0_dp) erf_val = -erf_val
   END FUNCTION ERF_APPROX

   !----------------------------------------------------------------------------
   ! Contribuição ao VaR do componente i: CVaR_incremental_i
   ! dVaR/dw_i ≈ sigma_i * rho_{i,p} / sigma_p * VaR
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION contribuicao_var(peso_i, sigma_i, rho_ip, sigma_port, &
      var_total) RESULT(cvar_i)
      REAL(dp), INTENT(IN) :: peso_i, sigma_i, rho_ip, sigma_port, var_total
      cvar_i = peso_i * sigma_i * rho_ip / (sigma_port + 1.0D-15) * var_total
   END FUNCTION contribuicao_var

END MODULE mod_var
