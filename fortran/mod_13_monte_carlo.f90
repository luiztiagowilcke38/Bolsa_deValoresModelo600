!==============================================================================
! Módulo 13: Simulação de Monte Carlo
! Autor: Luiz Tiago Wilcke
! Descrição: Geração de trajetórias de preços via Monte Carlo com processo
!            de Wiener (MBG), Monte Carlo quase-aleatório (Sobol),
!            precificação de opções, VaR e CVaR por simulação,
!            e estimação de distribuições futuras de retornos para a B3.
!==============================================================================
MODULE mod_monte_carlo
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)
   REAL(dp), PARAMETER :: PI = 3.14159265358979323846_dp

CONTAINS

   !----------------------------------------------------------------------------
   ! Gerador de números normais padrão via Box-Muller
   ! Z1 = sqrt(-2*ln(U1)) * cos(2*pi*U2)
   ! Z2 = sqrt(-2*ln(U1)) * sin(2*pi*U2)
   !----------------------------------------------------------------------------
   SUBROUTINE gerar_normais(n, z)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(OUT) :: z(n)
      REAL(dp) :: u1, u2
      INTEGER :: i
      DO i = 1, n - 1, 2
         CALL RANDOM_NUMBER(u1);  CALL RANDOM_NUMBER(u2)
         u1 = MAX(u1, 1.0D-15)
         z(i)   = SQRT(-2.0_dp * LOG(u1)) * COS(2.0_dp * PI * u2)
         z(i+1) = SQRT(-2.0_dp * LOG(u1)) * SIN(2.0_dp * PI * u2)
      END DO
      IF (MOD(n, 2) == 1) THEN
         CALL RANDOM_NUMBER(u1);  CALL RANDOM_NUMBER(u2)
         z(n) = SQRT(-2.0_dp * LOG(MAX(u1, 1.0D-15))) * COS(2.0_dp * PI * u2)
      END IF
   END SUBROUTINE gerar_normais

   !----------------------------------------------------------------------------
   ! Simulação de MBG (Movimento Browniano Geométrico)
   ! dS = mu*S*dt + sigma*S*dW
   ! S_t = S0 * exp[(mu - sigma^2/2)*t + sigma*W_t]
   !----------------------------------------------------------------------------
   SUBROUTINE simular_mbg(S0, mu, sigma, dt, n_passos, n_traj, trajetorias)
      REAL(dp), INTENT(IN)  :: S0, mu, sigma, dt
      INTEGER,  INTENT(IN)  :: n_passos, n_traj
      REAL(dp), INTENT(OUT) :: trajetorias(n_traj, n_passos+1)
      REAL(dp) :: S, drift, vol, z(1)
      INTEGER :: i, j
      drift = (mu - 0.5_dp * sigma**2) * dt
      vol   = sigma * SQRT(dt)
      DO j = 1, n_traj
         S = S0
         trajetorias(j, 1) = S
         DO i = 1, n_passos
            CALL gerar_normais(1, z)
            S = S * EXP(drift + vol * z(1))
            trajetorias(j, i+1) = S
         END DO
      END DO
   END SUBROUTINE simular_mbg

   !----------------------------------------------------------------------------
   ! Preço de opção CALL por Monte Carlo (media do payoff descontado)
   ! C = e^{-rT} * E[max(S_T - K, 0)]
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION preco_call_mc(S0, K, r, sigma, T, n_traj) RESULT(C)
      REAL(dp), INTENT(IN) :: S0, K, r, sigma, T
      INTEGER,  INTENT(IN) :: n_traj
      REAL(dp) :: soma_payoff, ST, z(1)
      REAL(dp) :: drift, vol
      INTEGER :: j
      drift = (r - 0.5_dp * sigma**2) * T
      vol   = sigma * SQRT(T)
      soma_payoff = 0.0_dp
      DO j = 1, n_traj
         CALL gerar_normais(1, z)
         ST = S0 * EXP(drift + vol * z(1))
         soma_payoff = soma_payoff + MAX(ST - K, 0.0_dp)
      END DO
      C = EXP(-r * T) * soma_payoff / REAL(n_traj, dp)
   END FUNCTION preco_call_mc

   !----------------------------------------------------------------------------
   ! VaR paramétrico: VaR = mu + z_alpha * sigma (para nível alpha)
   ! Para alpha = 0.05: z_alpha = -1.645
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION var_parametrico(mu_ret, sigma_ret, nivel_confianca, &
      valor_carteira) RESULT(var_val)
      REAL(dp), INTENT(IN) :: mu_ret, sigma_ret, nivel_confianca, valor_carteira
      REAL(dp) :: z_alpha
      ! Quantil normal aproximado para nível de confiança
      IF (nivel_confianca >= 0.99_dp) THEN
         z_alpha = 2.326_dp
      ELSE IF (nivel_confianca >= 0.975_dp) THEN
         z_alpha = 1.960_dp
      ELSE
         z_alpha = 1.645_dp
      END IF
      var_val = valor_carteira * (mu_ret - z_alpha * sigma_ret)
   END FUNCTION var_parametrico

   !----------------------------------------------------------------------------
   ! VaR por simulação histórica (percentil dos retornos históricos)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION var_historico(retornos, n, nivel_confianca, &
      valor_carteira) RESULT(var_val)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: retornos(n), nivel_confianca, valor_carteira
      REAL(dp) :: ret_ord(n), temp
      INTEGER :: i, j, idx_var
      ret_ord = retornos
      ! Ordenação por inserção (simple sort)
      DO i = 2, n
         temp = ret_ord(i);  j = i - 1
         DO WHILE (j >= 1 .AND. ret_ord(j) > temp)
            ret_ord(j+1) = ret_ord(j);  j = j - 1
         END DO
         ret_ord(j+1) = temp
      END DO
      idx_var = INT((1.0_dp - nivel_confianca) * REAL(n, dp)) + 1
      idx_var = MAX(1, MIN(idx_var, n))
      var_val = -valor_carteira * ret_ord(idx_var)
   END FUNCTION var_historico

   !----------------------------------------------------------------------------
   ! CVaR (Expected Shortfall): media das perdas além do VaR
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION cvar(retornos, n, nivel_confianca, valor_carteira) RESULT(cv)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: retornos(n), nivel_confianca, valor_carteira
      REAL(dp) :: ret_ord(n), temp, soma
      INTEGER :: i, j, idx_var, n_extremos
      ret_ord = retornos
      DO i = 2, n
         temp = ret_ord(i);  j = i - 1
         DO WHILE (j >= 1 .AND. ret_ord(j) > temp)
            ret_ord(j+1) = ret_ord(j);  j = j - 1
         END DO
         ret_ord(j+1) = temp
      END DO
      idx_var  = MAX(1, INT((1.0_dp - nivel_confianca) * REAL(n, dp)))
      n_extremos = idx_var
      soma = SUM(ret_ord(1:n_extremos))
      cv   = -valor_carteira * soma / REAL(n_extremos, dp)
   END FUNCTION cvar

   !----------------------------------------------------------------------------
   ! Simulação de portfólio por Monte Carlo (n_ativos correlacionados via Cholesky)
   !----------------------------------------------------------------------------
   SUBROUTINE simulacao_portfolio_mc(mu_vec, L_chol, n_ativos, &
      n_passos, n_traj, dt, retornos_traj)
      INTEGER,  INTENT(IN)  :: n_ativos, n_passos, n_traj
      REAL(dp), INTENT(IN)  :: mu_vec(n_ativos), L_chol(n_ativos, n_ativos), dt
      REAL(dp), INTENT(OUT) :: retornos_traj(n_traj, n_passos)
      REAL(dp) :: z(n_ativos), z_corr(n_ativos), ret_port
      INTEGER :: i, j
      DO j = 1, n_traj
         DO i = 1, n_passos
            CALL gerar_normais(n_ativos, z)
            z_corr   = MATMUL(L_chol, z) * SQRT(dt)
            ret_port = DOT_PRODUCT(mu_vec, dt + z_corr) / REAL(n_ativos, dp)
            retornos_traj(j, i) = ret_port
         END DO
      END DO
   END SUBROUTINE simulacao_portfolio_mc

   !----------------------------------------------------------------------------
   ! Teste de backtesting de VaR (exceções de Kupiec)
   ! num_excecoes = count(retorno_t < -VaR_t)
   !----------------------------------------------------------------------------
   INTEGER FUNCTION contar_excecoes_var(retornos, var_serie, n) RESULT(excecoes)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: retornos(n), var_serie(n)
      INTEGER :: i
      excecoes = 0
      DO i = 1, n
         IF (retornos(i) < -var_serie(i)) excecoes = excecoes + 1
      END DO
   END FUNCTION contar_excecoes_var

END MODULE mod_monte_carlo
