!==============================================================================
! Módulo 03: Modelo de Black-Scholes
! Autor: Luiz Tiago Wilcke
! Descrição: Implementação completa do modelo Black-Scholes-Merton para
!            precificação de opções europeias e gregas.
!            Fórmula analítica: C = S*N(d1) - K*e^{-rT}*N(d2)
!            d1 = [ln(S/K) + (r + sigma^2/2)*T] / (sigma*sqrt(T))
!            d2 = d1 - sigma*sqrt(T)
!==============================================================================
MODULE mod_black_scholes
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)
   REAL(dp), PARAMETER :: PI = 3.14159265358979323846_dp

CONTAINS

   !----------------------------------------------------------------------------
   ! Função de distribuição normal cumulativa N(x) via aproximação polinomial
   ! de Abramowitz & Stegun (erro < 7.5e-8)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION normal_cumulativa(x) RESULT(N)
      REAL(dp), INTENT(IN) :: x
      REAL(dp) :: t, z, polinomio
      REAL(dp), PARAMETER :: b1 = 0.31938153_dp,  b2 = -0.356563782_dp
      REAL(dp), PARAMETER :: b3 = 1.781477937_dp, b4 = -1.821255978_dp
      REAL(dp), PARAMETER :: b5 = 1.330274429_dp, p  = 0.2316419_dp
      z = ABS(x) / SQRT(2.0_dp)
      t = 1.0_dp / (1.0_dp + p * ABS(x))
      polinomio = t * (b1 + t * (b2 + t * (b3 + t * (b4 + t * b5))))
      N = 1.0_dp - (EXP(-x*x/2.0_dp) / SQRT(2.0_dp * PI)) * polinomio
      IF (x < 0.0_dp) N = 1.0_dp - N
   END FUNCTION normal_cumulativa

   !----------------------------------------------------------------------------
   ! Densidade da distribuição normal padrão n(x)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION densidade_normal(x) RESULT(phi)
      REAL(dp), INTENT(IN) :: x
      phi = EXP(-0.5_dp * x * x) / SQRT(2.0_dp * PI)
   END FUNCTION densidade_normal

   !----------------------------------------------------------------------------
   ! Parâmetros d1 e d2 do modelo de Black-Scholes
   !----------------------------------------------------------------------------
   SUBROUTINE calcular_d1_d2(S, K, r, sigma, T, d1, d2)
      REAL(dp), INTENT(IN)  :: S, K, r, sigma, T
      REAL(dp), INTENT(OUT) :: d1, d2
      REAL(dp) :: sigma_sqrt_T
      sigma_sqrt_T = sigma * SQRT(T)
      d1 = (LOG(S / K) + (r + 0.5_dp * sigma**2) * T) / (sigma_sqrt_T + 1.0D-15)
      d2 = d1 - sigma_sqrt_T
   END SUBROUTINE calcular_d1_d2

   !----------------------------------------------------------------------------
   ! Preço de opção de compra (CALL europeia) - Black-Scholes
   ! C = S*N(d1) - K*e^{-rT}*N(d2)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION preco_call(S, K, r, sigma, T) RESULT(C)
      REAL(dp), INTENT(IN) :: S, K, r, sigma, T
      REAL(dp) :: d1, d2
      CALL calcular_d1_d2(S, K, r, sigma, T, d1, d2)
      C = S * normal_cumulativa(d1) - K * EXP(-r * T) * normal_cumulativa(d2)
   END FUNCTION preco_call

   !----------------------------------------------------------------------------
   ! Preço de opção de venda (PUT europeia) - Black-Scholes
   ! P = K*e^{-rT}*N(-d2) - S*N(-d1)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION preco_put(S, K, r, sigma, T) RESULT(P)
      REAL(dp), INTENT(IN) :: S, K, r, sigma, T
      REAL(dp) :: d1, d2
      CALL calcular_d1_d2(S, K, r, sigma, T, d1, d2)
      P = K * EXP(-r * T) * normal_cumulativa(-d2) - S * normal_cumulativa(-d1)
   END FUNCTION preco_put

   !----------------------------------------------------------------------------
   ! GREGAS das opções (sensibilidades)
   !----------------------------------------------------------------------------

   ! Delta (sensibilidade ao preço do ativo): dV/dS
   REAL(dp) FUNCTION delta_call(S, K, r, sigma, T) RESULT(delta)
      REAL(dp), INTENT(IN) :: S, K, r, sigma, T
      REAL(dp) :: d1, d2
      CALL calcular_d1_d2(S, K, r, sigma, T, d1, d2)
      delta = normal_cumulativa(d1)
   END FUNCTION delta_call

   REAL(dp) FUNCTION delta_put(S, K, r, sigma, T) RESULT(delta)
      REAL(dp), INTENT(IN) :: S, K, r, sigma, T
      REAL(dp) :: d1, d2
      CALL calcular_d1_d2(S, K, r, sigma, T, d1, d2)
      delta = normal_cumulativa(d1) - 1.0_dp
   END FUNCTION delta_put

   ! Gamma (convexidade): d²V/dS²  (igual para call e put)
   REAL(dp) FUNCTION gamma_opcao(S, K, r, sigma, T) RESULT(gam)
      REAL(dp), INTENT(IN) :: S, K, r, sigma, T
      REAL(dp) :: d1, d2
      CALL calcular_d1_d2(S, K, r, sigma, T, d1, d2)
      gam = densidade_normal(d1) / (S * sigma * SQRT(T) + 1.0D-15)
   END FUNCTION gamma_opcao

   ! Vega (sensibilidade à volatilidade): dV/d(sigma)
   REAL(dp) FUNCTION vega_opcao(S, K, r, sigma, T) RESULT(vega)
      REAL(dp), INTENT(IN) :: S, K, r, sigma, T
      REAL(dp) :: d1, d2
      CALL calcular_d1_d2(S, K, r, sigma, T, d1, d2)
      vega = S * densidade_normal(d1) * SQRT(T)
   END FUNCTION vega_opcao

   ! Theta (decaimento temporal): dV/dt (call)
   REAL(dp) FUNCTION theta_call(S, K, r, sigma, T) RESULT(theta)
      REAL(dp), INTENT(IN) :: S, K, r, sigma, T
      REAL(dp) :: d1, d2
      CALL calcular_d1_d2(S, K, r, sigma, T, d1, d2)
      theta = -(S * densidade_normal(d1) * sigma) / (2.0_dp * SQRT(T) + 1.0D-15) &
         - r * K * EXP(-r * T) * normal_cumulativa(d2)
   END FUNCTION theta_call

   ! Rho (sensibilidade à taxa de juros): dV/dr (call)
   REAL(dp) FUNCTION rho_call(S, K, r, sigma, T) RESULT(rho)
      REAL(dp), INTENT(IN) :: S, K, r, sigma, T
      REAL(dp) :: d1, d2
      CALL calcular_d1_d2(S, K, r, sigma, T, d1, d2)
      rho = K * T * EXP(-r * T) * normal_cumulativa(d2)
   END FUNCTION rho_call

   !----------------------------------------------------------------------------
   ! Volatilidade implícita via método de Newton-Raphson
   ! Resolve: BS(sigma) - preco_mercado = 0
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION volatilidade_implicita_call(S, K, r, T, preco_mercado) &
      RESULT(sigma_impl)
      REAL(dp), INTENT(IN) :: S, K, r, T, preco_mercado
      REAL(dp) :: sigma, bs_preco, vega_val, erro
      INTEGER :: iter
      INTEGER, PARAMETER :: MAX_ITER = 200
      REAL(dp), PARAMETER :: TOLERANCIA = 1.0D-8
      ! Estimativa inicial de Manaster-Koehler
      sigma = SQRT(ABS(2.0_dp * ABS(LOG(S/K) + r*T) / T))
      IF (sigma < 0.01_dp) sigma = 0.20_dp
      DO iter = 1, MAX_ITER
         bs_preco = preco_call(S, K, r, sigma, T)
         erro = bs_preco - preco_mercado
         IF (ABS(erro) < TOLERANCIA) EXIT
         vega_val = vega_opcao(S, K, r, sigma, T)
         IF (ABS(vega_val) < 1.0D-12) EXIT
         sigma = sigma - erro / vega_val
         IF (sigma < 0.001_dp) sigma = 0.001_dp
         IF (sigma > 10.0_dp)  sigma = 10.0_dp
      END DO
      sigma_impl = sigma
   END FUNCTION volatilidade_implicita_call

   !----------------------------------------------------------------------------
   ! Superfície de volatilidade: calcula vol. implícita para grade de K e T
   !----------------------------------------------------------------------------
   SUBROUTINE superficie_volatilidade(S, r, Ks, Ts, precos_mercado, &
      nK, nT, vol_surface)
      REAL(dp), INTENT(IN) :: S, r, Ks(nK), Ts(nT), precos_mercado(nK, nT)
      INTEGER,  INTENT(IN) :: nK, nT
      REAL(dp), INTENT(OUT) :: vol_surface(nK, nT)
      INTEGER :: i, j
      DO i = 1, nK
         DO j = 1, nT
            vol_surface(i,j) = volatilidade_implicita_call(S, Ks(i), r, &
               Ts(j), precos_mercado(i,j))
         END DO
      END DO
   END SUBROUTINE superficie_volatilidade

   !----------------------------------------------------------------------------
   ! Paridade Put-Call europeia: C - P = S - K*e^{-rT}
   ! Verifica arbitragem entre call e put
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION desvio_paridade_put_call(S, K, r, T, C_obs, P_obs) &
      RESULT(desvio)
      REAL(dp), INTENT(IN) :: S, K, r, T, C_obs, P_obs
      desvio = (C_obs - P_obs) - (S - K * EXP(-r * T))
   END FUNCTION desvio_paridade_put_call

   !----------------------------------------------------------------------------
   ! Preço de forward (contrato futuro) no modelo log-normal
   ! F = S * e^{(r - q)*T}, q = taxa de dividendos contínuos
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION preco_forward(S, r, q, T) RESULT(F)
      REAL(dp), INTENT(IN) :: S, r, q, T
      F = S * EXP((r - q) * T)
   END FUNCTION preco_forward

   !----------------------------------------------------------------------------
   ! Retorna todas as gregas em um vetor (Delta, Gamma, Vega, Theta, Rho)
   ! para call europeia
   !----------------------------------------------------------------------------
   SUBROUTINE todas_gregas_call(S, K, r, sigma, T, gregas)
      REAL(dp), INTENT(IN)  :: S, K, r, sigma, T
      REAL(dp), INTENT(OUT) :: gregas(5)
      gregas(1) = delta_call(S, K, r, sigma, T)
      gregas(2) = gamma_opcao(S, K, r, sigma, T)
      gregas(3) = vega_opcao( S, K, r, sigma, T)
      gregas(4) = theta_call(S, K, r, sigma, T)
      gregas(5) = rho_call(  S, K, r, sigma, T)
   END SUBROUTINE todas_gregas_call

END MODULE mod_black_scholes
