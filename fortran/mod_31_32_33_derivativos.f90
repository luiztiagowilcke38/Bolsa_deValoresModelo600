!==============================================================================
! Módulo 31-32-33: Derivativos, Opções Exóticas e Arbitragem
! Autor: Luiz Tiago Wilcke
! Descrição: Precificação de derivativos complexos: futuros DI, swaps de
!            taxa de juros (IRS), opções exóticas (asiáticas, barreira,
!            lookback, digitais) e detecção de oportunidades de arbitragem
!            no mercado de derivativos da B3.
!==============================================================================
MODULE mod_derivativos
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)
   REAL(dp), PARAMETER :: PI = 3.14159265358979323846_dp

CONTAINS

   !----------------------------------------------------------------------------
   ! Preço de contrato futuro DI (taxa de juros brasileira)
   ! F_DI = 100000 / (1 + taxa_DI)^{n/252}
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION preco_futuro_di(taxa_di, n_dias_uteis) RESULT(F_DI)
      REAL(dp), INTENT(IN) :: taxa_di
      INTEGER,  INTENT(IN) :: n_dias_uteis
      F_DI = 100000.0_dp / (1.0_dp + taxa_di) ** (REAL(n_dias_uteis, dp) / 252.0_dp)
   END FUNCTION preco_futuro_di

   !----------------------------------------------------------------------------
   ! Valor Presente de um Swap de Taxa de Juros (IRS: pré vs CDI)
   ! VP_perna_fixa = N * taxa_fixa * sum(P(0,Ti))
   ! VP_perna_flutuante = N * (P(0,T0) - P(0,Tn))
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION valor_swap_irs(N_nominal, taxa_fixa, P_zcb, n_periodos, &
      P_t0, P_tn) RESULT(VP_swap)
      REAL(dp), INTENT(IN) :: N_nominal, taxa_fixa, P_zcb(n_periodos)
      REAL(dp), INTENT(IN) :: P_t0, P_tn
      INTEGER,  INTENT(IN) :: n_periodos
      REAL(dp) :: VP_fixa, VP_flutuante
      VP_fixa      = N_nominal * taxa_fixa * SUM(P_zcb)
      VP_flutuante = N_nominal * (P_t0 - P_tn)
      VP_swap      = VP_fixa - VP_flutuante   ! perspectiva do pagador da taxa fixa
   END FUNCTION valor_swap_irs

   !============ OPÇÕES EXÓTICAS ===============================================

   !----------------------------------------------------------------------------
   ! Opção Asiática (média aritmética) via Monte Carlo
   ! Payoff_call = max(A_T - K, 0);  A_T = (1/n) * sum S_t
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION call_asiatica_mc(S0, K, r, sigma, T, n_passos, n_traj) &
      RESULT(C_asia)
      REAL(dp), INTENT(IN) :: S0, K, r, sigma, T
      INTEGER,  INTENT(IN) :: n_passos, n_traj
      REAL(dp) :: S, media_traj, soma_payoff, dt, drift, vol, u1, u2, z
      INTEGER :: i, j
      dt    = T / REAL(n_passos, dp)
      drift = (r - 0.5_dp * sigma**2) * dt
      vol   = sigma * SQRT(dt)
      soma_payoff = 0.0_dp
      DO j = 1, n_traj
         S = S0;  media_traj = 0.0_dp
         DO i = 1, n_passos
            CALL RANDOM_NUMBER(u1);  CALL RANDOM_NUMBER(u2)
            u1 = MAX(u1, 1.0D-15)
            z  = SQRT(-2.0_dp*LOG(u1)) * COS(2.0_dp*PI*u2)
            S  = S * EXP(drift + vol * z)
            media_traj = media_traj + S
         END DO
         media_traj  = media_traj / REAL(n_passos, dp)
         soma_payoff = soma_payoff + MAX(media_traj - K, 0.0_dp)
      END DO
      C_asia = EXP(-r * T) * soma_payoff / REAL(n_traj, dp)
   END FUNCTION call_asiatica_mc

   !----------------------------------------------------------------------------
   ! Opção de Barreira (Down-and-Out Call) via Monte Carlo
   ! Payoff = max(S_T - K, 0) * I(S_t > B para todo t em [0,T])
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION call_barreira_down_out_mc(S0, K, B, r, sigma, T, &
      n_passos, n_traj) RESULT(C_bar)
      REAL(dp), INTENT(IN) :: S0, K, B, r, sigma, T
      INTEGER,  INTENT(IN) :: n_passos, n_traj
      REAL(dp) :: S, soma_payoff, dt, drift, vol, u1, u2, z
      LOGICAL :: toque_barreira
      INTEGER :: i, j
      dt    = T / REAL(n_passos, dp)
      drift = (r - 0.5_dp * sigma**2) * dt
      vol   = sigma * SQRT(dt)
      soma_payoff = 0.0_dp
      DO j = 1, n_traj
         S = S0;  toque_barreira = .FALSE.
         DO i = 1, n_passos
            CALL RANDOM_NUMBER(u1);  CALL RANDOM_NUMBER(u2)
            u1 = MAX(u1, 1.0D-15)
            z  = SQRT(-2.0_dp*LOG(u1)) * COS(2.0_dp*PI*u2)
            S  = S * EXP(drift + vol * z)
            IF (S <= B) THEN
               toque_barreira = .TRUE.;  EXIT
            END IF
         END DO
         IF (.NOT. toque_barreira) &
            soma_payoff = soma_payoff + MAX(S - K, 0.0_dp)
      END DO
      C_bar = EXP(-r * T) * soma_payoff / REAL(n_traj, dp)
   END FUNCTION call_barreira_down_out_mc

   !----------------------------------------------------------------------------
   ! Opção Digital (Cash-or-Nothing Call)
   ! Payoff = Q se S_T > K, ou 0 caso contrário
   ! C_dig = Q * e^{-rT} * N(d2)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION call_digital(S, K, r, sigma, T, Q) RESULT(C_dig)
      REAL(dp), INTENT(IN) :: S, K, r, sigma, T, Q
      REAL(dp) :: d2, N_d2
      d2  = (LOG(S/K) + (r - 0.5_dp*sigma**2)*T) / (sigma*SQRT(T) + 1.0D-15)
      N_d2 = 0.5_dp * (1.0_dp + ERF(d2 / SQRT(2.0_dp)))
      C_dig = Q * EXP(-r * T) * N_d2
   END FUNCTION call_digital

   !============ ARBITRAGEM ====================================================

   !----------------------------------------------------------------------------
   ! Detecção de arbitragem de paridade put-call
   ! Arbitragem existe se |C - P - S + K*e^{-rT}| > bid-ask spread
   !----------------------------------------------------------------------------
   LOGICAL FUNCTION arbitragem_paridade_pc(C, P, S, K, r, T, spread_bid_ask)
      REAL(dp), INTENT(IN) :: C, P, S, K, r, T, spread_bid_ask
      REAL(dp) :: desvio
      desvio = ABS(C - P - S + K * EXP(-r * T))
      arbitragem_paridade_pc = (desvio > spread_bid_ask)
   END FUNCTION arbitragem_paridade_pc

   !----------------------------------------------------------------------------
   ! Estratégia Box Spread (arbitragem sem risco em opções)
   ! Box = (K2-K1)*e^{-rT}  (sempre)
   ! Lucro = Box_mercado - Box_teórico (se diferente)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION lucro_box_spread(C1, P1, C2, P2, K1, K2, r, T) &
      RESULT(lucro)
      REAL(dp), INTENT(IN) :: C1, P1, C2, P2, K1, K2, r, T
      REAL(dp) :: custo_box, valor_teorico
      custo_box    = (C1 - C2) - (P1 - P2)
      valor_teorico = (K2 - K1) * EXP(-r * T)
      lucro = valor_teorico - custo_box
   END FUNCTION lucro_box_spread

   !----------------------------------------------------------------------------
   ! Função de erro para uso interno (ERF aproximado)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION ERF(x) RESULT(erf_val)
      REAL(dp), INTENT(IN) :: x
      REAL(dp) :: t, abs_x
      abs_x = ABS(x)
      t = 1.0_dp / (1.0_dp + 0.3275911_dp * abs_x)
      erf_val = 1.0_dp - (0.254829592_dp + t*(-0.284496736_dp + &
         t*(1.421413741_dp + t*(-1.453152027_dp + t*1.061405429_dp)))) * &
         t * EXP(-abs_x**2)
      IF (x < 0.0_dp) erf_val = -erf_val
   END FUNCTION ERF

END MODULE mod_derivativos
