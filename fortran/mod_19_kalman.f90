!==============================================================================
! Módulo 19: Filtro de Kalman para Rastreamento de Tendência
! Autor: Luiz Tiago Wilcke
! Descrição: Filtro de Kalman linear e não-linear (EKF/UKF simplificado)
!            para estimação do estado latente (tendência, nível) de preços
!            e spreads de juros na B3. Inclui filtro suavizador de Kalman.
!==============================================================================
MODULE mod_kalman
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! Filtro de Kalman escalar (modelo AR(1) para estado latente)
   ! Estado: x_t = phi*x_{t-1} + w_t,    w_t ~ N(0, Q)
   ! Observação: y_t = x_t + v_t,        v_t ~ N(0, R)
   !----------------------------------------------------------------------------
   SUBROUTINE filtro_kalman_escalar(observacoes, phi, Q, R, x0, P0, &
      estados_filtrados, P_filtrado)
      REAL(dp), INTENT(IN)  :: observacoes(:), phi, Q, R, x0, P0
      REAL(dp), INTENT(OUT) :: estados_filtrados(:), P_filtrado(:)
      REAL(dp) :: x_prev, P_prev, x_pred, P_pred, K_gain, inovacao
      INTEGER :: n, t
      n = SIZE(observacoes)
      x_prev = x0;  P_prev = P0
      DO t = 1, n
         ! Predição:  x_{t|t-1} = phi * x_{t-1|t-1}
         x_pred = phi * x_prev
         P_pred = phi**2 * P_prev + Q
         ! Ganho de Kalman: K_t = P_{t|t-1} / (P_{t|t-1} + R)
         K_gain = P_pred / (P_pred + R + 1.0D-15)
         ! Atualização:
         inovacao = observacoes(t) - x_pred
         x_prev = x_pred + K_gain * inovacao
         P_prev = (1.0_dp - K_gain) * P_pred
         estados_filtrados(t) = x_prev
         P_filtrado(t)         = P_prev
      END DO
   END SUBROUTINE filtro_kalman_escalar

   !----------------------------------------------------------------------------
   ! Filtro de Kalman multivariado (modelo espaço-estado geral)
   ! Estado: x_t = F*x_{t-1} + B*u_t + w_t   w_t ~ N(0, Q)
   ! Obs:    y_t = H*x_t + v_t                v_t ~ N(0, R)
   !----------------------------------------------------------------------------
   SUBROUTINE filtro_kalman_mv(observacoes, F, H, Q_mat, R_mat, x0, P0, &
      n_obs, n_estado, n_obs_dim, &
      x_filtrado, P_filtrado_diag)
      INTEGER,  INTENT(IN)  :: n_obs, n_estado, n_obs_dim
      REAL(dp), INTENT(IN)  :: observacoes(n_obs, n_obs_dim)
      REAL(dp), INTENT(IN)  :: F(n_estado, n_estado), H(n_obs_dim, n_estado)
      REAL(dp), INTENT(IN)  :: Q_mat(n_estado, n_estado), R_mat(n_obs_dim, n_obs_dim)
      REAL(dp), INTENT(IN)  :: x0(n_estado), P0(n_estado, n_estado)
      REAL(dp), INTENT(OUT) :: x_filtrado(n_obs, n_estado)
      REAL(dp), INTENT(OUT) :: P_filtrado_diag(n_obs, n_estado)
      REAL(dp) :: x_k(n_estado), P_k(n_estado, n_estado)
      REAL(dp) :: x_pred(n_estado), P_pred(n_estado, n_estado)
      REAL(dp) :: S_mat(n_obs_dim, n_obs_dim), K_mat(n_estado, n_obs_dim)
      REAL(dp) :: HP(n_obs_dim, n_estado), S_inv(n_obs_dim, n_obs_dim)
      REAL(dp) :: inovacao(n_obs_dim), KH(n_estado, n_estado), I_mat(n_estado, n_estado)
      INTEGER :: t, i
      ! Identidade
      I_mat = 0.0_dp
      DO i = 1, n_estado
         I_mat(i,i) = 1.0_dp
      END DO
      x_k = x0;  P_k = P0
      DO t = 1, n_obs
         ! Predição
         x_pred = MATMUL(F, x_k)
         P_pred = MATMUL(MATMUL(F, P_k), TRANSPOSE(F)) + Q_mat
         ! Inovação
         inovacao = observacoes(t,:) - MATMUL(H, x_pred)
         HP  = MATMUL(H, P_pred)
         S_mat = MATMUL(HP, TRANSPOSE(H)) + R_mat
         ! Inversão simples (diagonal para eficiência quando n_obs_dim=1)
         IF (n_obs_dim == 1) THEN
            S_inv(1,1) = 1.0_dp / MAX(S_mat(1,1), 1.0D-15)
         ELSE
            S_inv = S_mat   ! Simplificado — em produção usar LAPACK
            DO i = 1, n_obs_dim
               IF (ABS(S_inv(i,i)) > 1.0D-15) S_inv(i,i) = 1.0_dp / S_inv(i,i)
            END DO
         END IF
         ! Ganho de Kalman: K = P_pred * H' * S^{-1}
         K_mat = MATMUL(MATMUL(P_pred, TRANSPOSE(H)), S_inv)
         ! Atualização
         x_k = x_pred + MATMUL(K_mat, inovacao)
         KH  = MATMUL(K_mat, H)
         P_k = MATMUL(I_mat - KH, P_pred)
         x_filtrado(t,:) = x_k
         DO i = 1, n_estado
            P_filtrado_diag(t, i) = P_k(i, i)
         END DO
      END DO
   END SUBROUTINE filtro_kalman_mv

   !----------------------------------------------------------------------------
   ! Suavizador de Kalman (Rauch-Tung-Striebel) — pass backward
   ! x_{t|T} = x_{t|t} + G_t*(x_{t+1|T} - F*x_{t|t})
   ! G_t = P_{t|t}*F'*P_{t+1|t}^{-1}
   !----------------------------------------------------------------------------
   SUBROUTINE suavizador_kalman(x_filtrado, P_diag, F, n, n_estado, x_suavizado)
      INTEGER,  INTENT(IN)  :: n, n_estado
      REAL(dp), INTENT(IN)  :: x_filtrado(n, n_estado), P_diag(n, n_estado)
      REAL(dp), INTENT(IN)  :: F(n_estado, n_estado)
      REAL(dp), INTENT(OUT) :: x_suavizado(n, n_estado)
      REAL(dp) :: G(n_estado), diff(n_estado)
      INTEGER :: t, i
      x_suavizado(n,:) = x_filtrado(n,:)
      DO t = n-1, 1, -1
         DO i = 1, n_estado
            IF (ABS(P_diag(t+1, i)) > 1.0D-15) THEN
               G(i) = P_diag(t, i) * F(i,i) / P_diag(t+1, i)
            ELSE
               G(i) = 0.0_dp
            END IF
            diff(i) = x_suavizado(t+1,i) - DOT_PRODUCT(F(i,:), x_filtrado(t,:))
            x_suavizado(t,i) = x_filtrado(t,i) + G(i) * diff(i)
         END DO
      END DO
   END SUBROUTINE suavizador_kalman

   !----------------------------------------------------------------------------
   ! Kalman adaptativo: ajusta Q e R via critério de inovação
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION ajustar_ruido_processo(inovacao_quadratica, S_escalar, &
      Q_atual) RESULT(Q_novo)
      REAL(dp), INTENT(IN) :: inovacao_quadratica, S_escalar, Q_atual
      Q_novo = Q_atual + 0.1_dp * (inovacao_quadratica - S_escalar)
      Q_novo = MAX(Q_novo, 1.0D-8)
   END FUNCTION ajustar_ruido_processo

END MODULE mod_kalman
