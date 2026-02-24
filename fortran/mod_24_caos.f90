!==============================================================================
! Módulo 24: Teoria do Caos em Mercados Financeiros
! Autor: Luiz Tiago Wilcke
! Descrição: Expoente de Lyapunov, dimensão de correlação de Grassberger-
!            Procaccia, atratores de Lorenz financeiro, reconstrução de
!            espaço de fases (teorema de Takens) e análise de predibilidade
!            caótica para séries financeiras da B3.
!==============================================================================
MODULE mod_caos
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! Reconstrução do espaço de fases pelo teorema de Takens
   ! y_i = [x_i, x_{i+tau}, x_{i+2*tau}, ..., x_{i+(m-1)*tau}]
   !----------------------------------------------------------------------------
   SUBROUTINE espaco_fases(serie, n, m_embed, tau, Y, n_Y)
      INTEGER,  INTENT(IN)  :: n, m_embed, tau
      REAL(dp), INTENT(IN)  :: serie(n)
      REAL(dp), INTENT(OUT) :: Y(n, m_embed)
      INTEGER,  INTENT(OUT) :: n_Y
      INTEGER :: i, j
      n_Y = n - (m_embed - 1) * tau
      IF (n_Y < 1) THEN
         n_Y = 0; RETURN
      END IF
      DO i = 1, n_Y
         DO j = 1, m_embed
            Y(i, j) = serie(i + (j-1) * tau)
         END DO
      END DO
   END SUBROUTINE espaco_fases

   !----------------------------------------------------------------------------
   ! Expoente de Lyapunov máximo via algoritmo de Rosenstein
   ! lambda = (1/t) * <ln|dJ_t|>
   ! Mede taxa de separação de trajetórias próximas (positivo = caos)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION expoente_lyapunov(serie, n, m_embed, tau, dt) RESULT(lambda)
      INTEGER,  INTENT(IN) :: n, m_embed, tau
      REAL(dp), INTENT(IN) :: serie(n), dt
      REAL(dp) :: Y(n, m_embed), dist_inicial, dist_t, soma_log
      INTEGER :: n_Y, i, j, vizinho, idx_viz
      REAL(dp) :: dist, dist_min
      CALL espaco_fases(serie, n, m_embed, tau, Y, n_Y)
      IF (n_Y < 4) THEN
         lambda = 0.0_dp;  RETURN
      END IF
      soma_log = 0.0_dp
      DO i = 1, n_Y - 1
         dist_min = HUGE(1.0_dp);  idx_viz = 0
         DO j = 1, n_Y - 1
            IF (ABS(i - j) > 1) THEN   ! evita vizinhos temporais
               dist = SQRT(SUM((Y(i,:) - Y(j,:))**2))
               IF (dist < dist_min) THEN
                  dist_min = dist;  idx_viz = j
               END IF
            END IF
         END DO
         IF (idx_viz > 0 .AND. idx_viz < n_Y .AND. i < n_Y) THEN
            dist_t = SQRT(SUM((Y(i+1,:) - Y(idx_viz+1,:))**2))
            IF (dist_min > 1.0D-15 .AND. dist_t > 0.0_dp) &
               soma_log = soma_log + LOG(dist_t / dist_min)
         END IF
      END DO
      lambda = soma_log / (REAL(n_Y - 1, dp) * dt)
   END FUNCTION expoente_lyapunov

   !----------------------------------------------------------------------------
   ! Integral de correlação de Grassberger-Procaccia
   ! C(r) = (2/N^2) * count(|y_i - y_j| < r)
   !----------------------------------------------------------------------------
   SUBROUTINE integral_correlacao(Y, n_Y, m_embed, r_vals, n_r, C_r)
      INTEGER,  INTENT(IN)  :: n_Y, m_embed, n_r
      REAL(dp), INTENT(IN)  :: Y(n_Y, m_embed), r_vals(n_r)
      REAL(dp), INTENT(OUT) :: C_r(n_r)
      INTEGER :: i, j, k
      REAL(dp) :: dist
      C_r = 0.0_dp
      DO i = 1, n_Y
         DO j = i+1, n_Y
            dist = SQRT(SUM((Y(i,:) - Y(j,:))**2))
            DO k = 1, n_r
               IF (dist < r_vals(k)) C_r(k) = C_r(k) + 1.0_dp
            END DO
         END DO
      END DO
      C_r = 2.0_dp * C_r / REAL(n_Y * (n_Y - 1), dp)
   END SUBROUTINE integral_correlacao

   !----------------------------------------------------------------------------
   ! Dimensão de correlação (inclinação de log(C(r)) vs log(r))
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION dimensao_correlacao(C_r, r_vals, n_r) RESULT(D_corr)
      INTEGER,  INTENT(IN) :: n_r
      REAL(dp), INTENT(IN) :: C_r(n_r), r_vals(n_r)
      REAL(dp) :: sx, sy, sxy, sx2, log_r, log_c
      INTEGER :: k, n_pts
      sx = 0.0_dp; sy = 0.0_dp; sxy = 0.0_dp; sx2 = 0.0_dp; n_pts = 0
      DO k = 1, n_r
         IF (C_r(k) > 1.0D-15 .AND. r_vals(k) > 1.0D-15) THEN
            log_r = LOG(r_vals(k))
            log_c = LOG(C_r(k))
            sx = sx + log_r; sy = sy + log_c
            sxy = sxy + log_r * log_c; sx2 = sx2 + log_r**2
            n_pts = n_pts + 1
         END IF
      END DO
      IF (n_pts > 1) THEN
         D_corr = (REAL(n_pts,dp)*sxy - sx*sy) / (REAL(n_pts,dp)*sx2 - sx**2 + 1.0D-15)
      ELSE
         D_corr = 0.0_dp
      END IF
   END FUNCTION dimensao_correlacao

   !----------------------------------------------------------------------------
   ! Mapa Logístico Financeiro — demonstração de comportamento caótico
   ! x_{t+1} = r * x_t * (1 - x_t)
   ! Para r > 3.57: caos
   !----------------------------------------------------------------------------
   SUBROUTINE mapa_logistico(x0, r_param, n_iter, trajetoria)
      REAL(dp), INTENT(IN)  :: x0, r_param
      INTEGER,  INTENT(IN)  :: n_iter
      REAL(dp), INTENT(OUT) :: trajetoria(n_iter+1)
      INTEGER :: i
      trajetoria(1) = x0
      DO i = 1, n_iter
         trajetoria(i+1) = r_param * trajetoria(i) * (1.0_dp - trajetoria(i))
      END DO
   END SUBROUTINE mapa_logistico

   !----------------------------------------------------------------------------
   ! Análise de recorrência — matriz booleana de recorrência em espaço de fases
   ! R(i,j) = 1 se |y_i - y_j| < epsilon
   !----------------------------------------------------------------------------
   SUBROUTINE matriz_recorrencia(Y, n_Y, m_embed, epsilon, R)
      INTEGER,  INTENT(IN)  :: n_Y, m_embed
      REAL(dp), INTENT(IN)  :: Y(n_Y, m_embed), epsilon
      LOGICAL,  INTENT(OUT) :: R(n_Y, n_Y)
      REAL(dp) :: dist
      INTEGER :: i, j
      DO i = 1, n_Y
         DO j = 1, n_Y
            dist = SQRT(SUM((Y(i,:) - Y(j,:))**2))
            R(i,j) = (dist < epsilon)
         END DO
      END DO
   END SUBROUTINE matriz_recorrencia

   !----------------------------------------------------------------------------
   ! Taxa de recorrência (fração de pontos na caixa epsilon)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION taxa_recorrencia(R, n_Y) RESULT(RR)
      INTEGER, INTENT(IN) :: n_Y
      LOGICAL, INTENT(IN) :: R(n_Y, n_Y)
      INTEGER :: i, j, count_R
      count_R = 0
      DO i = 1, n_Y
         DO j = 1, n_Y
            IF (R(i,j)) count_R = count_R + 1
         END DO
      END DO
      RR = REAL(count_R, dp) / REAL(n_Y * n_Y, dp)
   END FUNCTION taxa_recorrencia

END MODULE mod_caos
