!==============================================================================
! Módulo 23: Geometria Fractal de Mercados
! Autor: Luiz Tiago Wilcke
! Descrição: Dimensão fractal de Hausdorff, expoente de Hurst (R/S),
!            DFA (Detrended Fluctuation Analysis), multifractais e
!            análise de persistência/antipersistência para mercados da B3.
!==============================================================================
MODULE mod_fractal
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! Análise R/S de Hurst
   ! Parte a série em sub-séries de comprimento n, calcula R/S para cada
   ! Expoente de Hurst: E[R/S] ~ c * n^H
   ! H > 0.5: persistente (tendência), H < 0.5: antipersistente (reversão)
   !----------------------------------------------------------------------------
   SUBROUTINE expoente_hurst_rs(serie, n, n_escalas, H_hurst)
      INTEGER,  INTENT(IN)  :: n, n_escalas
      REAL(dp), INTENT(IN)  :: serie(n)
      REAL(dp), INTENT(OUT) :: H_hurst
      REAL(dp) :: escala, log_rs_medio, log_escala_medio
      REAL(dp) :: soma_x, soma_y, soma_xy, soma_x2
      REAL(dp) :: log_escala(n_escalas), log_rs(n_escalas)
      INTEGER :: k, n_esc, m_sub, n_sub, i_sub, j
      REAL(dp) :: sub_serie(n), media_sub, desvio_acum(n), range_sub, sigma_sub
      DO k = 1, n_escalas
         n_esc = INT(REAL(n, dp) * REAL(k, dp) / REAL(n_escalas, dp))
         IF (n_esc < 4) n_esc = 4
         n_sub  = n / n_esc   ! número de sub-séries
         IF (n_sub < 1) n_sub = 1
         log_rs_medio = 0.0_dp
         DO m_sub = 1, n_sub
            i_sub = (m_sub - 1) * n_esc + 1
            sub_serie(1:n_esc) = serie(i_sub:MIN(i_sub+n_esc-1, n))
            media_sub = SUM(sub_serie(1:n_esc)) / REAL(n_esc, dp)
            ! Série de desvios acumulados
            desvio_acum(1) = sub_serie(1) - media_sub
            DO j = 2, n_esc
               desvio_acum(j) = desvio_acum(j-1) + sub_serie(j) - media_sub
            END DO
            range_sub = MAXVAL(desvio_acum(1:n_esc)) - MINVAL(desvio_acum(1:n_esc))
            sigma_sub = SQRT(SUM((sub_serie(1:n_esc) - media_sub)**2) / REAL(n_esc, dp))
            IF (sigma_sub > 1.0D-15) &
               log_rs_medio = log_rs_medio + LOG(range_sub / sigma_sub)
         END DO
         log_rs(k)     = log_rs_medio / REAL(n_sub, dp)
         log_escala(k) = LOG(REAL(n_esc, dp))
      END DO
      ! Ajuste linear: log(R/S) = H * log(n) + c
      soma_x = SUM(log_escala(1:n_escalas));  soma_y = SUM(log_rs(1:n_escalas))
      soma_xy = DOT_PRODUCT(log_escala(1:n_escalas), log_rs(1:n_escalas))
      soma_x2 = SUM(log_escala(1:n_escalas)**2)
      H_hurst = (REAL(n_escalas,dp)*soma_xy - soma_x*soma_y) / &
         (REAL(n_escalas,dp)*soma_x2 - soma_x**2 + 1.0D-15)
   END SUBROUTINE expoente_hurst_rs

   !----------------------------------------------------------------------------
   ! DFA — Detrended Fluctuation Analysis
   ! Calcula F(s) para cada escala s e estima expoente alpha via regressão
   !----------------------------------------------------------------------------
   SUBROUTINE dfa(serie, n, n_min, n_max, alpha_dfa)
      INTEGER,  INTENT(IN)  :: n, n_min, n_max
      REAL(dp), INTENT(IN)  :: serie(n)
      REAL(dp), INTENT(OUT) :: alpha_dfa
      REAL(dp) :: integral(n), media_glob
      REAL(dp) :: F2, y_fit, b0, b1, r2_dummy, ep_dummy
      REAL(dp) :: log_s_vals(50), log_F_vals(50)
      REAL(dp) :: sx, sy, sxy, sx2
      INTEGER :: s, n_box, i_box, n_seg, k, i, n_pts
      ! Série integrada: Y_i = sum_{j=1}^{i} (x_j - mean)
      media_glob = SUM(serie) / REAL(n, dp)
      integral(1) = serie(1) - media_glob
      DO i = 2, n
         integral(i) = integral(i-1) + serie(i) - media_glob
      END DO
      n_pts = 0
      s = n_min
      DO WHILE (s <= n_max .AND. n_pts < 50)
         n_seg = n / s
         IF (n_seg < 1) EXIT
         F2 = 0.0_dp
         DO i_box = 0, n_seg - 1
            ! Detrend por OLS linear dentro do segmento
            sx = 0.0_dp; sy = 0.0_dp; sxy = 0.0_dp; sx2 = 0.0_dp
            DO k = 1, s
               i = i_box * s + k
               sx  = sx  + REAL(k, dp)
               sy  = sy  + integral(i)
               sxy = sxy + REAL(k, dp) * integral(i)
               sx2 = sx2 + REAL(k, dp)**2
            END DO
            b1 = (REAL(s,dp)*sxy - sx*sy) / (REAL(s,dp)*sx2 - sx**2 + 1.0D-15)
            b0 = (sy - b1*sx) / REAL(s, dp)
            DO k = 1, s
               i = i_box * s + k
               y_fit = b0 + b1 * REAL(k, dp)
               F2 = F2 + (integral(i) - y_fit)**2
            END DO
         END DO
         n_pts = n_pts + 1
         log_s_vals(n_pts) = LOG(REAL(s, dp))
         log_F_vals(n_pts) = 0.5_dp * LOG(MAX(F2 / REAL(n_seg * s, dp), 1.0D-30))
         s = s * 2
      END DO
      ! Regressão log-log
      IF (n_pts > 1) THEN
         sx = SUM(log_s_vals(1:n_pts)); sy = SUM(log_F_vals(1:n_pts))
         sxy = DOT_PRODUCT(log_s_vals(1:n_pts), log_F_vals(1:n_pts))
         sx2 = SUM(log_s_vals(1:n_pts)**2)
         alpha_dfa = (REAL(n_pts,dp)*sxy - sx*sy) / &
            (REAL(n_pts,dp)*sx2 - sx**2 + 1.0D-15)
      ELSE
         alpha_dfa = 0.5_dp
      END IF
   END SUBROUTINE dfa

   !----------------------------------------------------------------------------
   ! Dimensão fractal de Higuchi — mede complexidade da série
   ! D_H = -dlog(L_m(k)) / dlog(k)   geralmente 1 < D_H < 2
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION dimensao_higuchi(serie, n, k_max) RESULT(D_H)
      INTEGER,  INTENT(IN) :: n, k_max
      REAL(dp), INTENT(IN) :: serie(n)
      REAL(dp) :: L_k(k_max), log_k(k_max), log_L(k_max)
      REAL(dp) :: sx, sy, sxy, sx2, L_m, soma
      INTEGER :: k, m, i_max, i
      DO k = 1, k_max
         L_k(k) = 0.0_dp
         DO m = 1, k
            i_max = INT((n - m) / k)
            soma  = 0.0_dp
            DO i = 1, i_max
               soma = soma + ABS(serie(m + i*k) - serie(m + (i-1)*k))
            END DO
            IF (i_max > 0) &
               L_m = soma * REAL(n - 1, dp) / REAL(k * i_max, dp) / REAL(k, dp)
            L_k(k) = L_k(k) + L_m
         END DO
         L_k(k) = L_k(k) / REAL(k, dp)
         log_k(k) = LOG(REAL(k, dp))
         log_L(k) = LOG(MAX(L_k(k), 1.0D-30))
      END DO
      sx = SUM(log_k); sy = SUM(log_L)
      sxy = DOT_PRODUCT(log_k, log_L); sx2 = SUM(log_k**2)
      D_H = -(REAL(k_max,dp)*sxy - sx*sy) / (REAL(k_max,dp)*sx2 - sx**2 + 1.0D-15)
   END FUNCTION dimensao_higuchi

   !----------------------------------------------------------------------------
   ! Classifica eficiência pelo Expoente de Hurst
   !----------------------------------------------------------------------------
   CHARACTER(LEN=20) FUNCTION interpretar_hurst(H) RESULT(classificacao)
      REAL(dp), INTENT(IN) :: H
      IF (H > 0.6_dp) THEN
         classificacao = "Persistente"
      ELSE IF (H < 0.4_dp) THEN
         classificacao = "Antipersistente"
      ELSE
         classificacao = "Ruido Browniano"
      END IF
   END FUNCTION interpretar_hurst

END MODULE mod_fractal
