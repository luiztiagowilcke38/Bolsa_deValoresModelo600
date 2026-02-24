!==============================================================================
! Módulo 22: Entropia e Informação em Mercados Financeiros
! Autor: Luiz Tiago Wilcke
! Descrição: Entropia de Shannon, entropia de permutação, entropia de amostra
!            (SampEn), transferência de entropia (causalidade de informação),
!            complexidade de Lempel-Ziv e eficiência informacional para
!            análise de mercados da B3.
!==============================================================================
MODULE mod_entropia
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! Entropia de Shannon a partir de frequências relativas (histograma)
   ! H = -sum p_i * log2(p_i)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION entropia_shannon(probabilidades, n) RESULT(H)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: probabilidades(n)
      INTEGER :: i
      H = 0.0_dp
      DO i = 1, n
         IF (probabilidades(i) > 1.0D-15) &
            H = H - probabilidades(i) * LOG(probabilidades(i)) / LOG(2.0_dp)
      END DO
   END FUNCTION entropia_shannon

   !----------------------------------------------------------------------------
   ! Histograma de retornos em n_bins classes de igual largura
   !----------------------------------------------------------------------------
   SUBROUTINE histograma_retornos(retornos, n, n_bins, freq_rel)
      INTEGER,  INTENT(IN)  :: n, n_bins
      REAL(dp), INTENT(IN)  :: retornos(n)
      REAL(dp), INTENT(OUT) :: freq_rel(n_bins)
      REAL(dp) :: vmin, vmax, largura
      INTEGER :: i, idx
      INTEGER :: contagem(n_bins)
      contagem = 0
      vmin = MINVAL(retornos);  vmax = MAXVAL(retornos)
      largura = (vmax - vmin) / REAL(n_bins, dp) + 1.0D-12
      DO i = 1, n
         idx = INT((retornos(i) - vmin) / largura) + 1
         idx = MAX(1, MIN(idx, n_bins))
         contagem(idx) = contagem(idx) + 1
      END DO
      DO i = 1, n_bins
         freq_rel(i) = REAL(contagem(i), dp) / REAL(n, dp)
      END DO
   END SUBROUTINE histograma_retornos

   !----------------------------------------------------------------------------
   ! Entropia de informação dos retornos (via histograma)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION entropia_retornos(retornos, n, n_bins) RESULT(H)
      INTEGER,  INTENT(IN) :: n, n_bins
      REAL(dp), INTENT(IN) :: retornos(n)
      REAL(dp) :: freq(n_bins)
      CALL histograma_retornos(retornos, n, n_bins, freq)
      H = entropia_shannon(freq, n_bins)
   END FUNCTION entropia_retornos

   !----------------------------------------------------------------------------
   ! Entropia de Permutação (Bandt & Pompe 2002) — ordem m
   ! Divide a série em padrões de ordem m e analisa a entropia das permutações
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION entropia_permutacao(serie, n, m) RESULT(H_p)
      INTEGER,  INTENT(IN) :: n, m
      REAL(dp), INTENT(IN) :: serie(n)
      INTEGER :: n_pad, fat_m, i, j, k, idx_pad
      REAL(dp) :: temp
      INTEGER :: padroes(n), contagem_pad(720)   ! max 6! = 720 permutações
      REAL(dp) :: freq_pad(720)
      REAL(dp) :: subseq(m)
      INTEGER :: perm(m), temp_int
      ! Fatorial de m
      fat_m = 1
      DO i = 1, m
         fat_m = fat_m * i
      END DO
      n_pad = n - m + 1
      contagem_pad = 0
      DO i = 1, n_pad
         subseq = serie(i:i+m-1)
         ! Gera padrão de ordenação (permutação)
         DO j = 1, m
            perm(j) = j
         END DO
         ! Bubble sort para obter ranking
         DO j = 1, m - 1
            DO k = 1, m - j
               IF (subseq(perm(k)) > subseq(perm(k+1))) THEN
                  temp_int  = perm(k);  perm(k) = perm(k+1);  perm(k+1) = temp_int
               END IF
            END DO
         END DO
         ! Codifica padrão como índice (Lehmer code simplificado)
         idx_pad = 1
         DO j = 1, m
            idx_pad = idx_pad + (perm(j) - 1) * (m - j + 1)
         END DO
         idx_pad = MAX(1, MIN(idx_pad, fat_m))
         contagem_pad(idx_pad) = contagem_pad(idx_pad) + 1
      END DO
      DO i = 1, fat_m
         freq_pad(i) = REAL(contagem_pad(i), dp) / REAL(n_pad, dp)
      END DO
      H_p = entropia_shannon(freq_pad(1:fat_m), fat_m)
      ! Normaliza pela entropia máxima log2(m!)
      H_p = H_p / (LOG(REAL(fat_m, dp)) / LOG(2.0_dp) + 1.0D-15)
   END FUNCTION entropia_permutacao

   !----------------------------------------------------------------------------
   ! Sample Entropy (SampEn) — mede complexidade/previsibilidade da série
   ! SampEn(m, r, N) = -ln(A/B)
   ! B = pares de templates de comprimento m que correspondem
   ! A = pares de templates de comprimento m+1 que correspondem
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION entropia_amostra(serie, n, m, r_tol) RESULT(SampEn)
      INTEGER,  INTENT(IN) :: n, m
      REAL(dp), INTENT(IN) :: serie(n), r_tol
      INTEGER :: A, B, i, j, k
      LOGICAL :: match_m, match_m1
      A = 0;  B = 0
      DO i = 1, n - m
         DO j = i + 1, n - m
            match_m = .TRUE.
            DO k = 0, m - 1
               IF (ABS(serie(i+k) - serie(j+k)) > r_tol) THEN
                  match_m = .FALSE.;  EXIT
               END IF
            END DO
            IF (match_m) THEN
               B = B + 1
               IF (ABS(serie(i+m) - serie(j+m)) <= r_tol) A = A + 1
            END IF
         END DO
      END DO
      IF (B > 0 .AND. A > 0) THEN
         SampEn = -LOG(REAL(A, dp) / REAL(B, dp))
      ELSE
         SampEn = 0.0_dp
      END IF
   END FUNCTION entropia_amostra

   !----------------------------------------------------------------------------
   ! Eficiência de Mercado (baseada na entropia)
   ! EME = H_observado / H_maximo (onde H_max = log2(n_bins))
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION eficiencia_mercado(H_obs, n_bins) RESULT(EME)
      REAL(dp), INTENT(IN) :: H_obs
      INTEGER,  INTENT(IN) :: n_bins
      REAL(dp) :: H_max
      H_max = LOG(REAL(n_bins, dp)) / LOG(2.0_dp)
      IF (H_max > 1.0D-10) THEN
         EME = H_obs / H_max
      ELSE
         EME = 0.0_dp
      END IF
   END FUNCTION eficiencia_mercado

   !----------------------------------------------------------------------------
   ! Informação Mútua entre dois ativos: I(X;Y) = H(X) + H(Y) - H(X,Y)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION informacao_mutua(ret_x, ret_y, n, n_bins) RESULT(IM)
      INTEGER,  INTENT(IN) :: n, n_bins
      REAL(dp), INTENT(IN) :: ret_x(n), ret_y(n)
      REAL(dp) :: H_X, H_Y, freq_xy(n_bins*n_bins), freq_x(n_bins), freq_y(n_bins)
      REAL(dp) :: vmin_x, vmax_x, vmin_y, vmax_y, larg_x, larg_y, H_XY
      INTEGER :: i, ix, iy, idx_joint
      CALL histograma_retornos(ret_x, n, n_bins, freq_x)
      CALL histograma_retornos(ret_y, n, n_bins, freq_y)
      H_X = entropia_shannon(freq_x, n_bins)
      H_Y = entropia_shannon(freq_y, n_bins)
      vmin_x = MINVAL(ret_x);  vmax_x = MAXVAL(ret_x)
      vmin_y = MINVAL(ret_y);  vmax_y = MAXVAL(ret_y)
      larg_x = (vmax_x - vmin_x) / REAL(n_bins, dp) + 1.0D-12
      larg_y = (vmax_y - vmin_y) / REAL(n_bins, dp) + 1.0D-12
      freq_xy = 0.0_dp
      DO i = 1, n
         ix = MAX(1, MIN(n_bins, INT((ret_x(i) - vmin_x) / larg_x) + 1))
         iy = MAX(1, MIN(n_bins, INT((ret_y(i) - vmin_y) / larg_y) + 1))
         idx_joint = (ix - 1) * n_bins + iy
         freq_xy(idx_joint) = freq_xy(idx_joint) + 1.0_dp / REAL(n, dp)
      END DO
      H_XY = entropia_shannon(freq_xy, n_bins*n_bins)
      IM = H_X + H_Y - H_XY
   END FUNCTION informacao_mutua

END MODULE mod_entropia
