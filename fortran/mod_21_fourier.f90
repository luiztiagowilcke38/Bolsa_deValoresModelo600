!==============================================================================
! Módulo 21: Análise de Fourier para Ciclos de Mercado
! Autor: Luiz Tiago Wilcke
! Descrição: Transformada Discreta de Fourier (DFT), análise espectral
!            de séries financeiras, detecção de ciclos dominantes,
!            filtragem por banda e reconstrução sintética de preços
!            para identificação de sazonalidades no mercado brasileiro.
!==============================================================================
MODULE mod_fourier
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)
   REAL(dp), PARAMETER :: PI = 3.14159265358979323846_dp

CONTAINS

   !----------------------------------------------------------------------------
   ! DFT (Transformada Discreta de Fourier) via definição direta O(N^2)
   ! X_k = sum_{n=0}^{N-1} x_n * exp(-2*pi*i*k*n/N)
   ! Armazena parte real e imaginária separadamente
   !----------------------------------------------------------------------------
   SUBROUTINE dft(x, n, X_real, X_imag)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: x(n)
      REAL(dp), INTENT(OUT) :: X_real(n), X_imag(n)
      INTEGER :: k, j
      REAL(dp) :: angulo
      DO k = 0, n-1
         X_real(k+1) = 0.0_dp;  X_imag(k+1) = 0.0_dp
         DO j = 0, n-1
            angulo = -2.0_dp * PI * REAL(k*j, dp) / REAL(n, dp)
            X_real(k+1) = X_real(k+1) + x(j+1) * COS(angulo)
            X_imag(k+1) = X_imag(k+1) + x(j+1) * SIN(angulo)
         END DO
      END DO
   END SUBROUTINE dft

   !----------------------------------------------------------------------------
   ! IDFT (Transformada Inversa Discreta de Fourier)
   ! x_n = (1/N) * sum_{k=0}^{N-1} X_k * exp(2*pi*i*k*n/N)
   !----------------------------------------------------------------------------
   SUBROUTINE idft(X_real, X_imag, n, x_reconstruido)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: X_real(n), X_imag(n)
      REAL(dp), INTENT(OUT) :: x_reconstruido(n)
      INTEGER :: k, j
      REAL(dp) :: soma, angulo
      DO j = 0, n-1
         soma = 0.0_dp
         DO k = 0, n-1
            angulo = 2.0_dp * PI * REAL(k*j, dp) / REAL(n, dp)
            soma = soma + X_real(k+1)*COS(angulo) - X_imag(k+1)*SIN(angulo)
         END DO
         x_reconstruido(j+1) = soma / REAL(n, dp)
      END DO
   END SUBROUTINE idft

   !----------------------------------------------------------------------------
   ! Módulo (amplitude espectral): |X_k| = sqrt(Re^2 + Im^2)
   !----------------------------------------------------------------------------
   SUBROUTINE amplitude_espectral(X_real, X_imag, n, amplitude)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: X_real(n), X_imag(n)
      REAL(dp), INTENT(OUT) :: amplitude(n)
      amplitude = SQRT(X_real**2 + X_imag**2)
   END SUBROUTINE amplitude_espectral

   !----------------------------------------------------------------------------
   ! Fase espectral: phi_k = arctan(Im / Re)
   !----------------------------------------------------------------------------
   SUBROUTINE fase_espectral(X_real, X_imag, n, fase)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: X_real(n), X_imag(n)
      REAL(dp), INTENT(OUT) :: fase(n)
      INTEGER :: k
      DO k = 1, n
         fase(k) = ATAN2(X_imag(k), X_real(k) + 1.0D-15)
      END DO
   END SUBROUTINE fase_espectral

   !----------------------------------------------------------------------------
   ! Frequência de amostragem → períodos em dias úteis
   ! freq_k = k / N  (ciclos por amostra)
   ! periodo_k = 1/freq_k = N/k (amostras por ciclo)
   !----------------------------------------------------------------------------
   SUBROUTINE calcular_periodos(n, freq, periodos)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(OUT) :: freq(n), periodos(n)
      INTEGER :: k
      DO k = 1, n
         freq(k) = REAL(k-1, dp) / REAL(n, dp)
         IF (k > 1) THEN
            periodos(k) = REAL(n, dp) / REAL(k-1, dp)
         ELSE
            periodos(k) = HUGE(1.0_dp)
         END IF
      END DO
   END SUBROUTINE calcular_periodos

   !----------------------------------------------------------------------------
   ! Identifica os K ciclos dominantes (maior amplitude espectral)
   !----------------------------------------------------------------------------
   SUBROUTINE ciclos_dominantes(amplitude, periodos, n, K, amp_dom, per_dom)
      INTEGER,  INTENT(IN)  :: n, K
      REAL(dp), INTENT(IN)  :: amplitude(n), periodos(n)
      REAL(dp), INTENT(OUT) :: amp_dom(K), per_dom(K)
      REAL(dp) :: amp_copia(n), max_amp
      INTEGER :: idx_max, i, j
      amp_copia = amplitude
      amp_copia(1) = 0.0_dp    ! ignora frequência zero (DC)
      DO i = 1, K
         max_amp = -1.0_dp;  idx_max = 1
         DO j = 2, n/2
            IF (amp_copia(j) > max_amp) THEN
               max_amp = amp_copia(j);  idx_max = j
            END IF
         END DO
         amp_dom(i) = max_amp
         per_dom(i) = periodos(idx_max)
         amp_copia(idx_max) = -1.0_dp   ! marca como usado
      END DO
   END SUBROUTINE ciclos_dominantes

   !----------------------------------------------------------------------------
   ! Filtro passa-baixa: zera frequências acima de freq_corte
   !----------------------------------------------------------------------------
   SUBROUTINE filtro_passa_baixa(X_real, X_imag, n, freq_corte)
      INTEGER,  INTENT(IN)    :: n
      REAL(dp), INTENT(INOUT) :: X_real(n), X_imag(n)
      REAL(dp), INTENT(IN)    :: freq_corte
      INTEGER :: k
      DO k = 1, n
         IF (REAL(k-1, dp) / REAL(n, dp) > freq_corte) THEN
            X_real(k) = 0.0_dp;  X_imag(k) = 0.0_dp
         END IF
      END DO
   END SUBROUTINE filtro_passa_baixa

   !----------------------------------------------------------------------------
   ! Reconstrução sintética usando apenas os K ciclos dominantes
   ! x_hat(t) = sum_{k=1}^{K} amp_k * cos(2*pi*t/per_k + fase_k)
   !----------------------------------------------------------------------------
   SUBROUTINE reconstruir_sintetico(amp_dom, per_dom, fase_dom, K, n, x_sint)
      INTEGER,  INTENT(IN)  :: K, n
      REAL(dp), INTENT(IN)  :: amp_dom(K), per_dom(K), fase_dom(K)
      REAL(dp), INTENT(OUT) :: x_sint(n)
      INTEGER :: i_t, i_k
      x_sint = 0.0_dp
      DO i_t = 1, n
         DO i_k = 1, K
            IF (per_dom(i_k) > 1.0D-10) &
               x_sint(i_t) = x_sint(i_t) + amp_dom(i_k) / REAL(n, dp) * &
               COS(2.0_dp * PI * REAL(i_t-1, dp) / per_dom(i_k) + fase_dom(i_k))
         END DO
      END DO
   END SUBROUTINE reconstruir_sintetico

   !----------------------------------------------------------------------------
   ! Janela de Hann (reduz vazamento espectral)
   ! w(n) = 0.5*(1 - cos(2*pi*n/N))
   !----------------------------------------------------------------------------
   SUBROUTINE janela_hann(sinal, n, sinal_janelado)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: sinal(n)
      REAL(dp), INTENT(OUT) :: sinal_janelado(n)
      INTEGER :: i
      DO i = 1, n
         sinal_janelado(i) = sinal(i) * 0.5_dp * (1.0_dp - &
            COS(2.0_dp * PI * REAL(i-1, dp) / REAL(n-1, dp)))
      END DO
   END SUBROUTINE janela_hann

END MODULE mod_fourier
