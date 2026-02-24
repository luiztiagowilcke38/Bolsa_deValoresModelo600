!==============================================================================
! Módulo 26: Equações de Difusão Financeira
! Autor: Luiz Tiago Wilcke
! Descrição: Solução numérica de equações de difusão-advecção para preços
!            de ativos: Fokker-Planck financeira, equação do calor (Black-
!            Scholes transformada), método ADI (direções alternadas),
!            e difusão anisotrópica para volatilidade local de superfícies.
!==============================================================================
MODULE mod_difusao
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)
   REAL(dp), PARAMETER :: PI = 3.14159265358979323846_dp

CONTAINS

   !----------------------------------------------------------------------------
   ! Equação de Fokker-Planck (evolução da densidade de probabilidade)
   ! dp/dt = -d/dx[mu(x)*p] + 0.5*d²/dx²[sigma²(x)*p]
   ! Solução numérica por diferenças finitas upwind + difusão centrada
   !----------------------------------------------------------------------------
   SUBROUTINE fokker_planck_1d(p0, x_grade, n_x, mu_drift, sigma_dif, &
      dt, n_passos, p_final)
      INTEGER,  INTENT(IN)  :: n_x, n_passos
      REAL(dp), INTENT(IN)  :: p0(n_x), x_grade(n_x), mu_drift, sigma_dif, dt
      REAL(dp), INTENT(OUT) :: p_final(n_x)
      REAL(dp) :: p(n_x), p_novo(n_x), dx, difusao, adveccao
      INTEGER :: i, t
      dx = (x_grade(n_x) - x_grade(1)) / REAL(n_x - 1, dp)
      p = p0
      DO t = 1, n_passos
         p_novo = p
         DO i = 2, n_x - 1
            ! Difusão (centrada): 0.5 * sigma² * d²p/dx²
            difusao = 0.5_dp * sigma_dif**2 * (p(i+1) - 2.0_dp*p(i) + p(i-1)) / dx**2
            ! Advecção upwind: mu * dp/dx
            IF (mu_drift >= 0.0_dp) THEN
               adveccao = mu_drift * (p(i) - p(i-1)) / dx
            ELSE
               adveccao = mu_drift * (p(i+1) - p(i)) / dx
            END IF
            p_novo(i) = p(i) + dt * (-adveccao + difusao)
         END DO
         ! Condições de contorno (reflexão)
         p_novo(1) = p_novo(2)
         p_novo(n_x) = p_novo(n_x - 1)
         ! Normaliza para manter densidade de probabilidade
         IF (SUM(p_novo) * dx > 1.0D-15) p_novo = p_novo / (SUM(p_novo) * dx)
         p = p_novo
      END DO
      p_final = p
   END SUBROUTINE fokker_planck_1d

   !----------------------------------------------------------------------------
   ! Black-Scholes via transformação em equação do calor
   ! Substituição: S = K*e^x, t = T-tau/(sigma^2/2), u = V/(K*e^x)
   ! Resulta em: du/dtau = d²u/dx²  (equação do calor padrão)
   !----------------------------------------------------------------------------
   SUBROUTINE equacao_calor_bs(u0, x_min, x_max, n_x, tau_T, n_t, u_final)
      INTEGER,  INTENT(IN)  :: n_x, n_t
      REAL(dp), INTENT(IN)  :: u0(n_x), x_min, x_max, tau_T
      REAL(dp), INTENT(OUT) :: u_final(n_x)
      REAL(dp) :: u(n_x), u_novo(n_x), dx, dt_calor, r_CFL
      INTEGER :: i, t
      dx = (x_max - x_min) / REAL(n_x - 1, dp)
      dt_calor = tau_T / REAL(n_t, dp)
      r_CFL = dt_calor / dx**2    ! deve ser <= 0.5 para estabilidade
      u = u0
      DO t = 1, n_t
         DO i = 2, n_x - 1
            u_novo(i) = u(i) + r_CFL * (u(i+1) - 2.0_dp*u(i) + u(i-1))
         END DO
         u_novo(1)   = 0.0_dp    ! condição de contorno para call: u=0 para S→0
         u_novo(n_x) = u_novo(n_x - 1)
         u = u_novo
      END DO
      u_final = u
   END SUBROUTINE equacao_calor_bs

   !----------------------------------------------------------------------------
   ! Método ADI (Alternating Direction Implicit) para EDP 2D (ex: Heston)
   ! passo X: implícito em x, explícito em y
   ! passo Y: implícito em y, explícito em x
   !----------------------------------------------------------------------------
   SUBROUTINE adi_2d(V, n_x, n_y, a_x, b_x, c_x, a_y, b_y, c_y, V_novo)
      INTEGER,  INTENT(IN)  :: n_x, n_y
      REAL(dp), INTENT(IN)  :: V(n_x, n_y)
      REAL(dp), INTENT(IN)  :: a_x(n_x), b_x(n_x), c_x(n_x)
      REAL(dp), INTENT(IN)  :: a_y(n_y), b_y(n_y), c_y(n_y)
      REAL(dp), INTENT(OUT) :: V_novo(n_x, n_y)
      REAL(dp) :: V_half(n_x, n_y), rhs_x(n_x), rhs_y(n_y), sol(MAX(n_x,n_y))
      INTEGER :: i, j
      ! Passo 1: sweep em x (resolve sistema tridiagonal para cada j fixo)
      DO j = 2, n_y - 1
         rhs_x = V(:,j)
         rhs_x(1)   = 0.0_dp
         rhs_x(n_x) = 0.0_dp
         CALL thomas_tridiag(a_x, b_x, c_x, rhs_x, sol(1:n_x), n_x)
         V_half(:,j) = sol(1:n_x)
      END DO
      V_half(:,1) = 0.0_dp;  V_half(:,n_y) = V(:,n_y)
      ! Passo 2: sweep em y
      DO i = 2, n_x - 1
         rhs_y = V_half(i,:)
         rhs_y(1)   = 0.0_dp
         rhs_y(n_y) = V_half(i, n_y)
         CALL thomas_tridiag(a_y, b_y, c_y, rhs_y, sol(1:n_y), n_y)
         V_novo(i,:) = sol(1:n_y)
      END DO
      V_novo(1,:) = 0.0_dp;  V_novo(n_x,:) = V_half(n_x,:)
   END SUBROUTINE adi_2d

   !----------------------------------------------------------------------------
   ! Algoritmo de Thomas (tridiagonal) — solver interno para ADI
   !----------------------------------------------------------------------------
   SUBROUTINE thomas_tridiag(a, b, c, d, x, n)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: a(n), b(n), c(n), d(n)
      REAL(dp), INTENT(OUT) :: x(n)
      REAL(dp) :: c2(n), d2(n), m
      INTEGER :: i
      c2(1) = c(1) / (b(1) + 1.0D-15)
      d2(1) = d(1) / (b(1) + 1.0D-15)
      DO i = 2, n
         m     = b(i) - a(i) * c2(i-1)
         IF (ABS(m) < 1.0D-15) m = 1.0D-15
         c2(i) = c(i) / m
         d2(i) = (d(i) - a(i) * d2(i-1)) / m
      END DO
      x(n) = d2(n)
      DO i = n-1, 1, -1
         x(i) = d2(i) - c2(i) * x(i+1)
      END DO
   END SUBROUTINE thomas_tridiag

   !----------------------------------------------------------------------------
   ! Densidade log-normal (distribuição dos preços no modelo BS)
   ! f(S) = (1/(S*sigma*sqrt(2*pi*T))) * exp(-(ln(S/S0)-(mu-sigma^2/2)*T)^2
   !                                           / (2*sigma^2*T))
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION densidade_lognormal(S, S0, mu, sigma, T) RESULT(f)
      REAL(dp), INTENT(IN) :: S, S0, mu, sigma, T
      REAL(dp) :: media_log, var_log
      IF (S <= 0.0_dp .OR. T <= 0.0_dp) THEN
         f = 0.0_dp;  RETURN
      END IF
      media_log = LOG(S0) + (mu - 0.5_dp * sigma**2) * T
      var_log   = sigma**2 * T
      f = EXP(-(LOG(S) - media_log)**2 / (2.0_dp * var_log)) / &
         (S * SQRT(2.0_dp * PI * var_log))
   END FUNCTION densidade_lognormal

END MODULE mod_difusao
