!==============================================================================
! Módulo 02: Equações Diferenciais Financeiras
! Autor: Luiz Tiago Wilcke
! Descrição: Implementa métodos numéricos para equações diferenciais
!            ordinárias (EDO) e parciais (EDP) aplicadas a modelos
!            financeiros. Inclui Runge-Kutta 4ª ordem, Euler, Crank-Nicolson
!            e métodos implícitos para equação de Black-Scholes.
!==============================================================================
MODULE mod_equacoes_diferenciais
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)
   REAL(dp), PARAMETER :: PI = 3.14159265358979323846_dp

CONTAINS

   !----------------------------------------------------------------------------
   ! Método de Euler explícito para EDO: dy/dt = f(t, y)
   ! y_{n+1} = y_n + h * f(t_n, y_n)
   !----------------------------------------------------------------------------
   SUBROUTINE euler_explicito(f, y0, t0, tf, n_passos, t_sol, y_sol)
      INTERFACE
         REAL(dp) FUNCTION f(t, y)
            import :: dp
            REAL(dp), INTENT(IN) :: t, y
         END FUNCTION f
      END INTERFACE
      REAL(dp), INTENT(IN)  :: y0, t0, tf
      INTEGER,  INTENT(IN)  :: n_passos
      REAL(dp), INTENT(OUT) :: t_sol(:), y_sol(:)
      REAL(dp) :: h, t_atual, y_atual
      INTEGER :: i
      h = (tf - t0) / REAL(n_passos, dp)
      t_atual = t0;  y_atual = y0
      t_sol(1) = t0; y_sol(1) = y0
      DO i = 1, n_passos
         y_atual = y_atual + h * f(t_atual, y_atual)
         t_atual = t_atual + h
         IF (i+1 <= SIZE(t_sol)) THEN
            t_sol(i+1) = t_atual
            y_sol(i+1) = y_atual
         END IF
      END DO
   END SUBROUTINE euler_explicito

   !----------------------------------------------------------------------------
   ! Runge-Kutta de 4ª ordem (RK4) para EDO escalar
   ! k1 = h*f(t_n, y_n)
   ! k2 = h*f(t_n + h/2, y_n + k1/2)
   ! k3 = h*f(t_n + h/2, y_n + k2/2)
   ! k4 = h*f(t_n + h, y_n + k3)
   ! y_{n+1} = y_n + (k1 + 2k2 + 2k3 + k4) / 6
   !----------------------------------------------------------------------------
   SUBROUTINE runge_kutta4(f, y0, t0, tf, n_passos, t_sol, y_sol)
      INTERFACE
         REAL(dp) FUNCTION f(t, y)
            import :: dp
            REAL(dp), INTENT(IN) :: t, y
         END FUNCTION f
      END INTERFACE
      REAL(dp), INTENT(IN)  :: y0, t0, tf
      INTEGER,  INTENT(IN)  :: n_passos
      REAL(dp), INTENT(OUT) :: t_sol(:), y_sol(:)
      REAL(dp) :: h, t_n, y_n, k1, k2, k3, k4
      INTEGER :: i
      h = (tf - t0) / REAL(n_passos, dp)
      t_n = t0;  y_n = y0
      t_sol(1) = t0; y_sol(1) = y0
      DO i = 1, n_passos
         k1 = h * f(t_n,         y_n)
         k2 = h * f(t_n + h/2.0_dp, y_n + k1/2.0_dp)
         k3 = h * f(t_n + h/2.0_dp, y_n + k2/2.0_dp)
         k4 = h * f(t_n + h,     y_n + k3)
         y_n = y_n + (k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4) / 6.0_dp
         t_n = t_n + h
         IF (i+1 <= SIZE(t_sol)) THEN
            t_sol(i+1) = t_n
            y_sol(i+1) = y_n
         END IF
      END DO
   END SUBROUTINE runge_kutta4

   !----------------------------------------------------------------------------
   ! Runge-Kutta 4ª ordem para sistemas de EDOs (vetor)
   ! dY/dt = F(t, Y), Y vetor de dimensão ndim
   !----------------------------------------------------------------------------
   SUBROUTINE rk4_sistema(F_sistema, Y0, t0, tf, n_passos, ndim, Y_sol)
      INTEGER, INTENT(IN) :: ndim, n_passos
      REAL(dp), INTENT(IN) :: Y0(ndim), t0, tf
      REAL(dp), INTENT(OUT) :: Y_sol(ndim, n_passos+1)
      INTERFACE
         SUBROUTINE F_sistema(t, Y, dYdt, ndim)
            import :: dp
            INTEGER, INTENT(IN) :: ndim
            REAL(dp), INTENT(IN) :: t, Y(ndim)
            REAL(dp), INTENT(OUT) :: dYdt(ndim)
         END SUBROUTINE F_sistema
      END INTERFACE
      REAL(dp) :: h, t_n
      REAL(dp) :: K1(ndim), K2(ndim), K3(ndim), K4(ndim), Y_n(ndim), aux(ndim)
      INTEGER :: i
      h = (tf - t0) / REAL(n_passos, dp)
      t_n = t0
      Y_n = Y0
      Y_sol(:, 1) = Y0
      DO i = 1, n_passos
         CALL F_sistema(t_n,           Y_n,             K1, ndim);  K1 = h * K1
         CALL F_sistema(t_n + h/2.0_dp, Y_n + K1/2.0_dp, K2, ndim);  K2 = h * K2
         CALL F_sistema(t_n + h/2.0_dp, Y_n + K2/2.0_dp, K3, ndim);  K3 = h * K3
         CALL F_sistema(t_n + h,        Y_n + K3,         K4, ndim);  K4 = h * K4
         Y_n = Y_n + (K1 + 2.0_dp*K2 + 2.0_dp*K3 + K4) / 6.0_dp
         t_n = t_n + h
         Y_sol(:, i+1) = Y_n
      END DO
   END SUBROUTINE rk4_sistema

   !----------------------------------------------------------------------------
   ! Esquema de Crank-Nicolson para EDP de Black-Scholes:
   ! dV/dt + (1/2)*sigma^2*S^2*d2V/dS2 + r*S*dV/dS - r*V = 0
   ! Discretização implícita (theta=0.5) - incondicionalmente estável
   !----------------------------------------------------------------------------
   SUBROUTINE crank_nicolson_black_scholes(sigma, taxa_juro, S_max, &
      K, T_exp, Ns, Nt, tipo_opcao, &
      grade_S, grade_V)
      REAL(dp), INTENT(IN) :: sigma, taxa_juro, S_max, K, T_exp
      INTEGER,  INTENT(IN) :: Ns, Nt
      CHARACTER(LEN=4), INTENT(IN) :: tipo_opcao    ! "CALL" ou "PUT"
      REAL(dp), INTENT(OUT) :: grade_S(Ns+1), grade_V(Ns+1)
      REAL(dp) :: dS, dt, theta
      REAL(dp) :: a_coef(Ns-1), b_coef(Ns-1), c_coef(Ns-1)
      REAL(dp) :: V_atual(Ns+1), V_novo(Ns-1)
      REAL(dp) :: rhs(Ns-1), alpha, beta_coef, gamma_coef
      INTEGER :: i, j
      dS = S_max / REAL(Ns, dp)
      dt = T_exp / REAL(Nt, dp)
      theta = 0.5_dp   ! Crank-Nicolson
      ! Grade de preços do ativo
      DO i = 0, Ns
         grade_S(i+1) = REAL(i, dp) * dS
      END DO
      ! Condição inicial (expiração)
      DO i = 0, Ns
         IF (tipo_opcao == "CALL") THEN
            V_atual(i+1) = MAX(grade_S(i+1) - K, 0.0_dp)
         ELSE
            V_atual(i+1) = MAX(K - grade_S(i+1), 0.0_dp)
         END IF
      END DO
      ! Iteração temporal (backwards)
      DO j = 1, Nt
         DO i = 2, Ns   ! índice interno i (S = (i-1)*dS)
            alpha      = 0.5_dp * dt * (sigma**2 * REAL(i-1,dp)**2 - taxa_juro * REAL(i-1,dp))
            beta_coef  = -dt * (sigma**2 * REAL(i-1,dp)**2 + taxa_juro)
            gamma_coef = 0.5_dp * dt * (sigma**2 * REAL(i-1,dp)**2 + taxa_juro * REAL(i-1,dp))
            a_coef(i-1) = -theta * alpha
            b_coef(i-1) =  1.0_dp - theta * beta_coef
            c_coef(i-1) = -theta * gamma_coef
            rhs(i-1) = (1.0_dp - theta) * alpha * V_atual(i-1) + &
               (1.0_dp + (1.0_dp - theta) * (-beta_coef)) * V_atual(i) + &
               (1.0_dp - theta) * gamma_coef * V_atual(i+1)
         END DO
         ! Solver tridiagonal (algoritmo de Thomas)
         CALL solver_tridiagonal(a_coef, b_coef, c_coef, rhs, V_novo, Ns-1)
         V_atual(1)    = 0.0_dp
         V_atual(2:Ns) = V_novo
         IF (tipo_opcao == "CALL") THEN
            V_atual(Ns+1) = S_max - K * EXP(-taxa_juro * REAL(Nt-j, dp) * dt)
         ELSE
            V_atual(Ns+1) = 0.0_dp
         END IF
      END DO
      grade_V = V_atual
   END SUBROUTINE crank_nicolson_black_scholes

   !----------------------------------------------------------------------------
   ! Solver tridiagonal pelo algoritmo de Thomas (LU decomposition)
   ! Resolve sistema Ax = d onde A é tridiagonal com diag a(i), b(i), c(i)
   !----------------------------------------------------------------------------
   SUBROUTINE solver_tridiagonal(a, b, c, d, x, n)
      INTEGER, INTENT(IN) :: n
      REAL(dp), INTENT(IN)  :: a(n), b(n), c(n), d(n)
      REAL(dp), INTENT(OUT) :: x(n)
      REAL(dp) :: c_prime(n), d_prime(n), m
      INTEGER :: i
      c_prime(1) = c(1) / b(1)
      d_prime(1) = d(1) / b(1)
      DO i = 2, n
         m = b(i) - a(i) * c_prime(i-1)
         IF (ABS(m) < 1.0D-15) m = 1.0D-15
         c_prime(i) = c(i) / m
         d_prime(i) = (d(i) - a(i) * d_prime(i-1)) / m
      END DO
      x(n) = d_prime(n)
      DO i = n-1, 1, -1
         x(i) = d_prime(i) - c_prime(i) * x(i+1)
      END DO
   END SUBROUTINE solver_tridiagonal

   !----------------------------------------------------------------------------
   ! Equação de difusão estocástica (processo de Itô):
   ! dX = mu*X*dt + sigma*X*dW
   ! Discretização de Euler-Maruyama
   !----------------------------------------------------------------------------
   SUBROUTINE euler_maruyama(X0, mu, sigma, dt, n_passos, n_trajetorias, X_traj)
      REAL(dp), INTENT(IN)  :: X0, mu, sigma, dt
      INTEGER,  INTENT(IN)  :: n_passos, n_trajetorias
      REAL(dp), INTENT(OUT) :: X_traj(n_trajetorias, n_passos+1)
      REAL(dp) :: X_atual, dW, ruido
      INTEGER :: i, j
      DO j = 1, n_trajetorias
         X_traj(j, 1) = X0
         X_atual = X0
         DO i = 1, n_passos
            CALL RANDOM_NUMBER(ruido)
            ! Box-Muller para gerar N(0,1)
            dW = SQRT(-2.0_dp * LOG(ruido + 1.0D-15)) * &
               COS(2.0_dp * PI * ruido) * SQRT(dt)
            X_atual = X_atual + mu * X_atual * dt + sigma * X_atual * dW
            X_traj(j, i+1) = X_atual
         END DO
      END DO
   END SUBROUTINE euler_maruyama

END MODULE mod_equacoes_diferenciais
