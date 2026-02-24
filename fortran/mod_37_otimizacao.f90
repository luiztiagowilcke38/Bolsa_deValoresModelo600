!==============================================================================
! Módulo 37: Otimização Numérica para Calibração de Modelos
! Autor: Luiz Tiago Wilcke
! Descrição: Algoritmos de otimização: gradiente descendente, Nelder-Mead
!            (simplex), algoritmo genético simplificado e BFGS quasi-Newton
!            para calibração de modelos financeiros aos dados da B3.
!==============================================================================
MODULE mod_otimizacao
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! Gradiente Descendente com momentum (SGD com momentum)
   ! theta_{t+1} = theta_t - lr*(1-beta)*grad - beta*delta_{t-1}
   !----------------------------------------------------------------------------
   SUBROUTINE gradiente_descendente(theta, grad_func, n_params, &
      lr, beta_momentum, n_iter, tol)
      INTEGER,  INTENT(IN)    :: n_params, n_iter
      REAL(dp), INTENT(INOUT) :: theta(n_params)
      REAL(dp), INTENT(IN)    :: lr, beta_momentum, tol
      INTERFACE
         SUBROUTINE grad_func(theta, n, grad)
            import :: dp
            INTEGER, INTENT(IN) :: n
            REAL(dp), INTENT(IN)  :: theta(n)
            REAL(dp), INTENT(OUT) :: grad(n)
         END SUBROUTINE grad_func
      END INTERFACE
      REAL(dp) :: grad(n_params), delta(n_params)
      INTEGER :: iter
      delta = 0.0_dp
      DO iter = 1, n_iter
         CALL grad_func(theta, n_params, grad)
         delta = beta_momentum * delta + (1.0_dp - beta_momentum) * grad
         theta = theta - lr * delta
         IF (SQRT(SUM(grad**2)) < tol) EXIT
      END DO
   END SUBROUTINE gradiente_descendente

   !----------------------------------------------------------------------------
   ! Nelder-Mead Simplex para otimização sem gradiente
   ! Minimiza f(x) para x em R^n
   !----------------------------------------------------------------------------
   SUBROUTINE nelder_mead(f, n, x0, x_otimo, f_otimo, n_iter, tol)
      INTEGER,  INTENT(IN)  :: n, n_iter
      REAL(dp), INTENT(IN)  :: x0(n), tol
      REAL(dp), INTENT(OUT) :: x_otimo(n), f_otimo
      INTERFACE
         REAL(dp) FUNCTION f(x, n)
            import :: dp
            INTEGER, INTENT(IN) :: n
            REAL(dp), INTENT(IN) :: x(n)
         END FUNCTION f
      END INTERFACE
      REAL(dp) :: simplex(n+1, n), f_vals(n+1), centroide(n)
      REAL(dp) :: x_r(n), x_e(n), x_c(n), f_r, f_e, f_c
      REAL(dp), PARAMETER :: alpha_nm = 1.0_dp, gamma_nm = 2.0_dp
      REAL(dp), PARAMETER :: rho_nm   = 0.5_dp, sigma_nm = 0.5_dp
      INTEGER :: i, j, iter, idx_max, idx_min
      REAL(dp) :: perturbacao
      ! Inicializa simplex ao redor de x0
      DO i = 1, n+1
         simplex(i,:) = x0
         IF (i > 1) THEN
            perturbacao = MERGE(0.05_dp*x0(i-1), 0.00025_dp, ABS(x0(i-1)) > 1.0D-10)
            simplex(i, i-1) = simplex(i, i-1) + perturbacao
         END IF
         f_vals(i) = f(simplex(i,:), n)
      END DO
      DO iter = 1, n_iter
         ! Ordena: f_vals(1) = mínimo, f_vals(n+1) = máximo
         CALL ordenar_simplex(f_vals, simplex, n+1, n)
         idx_min = 1;  idx_max = n+1
         IF (MAXVAL(ABS(f_vals - f_vals(1))) < tol) EXIT
         ! Centroide (sem o pior ponto)
         centroide = SUM(simplex(1:n,:), DIM=1) / REAL(n, dp)
         ! Reflexão: x_r = centroide + alpha*(centroide - x_max)
         x_r = centroide + alpha_nm*(centroide - simplex(idx_max,:))
         f_r = f(x_r, n)
         IF (f_r < f_vals(idx_min)) THEN
            ! Expansão
            x_e = centroide + gamma_nm*(x_r - centroide)
            f_e = f(x_e, n)
            IF (f_e < f_r) THEN
               simplex(idx_max,:) = x_e;  f_vals(idx_max) = f_e
            ELSE
               simplex(idx_max,:) = x_r;  f_vals(idx_max) = f_r
            END IF
         ELSE IF (f_r < f_vals(n)) THEN
            simplex(idx_max,:) = x_r;  f_vals(idx_max) = f_r
         ELSE
            ! Contração
            x_c = centroide + rho_nm*(simplex(idx_max,:) - centroide)
            f_c = f(x_c, n)
            IF (f_c < f_vals(idx_max)) THEN
               simplex(idx_max,:) = x_c;  f_vals(idx_max) = f_c
            ELSE
               ! Redução
               DO i = 2, n+1
                  simplex(i,:) = simplex(1,:) + sigma_nm*(simplex(i,:) - simplex(1,:))
                  f_vals(i)    = f(simplex(i,:), n)
               END DO
            END IF
         END IF
      END DO
      x_otimo = simplex(1,:);  f_otimo = f_vals(1)
   END SUBROUTINE nelder_mead

   !----------------------------------------------------------------------------
   ! Ordena simplex por ordem crescente de f_vals (bubble sort)
   !----------------------------------------------------------------------------
   SUBROUTINE ordenar_simplex(f_vals, simplex, m, n)
      INTEGER,  INTENT(IN)    :: m, n
      REAL(dp), INTENT(INOUT) :: f_vals(m), simplex(m, n)
      REAL(dp) :: temp_f, temp_x(n)
      INTEGER  :: i, j
      DO i = 1, m-1
         DO j = 1, m-i
            IF (f_vals(j) > f_vals(j+1)) THEN
               temp_f    = f_vals(j);    f_vals(j)    = f_vals(j+1); f_vals(j+1) = temp_f
               temp_x    = simplex(j,:); simplex(j,:) = simplex(j+1,:); simplex(j+1,:) = temp_x
            END IF
         END DO
      END DO
   END SUBROUTINE ordenar_simplex

   !----------------------------------------------------------------------------
   ! Algoritmo Genético simplificado (GA) — para otimização global
   !----------------------------------------------------------------------------
   SUBROUTINE algoritmo_genetico(f, n_params, x_min, x_max, n_pop, n_ger, &
      taxa_mutacao, x_melhor, f_melhor)
      INTEGER,  INTENT(IN)  :: n_params, n_pop, n_ger
      REAL(dp), INTENT(IN)  :: x_min(n_params), x_max(n_params), taxa_mutacao
      REAL(dp), INTENT(OUT) :: x_melhor(n_params), f_melhor
      INTERFACE
         REAL(dp) FUNCTION f(x, n)
            import :: dp
            INTEGER, INTENT(IN) :: n
            REAL(dp), INTENT(IN) :: x(n)
         END FUNCTION f
      END INTERFACE
      REAL(dp) :: populacao(n_pop, n_params), f_pop(n_pop), nova_pop(n_pop, n_params)
      REAL(dp) :: u, filho(n_params)
      INTEGER :: ger, i, j, idx1, idx2
      ! População inicial uniforme
      DO i = 1, n_pop
         DO j = 1, n_params
            CALL RANDOM_NUMBER(u)
            populacao(i,j) = x_min(j) + u * (x_max(j) - x_min(j))
         END DO
         f_pop(i) = f(populacao(i,:), n_params)
      END DO
      x_melhor = populacao(1,:);  f_melhor = f_pop(1)
      DO ger = 1, n_ger
         DO i = 1, n_pop
            CALL RANDOM_NUMBER(u);  idx1 = INT(u * n_pop) + 1
            CALL RANDOM_NUMBER(u);  idx2 = INT(u * n_pop) + 1
            idx1 = MAX(1,MIN(idx1,n_pop));  idx2 = MAX(1,MIN(idx2,n_pop))
            ! Cruzamento (crossover uniforme)
            DO j = 1, n_params
               CALL RANDOM_NUMBER(u)
               IF (u < 0.5_dp) THEN
                  filho(j) = populacao(idx1,j)
               ELSE
                  filho(j) = populacao(idx2,j)
               END IF
               ! Mutação gaussiana
               CALL RANDOM_NUMBER(u)
               IF (u < taxa_mutacao) THEN
                  CALL RANDOM_NUMBER(u)
                  filho(j) = filho(j) + (u-0.5_dp)*2.0_dp*(x_max(j)-x_min(j))*0.1_dp
                  filho(j) = MAX(x_min(j), MIN(x_max(j), filho(j)))
               END IF
            END DO
            nova_pop(i,:) = filho
            f_pop(i)      = f(filho, n_params)
            IF (f_pop(i) < f_melhor) THEN
               f_melhor = f_pop(i);  x_melhor = filho
            END IF
         END DO
         populacao = nova_pop
      END DO
   END SUBROUTINE algoritmo_genetico

END MODULE mod_otimizacao
