!==============================================================================
! Módulo 10: Correlação e Dependência Estatística
! Autor: Luiz Tiago Wilcke
! Descrição: Coeficientes de correlação (Pearson, Spearman, Kendall),
!            correlação dinâmica (DCC), correlação rolante e testes de
!            significância para construção de portfólios na B3.
!==============================================================================
MODULE mod_correlacao
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! Correlação de Pearson: rho = Cov(X,Y) / (sigma_X * sigma_Y)
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION pearson(x, y) RESULT(rho)
      REAL(dp), INTENT(IN) :: x(:), y(:)
      REAL(dp) :: mx, my, num, den_x, den_y
      INTEGER :: n, i
      n = MIN(SIZE(x), SIZE(y))
      mx = SUM(x(1:n)) / REAL(n, dp)
      my = SUM(y(1:n)) / REAL(n, dp)
      num = 0.0_dp; den_x = 0.0_dp; den_y = 0.0_dp
      DO i = 1, n
         num   = num   + (x(i) - mx) * (y(i) - my)
         den_x = den_x + (x(i) - mx)**2
         den_y = den_y + (y(i) - my)**2
      END DO
      IF (den_x > 1.0D-15 .AND. den_y > 1.0D-15) THEN
         rho = num / SQRT(den_x * den_y)
      ELSE
         rho = 0.0_dp
      END IF
   END FUNCTION pearson

   !----------------------------------------------------------------------------
   ! Rank de um vetor (utilizado para Spearman e Kendall)
   !----------------------------------------------------------------------------
   SUBROUTINE calcular_rank(x, ranks)
      REAL(dp), INTENT(IN)  :: x(:)
      REAL(dp), INTENT(OUT) :: ranks(:)
      INTEGER :: n, i, j
      REAL(dp) :: pos
      n = SIZE(x)
      DO i = 1, n
         pos = 1.0_dp
         DO j = 1, n
            IF (x(j) < x(i)) pos = pos + 1.0_dp
         END DO
         ranks(i) = pos   ! simplificado sem tratamento de empates
      END DO
   END SUBROUTINE calcular_rank

   !----------------------------------------------------------------------------
   ! Correlação de Spearman (baseada em ranks)
   ! rho_S = 1 - 6*sum(di^2) / (n*(n^2-1))
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION spearman(x, y) RESULT(rho_s)
      REAL(dp), INTENT(IN) :: x(:), y(:)
      REAL(dp), ALLOCATABLE :: rx(:), ry(:)
      REAL(dp) :: soma_d2
      INTEGER :: n, i
      n = MIN(SIZE(x), SIZE(y))
      ALLOCATE(rx(n), ry(n))
      CALL calcular_rank(x(1:n), rx)
      CALL calcular_rank(y(1:n), ry)
      soma_d2 = 0.0_dp
      DO i = 1, n
         soma_d2 = soma_d2 + (rx(i) - ry(i))**2
      END DO
      rho_s = 1.0_dp - 6.0_dp * soma_d2 / (REAL(n, dp) * (REAL(n, dp)**2 - 1.0_dp))
      DEALLOCATE(rx, ry)
   END FUNCTION spearman

   !----------------------------------------------------------------------------
   ! Correlação de Kendall-tau
   ! tau = (C - D) / (n*(n-1)/2)
   ! C = pares concordantes, D = pares discordantes
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION kendall_tau(x, y) RESULT(tau)
      REAL(dp), INTENT(IN) :: x(:), y(:)
      INTEGER :: n, i, j, concordantes, discordantes
      n = MIN(SIZE(x), SIZE(y))
      concordantes = 0; discordantes = 0
      DO i = 1, n - 1
         DO j = i+1, n
            IF ((x(j) - x(i)) * (y(j) - y(i)) > 0.0_dp) THEN
               concordantes = concordantes + 1
            ELSE IF ((x(j) - x(i)) * (y(j) - y(i)) < 0.0_dp) THEN
               discordantes = discordantes + 1
            END IF
         END DO
      END DO
      tau = REAL(concordantes - discordantes, dp) / &
         (REAL(n, dp) * REAL(n - 1, dp) / 2.0_dp)
   END FUNCTION kendall_tau

   !----------------------------------------------------------------------------
   ! Matriz de correlação de Pearson para N ativos
   !----------------------------------------------------------------------------
   SUBROUTINE matriz_correlacao(retornos, n_ativos, n_obs, mat_corr)
      INTEGER,  INTENT(IN)  :: n_ativos, n_obs
      REAL(dp), INTENT(IN)  :: retornos(n_obs, n_ativos)
      REAL(dp), INTENT(OUT) :: mat_corr(n_ativos, n_ativos)
      INTEGER :: i, j
      DO i = 1, n_ativos
         DO j = i, n_ativos
            mat_corr(i, j) = pearson(retornos(:, i), retornos(:, j))
            mat_corr(j, i) = mat_corr(i, j)
         END DO
      END DO
   END SUBROUTINE matriz_correlacao

   !----------------------------------------------------------------------------
   ! Correlação rolante (janela deslizante)
   !----------------------------------------------------------------------------
   SUBROUTINE correlacao_rolante(x, y, janela, corr_rol)
      REAL(dp), INTENT(IN)  :: x(:), y(:)
      INTEGER,  INTENT(IN)  :: janela
      REAL(dp), INTENT(OUT) :: corr_rol(:)
      INTEGER :: n, i
      n = MIN(SIZE(x), SIZE(y))
      corr_rol = 0.0_dp
      DO i = janela, n
         corr_rol(i) = pearson(x(i-janela+1:i), y(i-janela+1:i))
      END DO
   END SUBROUTINE correlacao_rolante

   !----------------------------------------------------------------------------
   ! Teste de significância da correlação de Pearson (aproximação t)
   ! t = rho * sqrt(n-2) / sqrt(1-rho^2)
   ! gl = n - 2
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION teste_t_correlacao(rho, n) RESULT(t_stat)
      REAL(dp), INTENT(IN) :: rho
      INTEGER,  INTENT(IN) :: n
      REAL(dp) :: denominador
      denominador = SQRT(MAX(1.0_dp - rho**2, 1.0D-15))
      t_stat = rho * SQRT(REAL(n - 2, dp)) / denominador
   END FUNCTION teste_t_correlacao

   !----------------------------------------------------------------------------
   ! Distância de correlação (métrica para clustering)
   ! d(i,j) = sqrt(2 * (1 - rho(i,j)))
   !----------------------------------------------------------------------------
   SUBROUTINE distancia_correlacao(mat_corr, n, mat_dist)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: mat_corr(n, n)
      REAL(dp), INTENT(OUT) :: mat_dist(n, n)
      INTEGER :: i, j
      DO i = 1, n
         DO j = 1, n
            mat_dist(i,j) = SQRT(2.0_dp * MAX(1.0_dp - mat_corr(i,j), 0.0_dp))
         END DO
      END DO
   END SUBROUTINE distancia_correlacao

   !----------------------------------------------------------------------------
   ! Autocorrelação de lag k para uma série
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION autocorrelacao(serie, lag) RESULT(acf_k)
      REAL(dp), INTENT(IN) :: serie(:)
      INTEGER,  INTENT(IN) :: lag
      REAL(dp) :: media, var, cov
      INTEGER :: n, i
      n = SIZE(serie)
      media = SUM(serie) / REAL(n, dp)
      var   = SUM((serie - media)**2) / REAL(n, dp)
      cov   = 0.0_dp
      DO i = 1, n - lag
         cov = cov + (serie(i) - media) * (serie(i+lag) - media)
      END DO
      cov = cov / REAL(n, dp)
      IF (ABS(var) > 1.0D-15) THEN
         acf_k = cov / var
      ELSE
         acf_k = 0.0_dp
      END IF
   END FUNCTION autocorrelacao

END MODULE mod_correlacao
