!==============================================================================
! Módulo 18: Redes Neurais para Previsão Financeira
! Autor: Luiz Tiago Wilcke
! Descrição: Rede neural feedforward (MLP) com backpropagation para previsão
!            de retornos e classificação de sinais de compra/venda em ações
!            da B3. Inclui funções de ativação, inicialização Xavier e
!            regularização L2 (weight decay).
!==============================================================================
MODULE mod_redes_neurais
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

   TYPE :: CamadaDensa
      REAL(dp), ALLOCATABLE :: pesos(:,:)    ! pesos(n_saida, n_entrada)
      REAL(dp), ALLOCATABLE :: bias(:)       ! bias(n_saida)
      REAL(dp), ALLOCATABLE :: grad_pesos(:,:)
      REAL(dp), ALLOCATABLE :: grad_bias(:)
      INTEGER :: n_entrada, n_saida
   END TYPE CamadaDensa

CONTAINS

   !----------------------------------------------------------------------------
   ! Inicialização Xavier/Glorot
   ! W ~ U[-sqrt(6/(n_in+n_out)), sqrt(6/(n_in+n_out))]
   !----------------------------------------------------------------------------
   SUBROUTINE inicializar_xavier(camada)
      TYPE(CamadaDensa), INTENT(INOUT) :: camada
      REAL(dp) :: limite, u
      INTEGER :: i, j
      limite = SQRT(6.0_dp / REAL(camada%n_entrada + camada%n_saida, dp))
      DO i = 1, camada%n_saida
         DO j = 1, camada%n_entrada
            CALL RANDOM_NUMBER(u)
            camada%pesos(i,j) = u * 2.0_dp * limite - limite
         END DO
         CALL RANDOM_NUMBER(u)
         camada%bias(i) = (u - 0.5_dp) * 0.1_dp
      END DO
      camada%grad_pesos = 0.0_dp
      camada%grad_bias  = 0.0_dp
   END SUBROUTINE inicializar_xavier

   !----------------------------------------------------------------------------
   ! Função de ativação ReLU: f(x) = max(0, x)
   !----------------------------------------------------------------------------
   ELEMENTAL REAL(dp) FUNCTION relu(x) RESULT(y)
      REAL(dp), INTENT(IN) :: x
      y = MAX(0.0_dp, x)
   END FUNCTION relu

   !----------------------------------------------------------------------------
   ! Derivada de ReLU: f'(x) = 1 se x>0, 0 cc
   !----------------------------------------------------------------------------
   ELEMENTAL REAL(dp) FUNCTION drelu(x) RESULT(dy)
      REAL(dp), INTENT(IN) :: x
      dy = MERGE(1.0_dp, 0.0_dp, x > 0.0_dp)
   END FUNCTION drelu

   !----------------------------------------------------------------------------
   ! Função Sigmoid: sigma(x) = 1 / (1 + e^{-x})
   !----------------------------------------------------------------------------
   ELEMENTAL REAL(dp) FUNCTION sigmoid(x) RESULT(y)
      REAL(dp), INTENT(IN) :: x
      y = 1.0_dp / (1.0_dp + EXP(-SIGN(MIN(ABS(x), 500.0_dp), x)))
   END FUNCTION sigmoid

   !----------------------------------------------------------------------------
   ! Tanh com clipping numérico
   !----------------------------------------------------------------------------
   ELEMENTAL REAL(dp) FUNCTION tanh_ativacao(x) RESULT(y)
      REAL(dp), INTENT(IN) :: x
      REAL(dp) :: xc
      xc = SIGN(MIN(ABS(x), 30.0_dp), x)
      y  = TANH(xc)
   END FUNCTION tanh_ativacao

   !----------------------------------------------------------------------------
   ! Forward pass de uma camada densa: y = f(W*x + b)
   !----------------------------------------------------------------------------
   SUBROUTINE forward_camada(camada, x_entrada, a_saida, tipo_ativacao)
      TYPE(CamadaDensa), INTENT(IN)  :: camada
      REAL(dp), INTENT(IN)  :: x_entrada(camada%n_entrada)
      REAL(dp), INTENT(OUT) :: a_saida(camada%n_saida)
      CHARACTER(LEN=6), INTENT(IN) :: tipo_ativacao  ! "relu", "sigm", "tanh", "linear"
      REAL(dp) :: z(camada%n_saida)
      INTEGER :: i
      z = MATMUL(camada%pesos, x_entrada) + camada%bias
      SELECT CASE (TRIM(tipo_ativacao))
       CASE ("relu")
         DO i = 1, camada%n_saida
            a_saida(i) = relu(z(i))
         END DO
       CASE ("sigm")
         DO i = 1, camada%n_saida
            a_saida(i) = sigmoid(z(i))
         END DO
       CASE ("tanh")
         DO i = 1, camada%n_saida
            a_saida(i) = tanh_ativacao(z(i))
         END DO
       CASE DEFAULT
         a_saida = z    ! linear (saída)
      END SELECT
   END SUBROUTINE forward_camada

   !----------------------------------------------------------------------------
   ! Erro quadrático médio (MSE): L = (1/n) * sum(y_hat - y)^2
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION mse(y_pred, y_true, n) RESULT(loss)
      INTEGER,  INTENT(IN) :: n
      REAL(dp), INTENT(IN) :: y_pred(n), y_true(n)
      loss = SUM((y_pred - y_true)**2) / REAL(n, dp)
   END FUNCTION mse

   !----------------------------------------------------------------------------
   ! Backpropagation para rede 2 camadas (entrada → hidden → saída)
   ! Atualiza gradientes acumulados em camadas
   !----------------------------------------------------------------------------
   SUBROUTINE backprop_2camadas(camada1, camada2, x_ent, h_out, y_out, &
      y_true, n_saida, eta, lambda_reg)
      TYPE(CamadaDensa), INTENT(INOUT) :: camada1, camada2
      REAL(dp), INTENT(IN) :: x_ent(:), h_out(:), y_out(:), y_true(:)
      INTEGER,  INTENT(IN) :: n_saida
      REAL(dp), INTENT(IN) :: eta, lambda_reg
      REAL(dp) :: delta_saida(n_saida), delta_hidden(camada1%n_saida)
      INTEGER :: i, j
      ! Gradiente na camada de saída (MSE + linear)
      DO i = 1, n_saida
         delta_saida(i) = 2.0_dp * (y_out(i) - y_true(i)) / REAL(n_saida, dp)
      END DO
      ! Backprop para camada2
      DO i = 1, camada2%n_saida
         DO j = 1, camada2%n_entrada
            camada2%grad_pesos(i,j) = camada2%grad_pesos(i,j) + &
               delta_saida(i) * h_out(j) + &
               lambda_reg * camada2%pesos(i,j)
         END DO
         camada2%grad_bias(i) = camada2%grad_bias(i) + delta_saida(i)
      END DO
      ! Propaga delta para camada oculta (ReLU)
      DO j = 1, camada1%n_saida
         delta_hidden(j) = DOT_PRODUCT(camada2%pesos(:,j), delta_saida) * drelu(h_out(j))
      END DO
      ! Backprop para camada1
      DO i = 1, camada1%n_saida
         DO j = 1, camada1%n_entrada
            camada1%grad_pesos(i,j) = camada1%grad_pesos(i,j) + &
               delta_hidden(i) * x_ent(j) + &
               lambda_reg * camada1%pesos(i,j)
         END DO
         camada1%grad_bias(i) = camada1%grad_bias(i) + delta_hidden(i)
      END DO
      ! Atualização SGD
      camada1%pesos = camada1%pesos - eta * camada1%grad_pesos
      camada1%bias  = camada1%bias  - eta * camada1%grad_bias
      camada2%pesos = camada2%pesos - eta * camada2%grad_pesos
      camada2%bias  = camada2%bias  - eta * camada2%grad_bias
      camada1%grad_pesos = 0.0_dp; camada1%grad_bias = 0.0_dp
      camada2%grad_pesos = 0.0_dp; camada2%grad_bias = 0.0_dp
   END SUBROUTINE backprop_2camadas

   !----------------------------------------------------------------------------
   ! Normalização de features (z-score): x_norm = (x - mu) / sigma
   !----------------------------------------------------------------------------
   SUBROUTINE normalizar_zscore(X, n, p, mu, sigma, X_norm)
      INTEGER,  INTENT(IN)  :: n, p
      REAL(dp), INTENT(IN)  :: X(n, p)
      REAL(dp), INTENT(OUT) :: mu(p), sigma(p), X_norm(n, p)
      INTEGER :: j
      DO j = 1, p
         mu(j)    = SUM(X(:,j)) / REAL(n, dp)
         sigma(j) = SQRT(SUM((X(:,j) - mu(j))**2) / REAL(n-1, dp))
         IF (sigma(j) < 1.0D-10) sigma(j) = 1.0_dp
         X_norm(:,j) = (X(:,j) - mu(j)) / sigma(j)
      END DO
   END SUBROUTINE normalizar_zscore

   !----------------------------------------------------------------------------
   ! Libera memória das camadas
   !----------------------------------------------------------------------------
   SUBROUTINE liberar_camada(camada)
      TYPE(CamadaDensa), INTENT(INOUT) :: camada
      IF (ALLOCATED(camada%pesos))      DEALLOCATE(camada%pesos)
      IF (ALLOCATED(camada%bias))       DEALLOCATE(camada%bias)
      IF (ALLOCATED(camada%grad_pesos)) DEALLOCATE(camada%grad_pesos)
      IF (ALLOCATED(camada%grad_bias))  DEALLOCATE(camada%grad_bias)
   END SUBROUTINE liberar_camada

END MODULE mod_redes_neurais
