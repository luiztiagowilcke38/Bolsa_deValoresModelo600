!==============================================================================
! Módulo 17: Cadeias de Markov para Regimes de Mercado
! Autor: Luiz Tiago Wilcke
! Descrição: Modelagem de regimes de mercado (bull/bear) via Cadeias de Markov,
!            modelo de Hamilton de mudança de regime (HMM simplificado),
!            estimação de probabilidades de transição e algoritmos de Viterbi
!            e Baum-Welch para identificação de estados latentes da B3.
!==============================================================================
MODULE mod_cadeias_markov
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)

CONTAINS

   !----------------------------------------------------------------------------
   ! Estimação da Matriz de Transição P por contagem de frequências
   ! P_{ij} = N_{ij} / sum_j N_{ij}
   !----------------------------------------------------------------------------
   SUBROUTINE estimar_matriz_transicao(estados, n, K, P)
      INTEGER,  INTENT(IN)  :: n, K       ! K = número de estados (ex: 2 = bull/bear)
      INTEGER,  INTENT(IN)  :: estados(n)
      REAL(dp), INTENT(OUT) :: P(K, K)
      REAL(dp) :: contagem(K, K)
      INTEGER :: t, i, j
      contagem = 0.0_dp
      DO t = 1, n - 1
         i = estados(t);  j = estados(t+1)
         IF (i >= 1 .AND. i <= K .AND. j >= 1 .AND. j <= K) &
            contagem(i, j) = contagem(i, j) + 1.0_dp
      END DO
      DO i = 1, K
         IF (SUM(contagem(i,:)) > 0.0_dp) THEN
            P(i,:) = contagem(i,:) / SUM(contagem(i,:))
         ELSE
            P(i,:) = 1.0_dp / REAL(K, dp)
         END IF
      END DO
   END SUBROUTINE estimar_matriz_transicao

   !----------------------------------------------------------------------------
   ! Distribuição estacionária (vetor pi tal que pi = pi*P)
   ! Resolve (P' - I)*pi = 0 com restrição sum(pi) = 1
   !----------------------------------------------------------------------------
   SUBROUTINE distribuicao_estacionaria(P, K, pi)
      INTEGER,  INTENT(IN)  :: K
      REAL(dp), INTENT(IN)  :: P(K, K)
      REAL(dp), INTENT(OUT) :: pi(K)
      REAL(dp) :: A(K+1, K), b(K+1), x(K)
      REAL(dp) :: ATA(K, K), ATb(K), fator
      INTEGER :: i, j, iter
      ! Monta sistema (P' - I)*pi = 0 mais sum(pi) = 1
      DO i = 1, K
         DO j = 1, K
            A(i,j) = MERGE(P(j,i) - 1.0_dp, P(j,i), i == j)
         END DO
         b(i) = 0.0_dp
      END DO
      DO j = 1, K
         A(K+1, j) = 1.0_dp
      END DO
      b(K+1) = 1.0_dp
      ! Resolve por mínimos quadrados (ATA * x = ATb)
      DO i = 1, K
         DO j = 1, K
            ATA(i,j) = DOT_PRODUCT(A(:,i), A(:,j))
         END DO
         ATb(i) = DOT_PRODUCT(A(:,i), b)
      END DO
      x = 1.0_dp / REAL(K, dp)   ! valor inicial
      DO iter = 1, 100
         DO i = 1, K
            x(i) = ATb(i)
            DO j = 1, K
               IF (j /= i) x(i) = x(i) - ATA(i,j) * x(j)
            END DO
            IF (ABS(ATA(i,i)) > 1.0D-15) x(i) = x(i) / ATA(i,i)
         END DO
      END DO
      WHERE (x < 0.0_dp) x = 0.0_dp
      pi = x / (SUM(x) + 1.0D-15)
   END SUBROUTINE distribuicao_estacionaria

   !----------------------------------------------------------------------------
   ! Probabilidade de estar no estado j em t+h dado estado i em t
   ! P^h(i,j) — elevar a matriz P à potência h
   !----------------------------------------------------------------------------
   SUBROUTINE potencia_matriz_markov(P, K, h, P_h)
      INTEGER,  INTENT(IN)  :: K, h
      REAL(dp), INTENT(IN)  :: P(K, K)
      REAL(dp), INTENT(OUT) :: P_h(K, K)
      REAL(dp) :: temp(K, K)
      INTEGER :: passo
      ! Inicializa P_h com identidade
      P_h = 0.0_dp
      DO passo = 1, K
         P_h(passo, passo) = 1.0_dp
      END DO
      DO passo = 1, h
         temp = MATMUL(P_h, P)
         P_h  = temp
      END DO
   END SUBROUTINE potencia_matriz_markov

   !----------------------------------------------------------------------------
   ! Algoritmo de Viterbi para encontrar sequência mais provável de estados
   ! Dados: observações Gaussianas X_t
   !----------------------------------------------------------------------------
   SUBROUTINE viterbi(observacoes, n, K, P, mu_estados, sigma_estados, &
      sequencia_otima)
      INTEGER,  INTENT(IN)  :: n, K
      REAL(dp), INTENT(IN)  :: observacoes(n), P(K, K)
      REAL(dp), INTENT(IN)  :: mu_estados(K), sigma_estados(K)
      INTEGER,  INTENT(OUT) :: sequencia_otima(n)
      REAL(dp) :: delta(n, K), pi_ini(K)
      INTEGER  :: psi(n, K)
      REAL(dp) :: log_emis, prob_max, prob_trial, PI_val
      INTEGER :: t, i, j, estado_max
      PI_val = 3.14159265358979323846_dp
      ! Distribuição inicial uniforme
      pi_ini = 1.0_dp / REAL(K, dp)
      ! Inicialização
      DO i = 1, K
         log_emis = -0.5_dp * LOG(2.0_dp * PI_val * sigma_estados(i)**2) &
            -0.5_dp * ((observacoes(1) - mu_estados(i))**2) / sigma_estados(i)**2
         delta(1, i) = LOG(pi_ini(i) + 1.0D-15) + log_emis
         psi(1, i)   = 0
      END DO
      ! Recursão
      DO t = 2, n
         DO j = 1, K
            prob_max = -HUGE(1.0_dp)
            stato_mx: DO i = 1, K
               prob_trial = delta(t-1, i) + LOG(MAX(P(i,j), 1.0D-15))
               IF (prob_trial > prob_max) THEN
                  prob_max  = prob_trial
                  estado_max = i
               END IF
            END DO stato_mx
            log_emis = -0.5_dp * LOG(2.0_dp * PI_val * sigma_estados(j)**2 + 1.0D-30) &
               -0.5_dp * ((observacoes(t) - mu_estados(j))**2) / &
               (sigma_estados(j)**2 + 1.0D-15)
            delta(t, j) = prob_max + log_emis
            psi(t, j)   = estado_max
         END DO
      END DO
      ! Backtracking
      prob_max = -HUGE(1.0_dp)
      sequencia_otima(n) = 1
      DO i = 1, K
         IF (delta(n, i) > prob_max) THEN
            prob_max = delta(n, i)
            sequencia_otima(n) = i
         END IF
      END DO
      DO t = n-1, 1, -1
         sequencia_otima(t) = psi(t+1, sequencia_otima(t+1))
      END DO
   END SUBROUTINE viterbi

   !----------------------------------------------------------------------------
   ! Classifica regime de mercado simples por limiar de retorno:
   ! +1 = bull (retorno > limiar), -1 = bear (retorno <= limiar)
   !----------------------------------------------------------------------------
   SUBROUTINE classificar_regime(retornos, n, limiar, regimes)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: retornos(n), limiar
      INTEGER,  INTENT(OUT) :: regimes(n)
      INTEGER :: t
      DO t = 1, n
         IF (retornos(t) > limiar) THEN
            regimes(t) = 1    ! bull
         ELSE
            regimes(t) = 2    ! bear
         END IF
      END DO
   END SUBROUTINE classificar_regime

   !----------------------------------------------------------------------------
   ! Tempo médio de permanência em cada estado (1/(1-P_ii))
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION tempo_medio_estado(P_ii) RESULT(tempo)
      REAL(dp), INTENT(IN) :: P_ii
      IF (P_ii < 1.0_dp) THEN
         tempo = 1.0_dp / (1.0_dp - P_ii)
      ELSE
         tempo = HUGE(1.0_dp)
      END IF
   END FUNCTION tempo_medio_estado

END MODULE mod_cadeias_markov
