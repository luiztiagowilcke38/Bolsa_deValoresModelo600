!==============================================================================
! Módulo 20: Análise de Ondas de Elliott
! Autor: Luiz Tiago Wilcke
! Descrição: Detecção de padrões de ondas de Elliott com base nos níveis
!            de Fibonacci, identificação de pivôs, sequências de impulso e
!            correção (5-3), retração de Fibonacci e projeções de alvo
!            para ações da B3.
!==============================================================================
MODULE mod_ondas
   IMPLICIT NONE
   INTEGER, PARAMETER :: dp = KIND(1.0D0)
   ! Razões de Fibonacci
   REAL(dp), PARAMETER :: FIB_236 = 0.236_dp
   REAL(dp), PARAMETER :: FIB_382 = 0.382_dp
   REAL(dp), PARAMETER :: FIB_500 = 0.500_dp
   REAL(dp), PARAMETER :: FIB_618 = 0.618_dp
   REAL(dp), PARAMETER :: FIB_786 = 0.786_dp
   REAL(dp), PARAMETER :: FIB_100 = 1.000_dp
   REAL(dp), PARAMETER :: FIB_1618 = 1.618_dp
   REAL(dp), PARAMETER :: FIB_2618 = 2.618_dp
   REAL(dp), PARAMETER :: FIB_4236 = 4.236_dp

CONTAINS

   !----------------------------------------------------------------------------
   ! Identifica pivôs locais (máximos e mínimos) em uma série de preços
   ! +1 = máximo local, -1 = mínimo local, 0 = não é pivô
   !----------------------------------------------------------------------------
   SUBROUTINE detectar_pivos(precos, n, janela, pivos)
      INTEGER,  INTENT(IN)  :: n, janela
      REAL(dp), INTENT(IN)  :: precos(n)
      INTEGER,  INTENT(OUT) :: pivos(n)
      INTEGER :: i, j
      LOGICAL :: e_max, e_min
      pivos = 0
      DO i = janela+1, n-janela
         e_max = .TRUE.;  e_min = .TRUE.
         DO j = -janela, janela
            IF (j /= 0) THEN
               IF (precos(i+j) >= precos(i)) e_max = .FALSE.
               IF (precos(i+j) <= precos(i)) e_min = .FALSE.
            END IF
         END DO
         IF (e_max) pivos(i) = 1
         IF (e_min) pivos(i) = -1
      END DO
   END SUBROUTINE detectar_pivos

   !----------------------------------------------------------------------------
   ! Calcula níveis de retração de Fibonacci entre dois preços
   ! Niveis(k) = preco_max - ratio_k * (preco_max - preco_min)
   !----------------------------------------------------------------------------
   SUBROUTINE niveis_fibonacci(preco_min, preco_max, niveis)
      REAL(dp), INTENT(IN)  :: preco_min, preco_max
      REAL(dp), INTENT(OUT) :: niveis(7)
      REAL(dp) :: amplitude
      amplitude = preco_max - preco_min
      niveis(1) = preco_max - FIB_236 * amplitude   ! 23.6%
      niveis(2) = preco_max - FIB_382 * amplitude   ! 38.2%
      niveis(3) = preco_max - FIB_500 * amplitude   ! 50.0%
      niveis(4) = preco_max - FIB_618 * amplitude   ! 61.8%
      niveis(5) = preco_max - FIB_786 * amplitude   ! 78.6%
      niveis(6) = preco_min                          ! 100%
      niveis(7) = preco_min - FIB_618 * amplitude   ! 161.8% (extensão)
   END SUBROUTINE niveis_fibonacci

   !----------------------------------------------------------------------------
   ! Projeção de alvo da Onda 3 a partir de Onda 1 e Onda 2
   ! Alvo_onda3 = inicio_onda3 + 1.618 * amplitude_onda1
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION alvo_onda3(inicio_onda3, amplitude_onda1) RESULT(alvo)
      REAL(dp), INTENT(IN) :: inicio_onda3, amplitude_onda1
      alvo = inicio_onda3 + FIB_1618 * amplitude_onda1
   END FUNCTION alvo_onda3

   !----------------------------------------------------------------------------
   ! Projeção de alvo da Onda 5 (extensão da onda completa)
   ! Alvo_onda5 = inicio_onda5 + 1.000 * amplitude_onda1
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION alvo_onda5(inicio_onda5, amplitude_onda1) RESULT(alvo)
      REAL(dp), INTENT(IN) :: inicio_onda5, amplitude_onda1
      alvo = inicio_onda5 + FIB_100 * amplitude_onda1
   END FUNCTION alvo_onda5

   !----------------------------------------------------------------------------
   ! Classifica sequência de ondas (impulso = 5 ondas, correção = 3 ondas)
   ! Retorna 5 para impulso, 3 para correção ou 0 para indefinido
   !----------------------------------------------------------------------------
   INTEGER FUNCTION classificar_onda(n_pivos) RESULT(tipo)
      INTEGER, INTENT(IN) :: n_pivos
      SELECT CASE (n_pivos)
       CASE (5)
         tipo = 5    ! impulso
       CASE (3)
         tipo = 3    ! correção
       CASE DEFAULT
         tipo = 0
      END SELECT
   END FUNCTION classificar_onda

   !----------------------------------------------------------------------------
   ! Verifica regras de Elliott para ondas de impulso:
   ! Regra 1: Onda 2 não pode retrair 100% da Onda 1
   ! Regra 2: Onda 3 não pode ser a mais curta das ondas 1, 3, 5
   ! Regra 3: Onda 4 não pode entrar no território da Onda 1
   !----------------------------------------------------------------------------
   LOGICAL FUNCTION validar_impulso(V1, V2, V3, V4, V5) RESULT(valido)
      REAL(dp), INTENT(IN) :: V1, V2, V3, V4, V5
      REAL(dp) :: amp1, amp2, amp3, amp5
      amp1 = ABS(V2 - V1)
      amp2 = ABS(V3 - V2)
      amp3 = ABS(V4 - V3)
      amp5 = ABS(V5 - V4)
      valido = (amp2 < amp1 * 1.0_dp) .AND.          &  ! Onda 2 < 100% de onda 1
         (amp3 >= MIN(amp1, amp5)) .AND.         &  ! Onda 3 não é a menor
         (V4 > V1)                                  ! Onda 4 não entra na onda 1 (alta)
   END FUNCTION validar_impulso

   !----------------------------------------------------------------------------
   ! Retração percentual de um movimento entre dois preços
   !----------------------------------------------------------------------------
   REAL(dp) FUNCTION retracao_percentual(preco_inicio, preco_fim, preco_atual) &
      RESULT(ret_pct)
      REAL(dp), INTENT(IN) :: preco_inicio, preco_fim, preco_atual
      REAL(dp) :: amplitude
      amplitude = ABS(preco_fim - preco_inicio)
      IF (amplitude > 1.0D-10) THEN
         ret_pct = ABS(preco_atual - preco_fim) / amplitude
      ELSE
         ret_pct = 0.0_dp
      END IF
   END FUNCTION retracao_percentual

   !----------------------------------------------------------------------------
   ! Calcula suporte e resistência por agrupamento de pivôs (clustering simples)
   !----------------------------------------------------------------------------
   SUBROUTINE suportes_resistencias(precos, pivos, n, tol_pct, niveis_sr, n_sr)
      INTEGER,  INTENT(IN)  :: n
      REAL(dp), INTENT(IN)  :: precos(n), tol_pct
      INTEGER,  INTENT(IN)  :: pivos(n)
      REAL(dp), INTENT(OUT) :: niveis_sr(50)
      INTEGER,  INTENT(OUT) :: n_sr
      REAL(dp) :: nivel_candidato
      LOGICAL  :: novo
      INTEGER :: i, j
      n_sr = 0
      DO i = 1, n
         IF (pivos(i) /= 0) THEN
            nivel_candidato = precos(i)
            novo = .TRUE.
            DO j = 1, n_sr
               IF (ABS(niveis_sr(j) - nivel_candidato) / (niveis_sr(j) + 1.0D-10) < tol_pct) THEN
                  niveis_sr(j) = (niveis_sr(j) + nivel_candidato) / 2.0_dp
                  novo = .FALSE.
                  EXIT
               END IF
            END DO
            IF (novo .AND. n_sr < 50) THEN
               n_sr = n_sr + 1
               niveis_sr(n_sr) = nivel_candidato
            END IF
         END IF
      END DO
   END SUBROUTINE suportes_resistencias

END MODULE mod_ondas
