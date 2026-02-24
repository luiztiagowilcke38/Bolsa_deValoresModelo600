!==============================================================================
! Módulo 01: Séries Temporais Financeiras
! Autor: Luiz Tiago Wilcke
! Descrição: Análise avançada de séries temporais para mercado financeiro
!            brasileiro (BOVESPA/B3). Implementa decomposição, autocorrelação,
!            testes de estacionariedade e transformações.
!==============================================================================
MODULE mod_series_temporais
  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = KIND(1.0D0)

  ! Estrutura para armazenar série temporal financeira
  TYPE :: SerieFinanceira
    REAL(dp), ALLOCATABLE :: preco(:)          ! Preços de fechamento
    REAL(dp), ALLOCATABLE :: retorno(:)         ! Retornos logarítmicos
    REAL(dp), ALLOCATABLE :: retorno_acum(:)    ! Retorno acumulado
    REAL(dp), ALLOCATABLE :: tendencia(:)       ! Componente de tendência
    REAL(dp), ALLOCATABLE :: sazonalidade(:)    ! Componente sazonal
    REAL(dp), ALLOCATABLE :: residuo(:)         ! Resíduo (ruído)
    REAL(dp), ALLOCATABLE :: acf(:)             ! Autocorrelação
    REAL(dp), ALLOCATABLE :: pacf(:)            ! Autocorrelação parcial
    CHARACTER(LEN=10), ALLOCATABLE :: datas(:)  ! Datas em formato YYYY-MM-DD
    INTEGER :: num_observacoes                  ! Número de observações
    CHARACTER(LEN=20) :: ticker                 ! Código do ativo (ex: PETR4)
  END TYPE SerieFinanceira

CONTAINS

  !----------------------------------------------------------------------------
  ! Calcula retornos logarítmicos da série de preços
  ! r_t = ln(P_t / P_{t-1})
  !----------------------------------------------------------------------------
  SUBROUTINE calcular_retornos_log(serie)
    TYPE(SerieFinanceira), INTENT(INOUT) :: serie
    INTEGER :: i, n
    n = serie%num_observacoes
    IF (.NOT. ALLOCATED(serie%retorno)) ALLOCATE(serie%retorno(n-1))
    DO i = 2, n
      IF (serie%preco(i-1) > 0.0_dp) THEN
        serie%retorno(i-1) = LOG(serie%preco(i) / serie%preco(i-1))
      ELSE
        serie%retorno(i-1) = 0.0_dp
      END IF
    END DO
  END SUBROUTINE calcular_retornos_log

  !----------------------------------------------------------------------------
  ! Calcula retorno acumulado: R_T = prod(1 + r_t) - 1
  !----------------------------------------------------------------------------
  SUBROUTINE calcular_retorno_acumulado(serie)
    TYPE(SerieFinanceira), INTENT(INOUT) :: serie
    INTEGER :: i, n
    n = SIZE(serie%retorno)
    IF (.NOT. ALLOCATED(serie%retorno_acum)) ALLOCATE(serie%retorno_acum(n))
    serie%retorno_acum(1) = serie%retorno(1)
    DO i = 2, n
      serie%retorno_acum(i) = (1.0_dp + serie%retorno_acum(i-1)) * &
                               (1.0_dp + serie%retorno(i)) - 1.0_dp
    END DO
  END SUBROUTINE calcular_retorno_acumulado

  !----------------------------------------------------------------------------
  ! Decomposição aditiva da série temporal:
  ! X_t = T_t + S_t + E_t
  ! Usa média móvel centrada para extrair tendência
  !----------------------------------------------------------------------------
  SUBROUTINE decomposicao_aditiva(serie, periodo)
    TYPE(SerieFinanceira), INTENT(INOUT) :: serie
    INTEGER, INTENT(IN) :: periodo    ! Período de sazonalidade
    INTEGER :: i, j, n, meio
    REAL(dp) :: soma
    n = serie%num_observacoes
    meio = periodo / 2
    IF (.NOT. ALLOCATED(serie%tendencia))   ALLOCATE(serie%tendencia(n))
    IF (.NOT. ALLOCATED(serie%sazonalidade)) ALLOCATE(serie%sazonalidade(n))
    IF (.NOT. ALLOCATED(serie%residuo))     ALLOCATE(serie%residuo(n))
    serie%tendencia = 0.0_dp
    ! Média móvel centrada para tendência
    DO i = meio+1, n-meio
      soma = 0.0_dp
      DO j = -meio, meio
        soma = soma + serie%preco(i+j)
      END DO
      serie%tendencia(i) = soma / REAL(periodo, dp)
    END DO
    ! Extrapolamos extremos com o valor mais próximo disponível
    DO i = 1, meio
      serie%tendencia(i) = serie%tendencia(meio+1)
    END DO
    DO i = n-meio+1, n
      serie%tendencia(i) = serie%tendencia(n-meio)
    END DO
    ! Sazonalidade: média dos desvios em relação à tendência por posição
    serie%sazonalidade = 0.0_dp
    DO i = 1, n
      serie%sazonalidade(i) = serie%preco(i) - serie%tendencia(i)
    END DO
    ! Resíduo
    DO i = 1, n
      serie%residuo(i) = serie%preco(i) - serie%tendencia(i) - serie%sazonalidade(i)
    END DO
  END SUBROUTINE decomposicao_aditiva

  !----------------------------------------------------------------------------
  ! Função de Autocorrelação (ACF)
  ! rho(k) = Cov(X_t, X_{t+k}) / Var(X_t)
  !----------------------------------------------------------------------------
  SUBROUTINE calcular_acf(serie, num_lags)
    TYPE(SerieFinanceira), INTENT(INOUT) :: serie
    INTEGER, INTENT(IN) :: num_lags
    INTEGER :: i, k, n
    REAL(dp) :: media, variancia, covariancia
    n = SIZE(serie%retorno)
    IF (.NOT. ALLOCATED(serie%acf)) ALLOCATE(serie%acf(num_lags))
    media = SUM(serie%retorno) / REAL(n, dp)
    variancia = 0.0_dp
    DO i = 1, n
      variancia = variancia + (serie%retorno(i) - media)**2
    END DO
    variancia = variancia / REAL(n, dp)
    DO k = 1, num_lags
      covariancia = 0.0_dp
      DO i = 1, n-k
        covariancia = covariancia + (serie%retorno(i) - media) * &
                                    (serie%retorno(i+k) - media)
      END DO
      covariancia = covariancia / REAL(n, dp)
      IF (variancia > 1.0D-15) THEN
        serie%acf(k) = covariancia / variancia
      ELSE
        serie%acf(k) = 0.0_dp
      END IF
    END DO
  END SUBROUTINE calcular_acf

  !----------------------------------------------------------------------------
  ! Teste de Ljung-Box para ruído branco
  ! Q = n(n+2) * sum_{k=1}^{m} rho(k)^2 / (n-k)
  ! H0: série é ruído branco (retornos não-correlacionados)
  !----------------------------------------------------------------------------
  REAL(dp) FUNCTION teste_ljung_box(serie, num_lags) RESULT(estatistica_q)
    TYPE(SerieFinanceira), INTENT(IN) :: serie
    INTEGER, INTENT(IN) :: num_lags
    INTEGER :: k, n
    REAL(dp) :: soma
    n = SIZE(serie%retorno)
    soma = 0.0_dp
    DO k = 1, num_lags
      IF (k <= SIZE(serie%acf)) THEN
        soma = soma + serie%acf(k)**2 / REAL(n - k, dp)
      END IF
    END DO
    estatistica_q = REAL(n * (n + 2), dp) * soma
  END FUNCTION teste_ljung_box

  !----------------------------------------------------------------------------
  ! Teste ADF simplificado de raiz unitária (Augmented Dickey-Fuller)
  ! Delta(y_t) = alpha + beta*t + gamma*y_{t-1} + sum phi_i*Delta(y_{t-i}) + e_t
  ! H0: gamma = 0 (raiz unitária — série não-estacionária)
  !----------------------------------------------------------------------------
  REAL(dp) FUNCTION estatistica_adf(retornos, num_lags) RESULT(tau)
    REAL(dp), INTENT(IN) :: retornos(:)
    INTEGER, INTENT(IN) :: num_lags
    INTEGER :: n, i, t
    REAL(dp) :: soma_xy, soma_xx, gamma_hat, erro_padrao
    REAL(dp), ALLOCATABLE :: delta_y(:), y_lag(:)
    n = SIZE(retornos)
    ALLOCATE(delta_y(n-1), y_lag(n-1))
    DO i = 1, n-1
      delta_y(i) = retornos(i+1) - retornos(i)
      y_lag(i)   = retornos(i)
    END DO
    soma_xy = 0.0_dp; soma_xx = 0.0_dp
    DO i = 1, n-1
      soma_xy = soma_xy + y_lag(i) * delta_y(i)
      soma_xx = soma_xx + y_lag(i)**2
    END DO
    IF (ABS(soma_xx) > 1.0D-15) THEN
      gamma_hat = soma_xy / soma_xx
    ELSE
      gamma_hat = 0.0_dp
    END IF
    erro_padrao = 0.0_dp
    DO i = 1, n-1
      erro_padrao = erro_padrao + (delta_y(i) - gamma_hat * y_lag(i))**2
    END DO
    erro_padrao = SQRT(erro_padrao / REAL(n - 2, dp)) / &
                  (SQRT(soma_xx) + 1.0D-15)
    tau = gamma_hat / (erro_padrao + 1.0D-15)
    DEALLOCATE(delta_y, y_lag)
  END FUNCTION estatistica_adf

  !----------------------------------------------------------------------------
  ! Normalização Min-Max da série de preços: x' = (x - min) / (max - min)
  !----------------------------------------------------------------------------
  SUBROUTINE normalizar_minmax(valores, valores_norm)
    REAL(dp), INTENT(IN)  :: valores(:)
    REAL(dp), INTENT(OUT) :: valores_norm(:)
    REAL(dp) :: vmin, vmax, amplitude
    vmin = MINVAL(valores)
    vmax = MAXVAL(valores)
    amplitude = vmax - vmin
    IF (amplitude > 1.0D-15) THEN
      valores_norm = (valores - vmin) / amplitude
    ELSE
      valores_norm = 0.0_dp
    END IF
  END SUBROUTINE normalizar_minmax

  !----------------------------------------------------------------------------
  ! Calcula assimetria (Skewness) da distribuição de retornos
  ! S = (1/n) * sum[(x - mu)^3] / sigma^3
  !----------------------------------------------------------------------------
  REAL(dp) FUNCTION assimetria(retornos) RESULT(skew)
    REAL(dp), INTENT(IN) :: retornos(:)
    INTEGER :: n, i
    REAL(dp) :: media, desvio, soma
    n = SIZE(retornos)
    media = SUM(retornos) / REAL(n, dp)
    desvio = 0.0_dp
    DO i = 1, n
      desvio = desvio + (retornos(i) - media)**2
    END DO
    desvio = SQRT(desvio / REAL(n, dp))
    soma = 0.0_dp
    DO i = 1, n
      soma = soma + ((retornos(i) - media) / (desvio + 1.0D-15))**3
    END DO
    skew = soma / REAL(n, dp)
  END FUNCTION assimetria

  !----------------------------------------------------------------------------
  ! Calcula curtose (Kurtosis) da distribuição de retornos
  ! K = (1/n) * sum[(x - mu)^4] / sigma^4 - 3 (excesso de curtose)
  !----------------------------------------------------------------------------
  REAL(dp) FUNCTION curtose(retornos) RESULT(kurt)
    REAL(dp), INTENT(IN) :: retornos(:)
    INTEGER :: n, i
    REAL(dp) :: media, desvio, soma
    n = SIZE(retornos)
    media = SUM(retornos) / REAL(n, dp)
    desvio = 0.0_dp
    DO i = 1, n
      desvio = desvio + (retornos(i) - media)**2
    END DO
    desvio = SQRT(desvio / REAL(n, dp))
    soma = 0.0_dp
    DO i = 1, n
      soma = soma + ((retornos(i) - media) / (desvio + 1.0D-15))**4
    END DO
    kurt = soma / REAL(n, dp) - 3.0_dp     ! excesso de curtose
  END FUNCTION curtose

  !----------------------------------------------------------------------------
  ! Libera memória alocada na estrutura
  !----------------------------------------------------------------------------
  SUBROUTINE liberar_serie(serie)
    TYPE(SerieFinanceira), INTENT(INOUT) :: serie
    IF (ALLOCATED(serie%preco))       DEALLOCATE(serie%preco)
    IF (ALLOCATED(serie%retorno))     DEALLOCATE(serie%retorno)
    IF (ALLOCATED(serie%retorno_acum)) DEALLOCATE(serie%retorno_acum)
    IF (ALLOCATED(serie%tendencia))   DEALLOCATE(serie%tendencia)
    IF (ALLOCATED(serie%sazonalidade)) DEALLOCATE(serie%sazonalidade)
    IF (ALLOCATED(serie%residuo))     DEALLOCATE(serie%residuo)
    IF (ALLOCATED(serie%acf))         DEALLOCATE(serie%acf)
    IF (ALLOCATED(serie%pacf))        DEALLOCATE(serie%pacf)
    IF (ALLOCATED(serie%datas))       DEALLOCATE(serie%datas)
  END SUBROUTINE liberar_serie

END MODULE mod_series_temporais
