# Ejecuta la cadena completa de la revisión R3 en orden (Windows PowerShell).
# Uso, desde la raíz del repositorio:  ./ejecutar_todo.ps1 [-SinReplica] [-Desde 18]
# Orden: 01–16 (R) · 17 réplica SBM R2 (Python) · 18 (R) · 19 SBM R3 (Python, ~2 h) · 20 brechas (R) · 21 Tabla 2 fiel a la R2 (R) · 22 E1 niveles igualados (R)
# Requisitos: R 4.5.2 y Python (ver ENTORNO.md); datos según DATOS.md. Recomendado PowerShell 7.2+.
# 45–48 y los opcionales (30–32, 40–44) se ejecutan a mano después (ver README.md).
# EN: runs the main R3 pipeline (01–22) in order on Windows; then run 46, 47, 48 and 45 by hand (see README.md).
param([switch]$SinReplica, [string]$Desde = "01",
      [string]$Rscript = "C:\Program Files\R\R-4.5.2\bin\x64\Rscript.exe",
      [string]$Python = "python")
$ErrorActionPreference = "Stop"
Set-Location (Join-Path $PSScriptRoot "code")
New-Item -ItemType Directory -Force -Path "analisis/output", "logs" | Out-Null
$pasos = @(
  @("01","R","01_datos_tidy"), @("02","R","02_reproducir_R2"), @("03","R","03_orden_sesiones"),
  @("04","R","04_produccion_sesion"), @("05","R","05_redes_igualadas_preparacion"), @("06","R","06_psicolinguistica_sesion"),
  @("07","R","07_redes_HC_bootstrap"), @("08","R","08_auditoria_ML"), @("09","R","09_auditoria_secuencia"),
  @("10","R","10_perfiles_posicion_corregidos"), @("11","R","11_HC_primera_sesion_influencia"),
  @("12","R","12_permutacion_secuencia_sensibilidad"), @("13","R","13_auditoria_metricas_historicas"),
  @("14","R","14_preparar_entrada_SBM"), @("15","R","15_ML_agrupado"), @("16","R","16_verificacion_integrada"),
  @("17","PY","17_replica_SBM_historico"), @("18","R","18_exportar_redes_SBM"),
  @("19","PY","19_SBM_redes_limpias"), @("20","R","20_tabla2_core_periferia"),
  @("21","R","21_tabla2_niveles_R2"), @("22","R","22_tabla2_niveles_igualadas"))
foreach ($p in $pasos) {
  $n, $lang, $nombre = $p
  if ($n -lt $Desde) { continue }
  if ($n -eq "17" -and $SinReplica) { continue }
  Write-Host "== $(Get-Date -Format HH:mm:ss) $nombre"
  if ($lang -eq "R") { & $Rscript --vanilla "analisis/R/$nombre.R" *> "logs/$nombre.log" }
  else { & $Python -I "analisis/python/$nombre.py" *> "logs/$nombre.log" }
  if ($LASTEXITCODE -ne 0) { throw "FALLO en $nombre (ver logs/)" }
}
Write-Host "== Comparación con la referencia"
& $Python -I analisis/python/98_comparar_outputs.py --informe logs/comparacion_referencia.csv
if ($LASTEXITCODE -ne 0) { Write-Warning "Hay diferencias con la referencia: revisar logs/comparacion_referencia.csv" }
