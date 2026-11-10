# GNSSVis 2.0 — atualização da simulação

Atualização do modelo de recepção e da simulação de 24 horas.

## O que está pronto

- API importável para executar e salvar a simulação, com download/cache de
  efemérides habilitável. Não há requisito de máquina offline.
- Efemérides aparentes topocêntricas do Sol por JPL DE421.
- Catálogo histórico expandido: 32 GPS com TLEs próximos à observação.
- SGP4 preservado; velocidade radial adicionada para Doppler físico.
- Espectros espalhados com perfis sinc² integrados na largura dos canais,
  combinando componentes estreitas/largas. L5 apenas nos blocos compatíveis.
- Feixe elíptico gaussiano de 50° × 40° em L1, com largura proporcional a 1/f.
  A emissão não é cortada na borda de meia potência.
- Componentes separados: ruído, satélites, Sol e potência total.
- Cinco gráficos, arrays comprimidos, passagens e metadados com proveniência.

## Cenário já executado

24 horas: 21/08/2026 às 03h até 22/08 às 03h, horário de Recife.
Banda: 1150–1600 MHz, incluindo GPS L1, L2 e L5.
Resolução: 60 segundos, 768 canais (~586 kHz por canal), float32.
Sol no eixo: 21/08 às 09h15 de Recife (12h15 UTC).
Apontamento fixo: azimute 59,819450°, elevação 52,016700°.

Foram encontradas 17 passagens de 16 satélites distintos na meia potência,
com até três GPS simultâneos no feixe em L1. Dois estão presentes no centro
do trânsito solar. Há até 15 GPS acima do horizonte; esse número não é o
número dentro do feixe.

O Doppler L1 de todos os GPS visíveis varia aproximadamente de −3,48 a
+3,77 kHz. Ele está calculado nos espectros, mas não deve ser artificialmente
ampliado para aparecer em uma imagem de 450 MHz; por isso há um painel em kHz.

O Sol quieto não foi tratado como fonte puramente síncrotron. O cenário padrão
usa um contínuo térmico ilustrativo. Existe uma componente não térmica ativa
opcional, com inclinação espectral de síncrotron opticamente fino, desativada
por padrão. Não foi inventado um evento solar na data escolhida.

Os ganhos e as amplitudes são ilustrativos: não representam medidas em dBm,
fluxo solar calibrado ou um padrão de antena medido. As coordenadas e o feixe
são configuráveis; nenhuma vinculação da GNSSVis ao BINGO é assumida.

## Arquivos entregues

- `GNSSVis_2_0_codigo_atualizado.zip`: projeto completo, fontes, referências
  históricas, novo cenário, exemplo, documentação e testes. Não inclui o
  kernel DE421 de ~17 MB; ele é baixado e armazenado na primeira execução.
- `GNSSVis_2_0_simulacao_24h.zip`: cinco PNGs, arrays NPZ, geometria,
  passagens CSV, cenário efetivo, metadados e resultado dos testes.
- `GNSSVis_2_0_waterfall_24h.png`: prévia principal em separado.

## Executar após extrair o código

Na pasta `GNSSVis-2.0`:

```bash
python -m venv .venv
source .venv/bin/activate
# No PowerShell: .venv\Scripts\Activate.ps1
pip install -e '.[dev]'
python -m gnssvis_simulation.realistic --config examples/scenarios/pibic_24h.json --output outputs/pibic_24h
python -m pytest -q
```

API:

```python
from pathlib import Path
from gnssvis_simulation import load_scenario, simulate_realistic, save_realistic

scenario = load_scenario('examples/scenarios/pibic_24h.json')
result = simulate_realistic(scenario)
save_realistic(result, Path('outputs/pibic_24h'))
```

O código pode usar internet para obter DE421 e reutiliza o arquivo depois.
Para utilizar um kernel existente, acrescente `--ephemeris caminho/de421.bsp`.
Os TLEs do exemplo já estão incluídos, com SHA-256 verificado.

Validação final: **110 testes passaram, 1 teste preexistente foi ignorado**.
Os sete controles de consistência do cenário completo também passaram.
O resultado detalhado está em `test_results.txt` no pacote da simulação.

O rótulo de desenvolvimento 2.0 foi preservado. Nenhuma versão foi publicada
no PyPI,  Consulte
`docs/REALISTIC_SIMULATION.md` para fórmulas, fontes e limitações.
